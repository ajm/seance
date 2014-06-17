import sys
import os
import collections
import logging
import operator

from sys import exit
from os.path import splitext,join,basename,exists
from glob import glob

from seance.sample import Sample, MetadataSample
from seance.filetypes import MetadataReader, DataFileError, FastqFile, SffFile
from seance.datatypes import SampleMetadata
from seance.filters import *
from seance.db import SequenceDB
from seance.progress import Progress
from seance.tools import Sff2Fastq, GetMID2, Pagan, BlastN, AmpliconNoise
from seance.cluster import Cluster
from seance.biom import BiomFile
from seance.heatmap import heatmap as phylogenetic_heatmap
from seance.wasabi import wasabi as view_in_wasabi


class WorkFlow(object) :
    def __init__(self, options) :
        self.options = options
        self.log = logging.getLogger('seance')
        self.seqdb = None

    def __filters(self, mid) :
        mf = MultiFilter()

        ## if we have trimmed the forward primer, then the mid will be trimmed 
        ## off with it and this check will always fail
        #if not self.options['clipprimers'] :
        #    mf.add(MidFilter(mid, self.options['miderrors']))


        # always remove the mid, because it is now left by denoising
        mf.add(MidFilter(mid, self.options['miderrors']))

        # check primer, remove if requested
        if self.options['forwardprimer'] is not None :
            mf.add(PrimerFilter(self.options['forwardprimer'], 2, self.options['clipprimers']))

        if self.options['length'] != None :
            mf.add(LengthFilter(self.options['length']))

        if self.options['removeambiguous'] :
            mf.add(AmbiguousFilter())

        # homopolymer and quality are taken care of
        if self.options['denoise'] :
            return mf

        if self.options['maxhomopolymer'] > 0 :
            mf.add(HomopolymerFilter(self.options['maxhomopolymer']))

        qual = self.options['quality-method']

        if qual == 'min' :
            mf.add(MinimumQualityFilter(self.options['quality']))

        elif qual == 'average' :
            mf.add(AverageQualityFilter(self.options['quality']))

        elif qual == 'window' :
                mf.add(
                        WindowedQualityFilter(
                            self.options['quality'], 
                            elf.options['windowlength'])
                      )

        return mf

    def __mid_sff(self, sff) :
        fastq = Sff2Fastq().run(sff, self.options['outdir'])
        mid = self.__mid_fastq(fastq)
        os.remove(fastq.get_filename())
        return mid

    def __mid_fastq(self, fastq) :
        #return GetMID(self.options['midlength']).run(fastq.get_filename())
        return GetMID2(self.options['midlength']).run(fastq)

    def __get_files(self, file_names) :
        for fname in file_names :
            root,ext = splitext(basename(fname))

            self.log.info("current file = %s" % fname)

            if ext == '.sff' :
                if os.path.exists(join(self.options['outdir'], basename(fname)+'.fasta.sample')) :
                    self.log.info("skipping %s (already preprocessed)" % fname)
                    continue

                sff = SffFile(fname)

                # annoyingly, mid is figured out twice if we are
                # doing denoising
                if self.options['denoise'] :
                    yield AmpliconNoise().run(sff,
                            self.options['outdir'],
                            self.options['forwardprimer'],
                            self.__mid_sff(sff),
                            self.options['miderrors'], 
                            self.options['maxhomopolymer']
                            )
                else :
                    yield Sff2Fastq().run(sff, 
                        self.options['outdir'])

            else :
                yield FastqFile(fname)

    def preprocess(self) :
        if len(self.options['input-files']) == 0 :
            self.log.error("nothing to do")
            sys.exit(1)

        self.seqdb = SequenceDB(preprocessed=False)

        p = Progress("Preprocessing", len(self.options['input-files']), self.options['verbose'])
        p.start()

        for f in self.__get_files(self.options['input-files']) :
            mid = self.__mid_fastq(f)

            sample = Sample(f, 
                        self.options['outdir'],
                        self.seqdb, 
                        self.__filters(mid),
                        chimeras=self.options['chimeras'])
            
            sample.print_sample()
            p.increment()

        p.end()

    def __preprocessed_samples(self) :
        if self.options['metadata'] is None :
            tmp = []
            for sample in glob(join(self.options['outdir'], '*.sample')) :
                md = SampleMetadata()
                md.defaults()
                s = basename(sample)
                md['file'] = s[:s.find('.')]
                tmp.append(MetadataSample(FastqFile(sample), self.options['outdir'], self.seqdb, md))
            return tmp

        mdr = MetadataReader(self.options['metadata'])
        mdr.process()

        tmp = []
        md_used = []

        for sample in glob(join(self.options['outdir'], '*.sample')) :
            md = mdr.get(basename(sample))

            if md == None :
                self.log.warn("metadata missing for %s ..." % basename(sample))
                md = SampleMetadata()
                md.defaults()
                s = basename(sample)
                md['file'] = s[:s.find('.')]
                #self.log.warn("skipping %s, metadata missing..." % basename(sample))
                #continue

            tmp.append(MetadataSample(FastqFile(sample), self.options['outdir'], self.seqdb, md))
            md_used.append(md['file'])

        # warn about used metadata
        for i in mdr.metadata.keys() :
            if i not in md_used :
                self.log.warn("preprocessed file for %s not found" % i)

        return sorted(tmp)

    def __get_cluster_input(self, samples, duplicate_threshold, sample_threshold, contamination_threshold) :
        ref_count = collections.Counter()
        cluster_input_keys = set()

        for sample in samples :
            sample.remove_less_than(contamination_threshold)

        # first collect keys for all sequences that fit number of reads
        for sample in samples :
            for key,freq in sample.seqcounts.most_common() :
                if self.seqdb.get(key).duplicates >= duplicate_threshold :
                    ref_count[key] += 1

        # then for these keys see how many samples they occurred in
        for key,freq in ref_count.most_common() :
            if freq < sample_threshold :
                break

            cluster_input_keys.add(key)

        # some samples contain reads you would not necessarily expect to see 
        # in any other samples (control samples) so add these separately, but
        # still respect duplicate_threshold
        for sample in samples :
            if sample.metadata['allow-singletons'] :
                count = 0

                for key,freq in sample.seqcounts.most_common() :
                    if self.seqdb.get(key).duplicates >= duplicate_threshold :
                        if key not in cluster_input_keys :
                            cluster_input_keys.add(key)
                            count += 1

                self.log.info("allowing singletons for (%s) - added %d sequences" % \
                        (sample.description(), count))

        return list(cluster_input_keys)

    def __get_cluster_input2(self, samples, duplicate_threshold, contamination_threshold) :
        seq2samp = collections.defaultdict(list)
        singletons = []

        for sample in samples :
            sample.remove_less_than(contamination_threshold)

        # first collect keys for all sequences that fit number of reads
        for index, sample in enumerate(samples) :
            is_singleton = sample.metadata['allow-singletons']

            for key,freq in sample.seqcounts.most_common() :
                if is_singleton or (self.seqdb.get(key).duplicates >= duplicate_threshold) :
                    seq2samp[key].append(index)

                if is_singleton :
                    singletons.append(key)

        return seq2samp, singletons

    def cluster(self) :
        # rebuild the database from preprocessed samples in outdir
        self.seqdb = SequenceDB(preprocessed=False) # setting this to true causes things to be overwritten if we merge mulitiple preprocessing steps
        samples = self.__preprocessed_samples()
        
        if self.seqdb.num_sequences() == 0 :
            self.log.error("no sequences loaded")
            exit(1)

        # to prove the rebuild of the database is the same
        #for sample in samples :
        #    sample.print_sample(extension=".rebuild")
        
        # get keys of sequences we want to cluster
#        input_keys = self.__get_cluster_input(samples,
#                                self.options['total-duplicate-threshold'],
#                                self.options['sample-threshold'],
#                                self.options['duplicate-threshold'])

        input_keys, singleton_keys = self.__get_cluster_input2(samples,
                                    self.options['total-duplicate-threshold'],
                                    self.options['duplicate-threshold'])

        # get some info
#        num_reads = sum([ self.seqdb.get(i).duplicates for i in input_keys ])
        num_reads = sum([ self.seqdb.get(i).duplicates for i in input_keys.keys() ])
        print "clustering %d/%d (%.2f%%) sequences (%d/%d (%.2f%%) reads)" % \
                        (len(input_keys), self.seqdb.num_sequences(), \
                        len(input_keys) * 100 / float(self.seqdb.num_sequences()), \
                        num_reads, self.seqdb.num_reads(), \
                        num_reads * 100 / float(self.seqdb.num_reads()))

        # clustering
        c = Cluster(self.seqdb, self.options['otu-similarity'], self.options['verbose'])
#        c.create_clusters(keys=input_keys, homopolymer_correction=not self.options['no-homopolymer-correction'])
        c.create_clusters2(keys=input_keys, 
                           homopolymer_correction=not self.options['no-homopolymer-correction'], 
                           singletons=singleton_keys,
                           sample_threshold=self.options['sample-threshold'])

        print "created %d clusters" % len(c)

        # output centroids to file
        # output biom file
        # run blast if necessary
        centroid_fname = self.options['cluster-fasta']
        biom_fname = self.options['cluster-biom']
        otu_names = {}

        # blast to get better names
        if self.options['label-centroids'] :
            self.log.info("getting OTU names...")
            otu_names = BlastN().get_names(self.__fasta(centroid_fname, c.centroids()), self.options['label-centroids'])

            if self.options['merge-blast-hits'] :
                self.log.info("merging clusters based on labels")
                c.merge(otu_names)


        self.__fasta(centroid_fname, c.centroids(), names=otu_names)
        self.__biom(biom_fname, samples, c, otu_names)
        
        return 0

    def __fasta(self, filename, keys, names=None) :
        f = open(filename, 'w')
        
        if names is None :
            for key in keys :
                print >> f, self.seqdb.get(key).fasta()
        else :
            for key in keys :
                s = self.seqdb.get(key)
                print >> f, ">seance%s %s" % (str(key), names.get(key, "unknown"))
                print >> f, s.sequence

        f.close()
        self.log.info("written %s" % filename)

        return filename

    def __biom(self, filename, samples, clustering, cluster_names) :
        centroids = clustering.centroids()
        all_keys = clustering.all()

        output_clusters = clustering.clusters
        output_samples = [ s for s in samples if s.contains(all_keys) ]
        output_otus = [ ("seance%s" % str(k), cluster_names.get(k, "unknown")) for k in centroids ]

        self.log.info("%d / %d samples have at least one sequence used in clustering" % \
                (len(output_samples), len(samples)))

        b = BiomFile()
        b.set_samples(output_samples)
        b.set_otus(output_otus)

        for sind,sample in enumerate(output_samples) :
            for cind,cluster in enumerate(output_clusters) :
                count = 0

                for read in cluster :
                    if read in sample :
                        count += sample.seqcounts[read]

                b.add_quantity(cind, sind, count)

        b.write_to(filename)
        self.log.info("written %s" % filename)

    def phylogeny(self) :
        num_sequences = self.__count(self.options['cluster-fasta'])
        p = Pagan()

        if not self.options['silva-fasta'] :
            self.log.info("aligning %s sequences with PAGAN ..." % (num_sequences))
            alignment,tree,xmlfile = p.phylogenetic_alignment(self.options['cluster-fasta'])
        else :
            self.log.info("aligning %s sequences with PAGAN against SILVA ..." % (num_sequences))
            alignment,tree,xmlfile = p.silva_phylogenetic_alignment(self.options['silva-fasta'], 
                                                                    self.options['silva-tree'], 
                                                                    self.options['cluster-fasta'])

        os.rename(alignment, self.options['phylogeny-fasta'])
        os.rename(tree,      self.options['phylogeny-tree'])
        os.rename(xmlfile,   self.options['phylogeny-xml'])

        self.log.info("created %s" % self.options['phylogeny-fasta'])
        self.log.info("created %s" % self.options['phylogeny-tree'])
        
        if xmlfile :
            self.log.info("created %s" % self.options['phylogeny-xml'])

        return 0

    def __count(self, fasta) :
        fq = FastqFile(fasta)
        fq.open()

        count = 0
        for seq in fq :
            count += 1

        fq.close()
        return count

    def heatmap(self) :
        self.log.info("creating heatmap using %s and %s" % (self.options['cluster-biom'], self.options['phylogeny-tree']))

        if self.options['heatmap-no-tree'] :
            self.options['phylogeny-tree'] = None
        
        if self.options['phylogeny-tree'] is not None and not exists(self.options['phylogeny-tree']) :
            self.log.warn("%s does not exist, drawing heatmap without tree" % self.options['phylogeny-tree'])
            self.options['phylogeny-tree'] = None

        phylogenetic_heatmap(self.options['cluster-biom'], 
                             tree=self.options['phylogeny-tree'], 
                             output=self.options['heatmap-pdf'],
                             include=self.options['heatmap-regex'],
                             output_tree=self.options['heatmap-out-tree'],
                             flip_tree=self.options['heatmap-flip-tree'],
                             scale=self.options['heatmap-scale'],
                             tree_height_blocks=self.options['heatmap-tree-height'],
                             label_clips=self.options['heatmap-label-clip'])

        print "wrote %s" % self.options['heatmap-pdf']
        return 0

    def wasabi(self) :
        return view_in_wasabi(self.options['phylogeny-xml'], 
                              basename(self.options['outdir']), 
                              self.options['wasabi-url'],
                              self.options['wasabi-user'])

