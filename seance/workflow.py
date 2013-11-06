import sys
import os
import collections
import logging
import operator

from os.path import splitext,join,basename
from glob import glob

from seance.sample import Sample, MetadataSample
from seance.filetypes import MetadataReader, DataFileError, FastqFile, SffFile
from seance.datatypes import SampleMetadata
from seance.filters import *
from seance.db import SequenceDB
from seance.progress import Progress
from seance.tools import Sff2Fastq, GetMID, Cutadapt, Pagan, BlastN, PyroNoise
from seance.cluster import Cluster
from seance.biom import BiomFile


class WorkFlow(object) :
    def __init__(self, options) :
        self.options = options
        self.log = logging.getLogger('seance')
        self.seqdb = None

    def __filters(self, mid) :
        mf = MultiFilter()

        # the mid may be trimmed off by mothur/pyronoise 
        # or during the removal of the forward adaptor
        mid_intact = not (self.options['denoise'] or (self.options['clipprimers'] and self.options['forwardprimer'] != None))

        # if we have trimmed the forward primer, then the mid will be trimmed 
        # off with it and this check will always fail
        if mid_intact :
            mf.add(MidFilter(mid, self.options['miderrors']))

        # if --denoise is used most of these will be None unless explicitly set
        if self.options['length'] != None :
            mf.add(LengthFilter(self.options['length']))

        if self.options['removeambiguous'] :
            mf.add(AmbiguousFilter())

        if self.options['maxhomopolymer'] > 0 :
            mf.add(HomopolymerFilter(self.options['maxhomopolymer']))

        if self.options['quality'] != None :
            if self.options['windowlength'] != None :
                mf.add(
                        WindowedQualityFilter(
                            self.options['quality'], 
                            self.options['windowlength'])
                      )
            else :
                mf.add(AverageQualityFilter(self.options['quality']))

        return mf

    def __mid_sff(self, sff) :
        fastq = Sff2Fastq().run(sff, self.options['outdir'])
        mid = self.__mid_fastq(fastq)
        os.remove(fastq.get_filename())
        return mid

    def __mid_fastq(self, fastq) :
        return GetMID(self.options['midlength']).run(fastq.get_filename())

    def __get_files(self, file_names) :
        for fname in file_names :
            root,ext = splitext(fname)

            self.log.info("current file = %s" % fname)

            if ext == '.sff' :
                sff = SffFile(fname)

                # annoyingly, mid is figured out twice if we are
                # doing denoising
                if self.options['denoise'] :
                    yield PyroNoise().run(sff,
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
            
            if self.options['clipprimers'] :
                f = Cutadapt().run(f, 
                        self.options['outdir'],
                        self.options['forwardprimer'], 
                        self.options['reverseprimer'])
            
            sample = Sample(f, 
                        self.seqdb, 
                        filters=self.__filters(mid), 
                        chimeras=self.options['chimeras'])
            
            sample.print_sample()
            p.increment()

        p.end()

    def __preprocessed_samples(self) :
        mdr = MetadataReader(self.options['metadata'])
        mdr.process()

        tmp = []
        md_used = []

        for sample in glob(join(self.options['outdir'], '*.sample')) :
            md = mdr.get(basename(sample))

            if md == None :
                self.log.warn("skipping %s, metadata missing..." % basename(sample))
                continue

            tmp.append(MetadataSample(FastqFile(sample), self.seqdb, md))
            md_used.append(md['file'])

        # warn about used metadata
        for i in mdr.metadata.keys() :
            if i not in md_used :
                self.log.warn("preprocessed file for %s not found" % i)

        return tmp

    def __get_cluster_input(self, samples, duplicate_threshold, sample_threshold, contamination_threshold) :
        ref_count = collections.Counter()
        cluster_input_keys = set()

        # first collect keys for all sequences that fit number of reads
        for sample in samples :
            for key,freq in sample.seqcounts.most_common() :
                if freq < contamination_threshold :
                    continue

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
                    if freq < contamination_threshold :
                        continue

                    if self.seqdb.get(key).duplicates >= duplicate_threshold :
                        cluster_input_keys.add(key)
                        count += 1

                self.log.info("allowing singletons for (%s) - added %d sequences" % \
                        (sample.description(), count))

        return list(cluster_input_keys)

    def cluster(self) :
        # rebuild the database from preprocessed samples in outdir
        self.seqdb = SequenceDB(preprocessed=True)
        samples = self.__preprocessed_samples()
        
        # to prove the rebuild of the database is the same
        #for sample in samples :
        #    sample.print_sample(extension=".rebuild")
        
        # get keys of sequences we want to cluster
        input_keys = self.__get_cluster_input(samples,
                                self.options['duplicate-threshold'],
                                self.options['sample-threshold'],
                                self.options['contamination-threshold'])

        # clustering
        c = Cluster(self.seqdb, self.options['otu-similarity'], self.options['verbose'])
        c.create_clusters(keys=input_keys, homopolymer_correction=not self.options['no-homopolymer-correction'])

        # output centroids to file
        # output biom file
        # run blast if necessary
        centroid_fname = self.options['cluster-fasta'] #join(self.options['outdir'], 'cluster_centroids.fa')
        biom_fname = self.options['cluster-biom'] #join(self.options['outdir'], 'seance.biom')
        otu_names = {}

        # blast to get better names
        if self.options['blast-centroids'] :
            self.log.info("running blastn to get OTU names...")
            otu_names = BlastN().get_names(
                    self.__fasta(centroid_fname, c.centroids()))

        self.__fasta(centroid_fname, c.centroids(), names=otu_names)
        self.__biom(biom_fname, samples, c, otu_names)

    def __fasta(self, filename, keys, names=None) :
        f = open(filename, 'w')

        if names is None :
            for key in keys :
                print >> f, self.seqdb.get(key).fasta()
        else :
            for key in keys :
                print >> f, ">%s" % names.get(key, "%s_unknown" % key)
                print >> f, self.seqdb.get(key).sequence

        f.close()
        self.log.info("written %s" % filename)

        return filename

    def __biom(self, filename, samples, clustering, cluster_names) :
        centroids = clustering.centroids()
        all_keys = clustering.all()

        output_clusters = clustering.clusters
        output_samples = [ s for s in samples if s.contains(all_keys) ]
        output_otus = [ cluster_names.get(k, "%s_unknown" % k) for k in centroids ]

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




#-------------------------------------------------------------------------------------------------

    def phylogeny(self) :
        self.__rebuild_database()

        phy_keys = self.__get_important_keys(self.options['read-threshold'], self.options['sample-threshold'])
        phy_fasta = self.__write_fasta(phy_keys, "reference_phyla.fasta")

        # create a reference phylogeny
        print "Aligning %s sequences with PAGAN%s ..." % (len(phy_keys), "" if len(phy_keys) < 50 else ", (this might take a while)")
        ref_alignment,ref_tree = Pagan().phylogenetic_alignment(phy_fasta)

        # perform phylogenetic placement of all the reads in a sample
        # probably best to do this in the Sample class
        for sample in self.samples :
            if len(sample) == 0 :
                continue

            print "\n\n\nPhylogenetic placement: %d sequences\n\n" % len(sample)
            queries = os.path.join(self.temp_directory, sample.sff.get_basename() + ".sample")
            placement = Pagan().phylogenetic_placement(ref_alignment, ref_tree, queries)

    def otu_phylogeny(self) :
        self.__rebuild_database()

        phy_keys = self.__get_important_keys(self.options['read-threshold'], self.options['sample-threshold'])

        # cluster everything
        c = Cluster(self.seqdb, self.options['otu-similarity'])
        c.create_clusters(keys=phy_keys)

        clust_keys = map(lambda x : x[0], c.clusters)
        phy_fasta = self.__write_fasta(clust_keys, "reference_phyla.fasta")

        # blast to get better names
        print "Running blastn to get OTU names..."
        otu_names = BlastN().get_names(phy_fasta)

        phy_fasta = self.__write_fasta(clust_keys, "reference_phyla.fasta", otu_names)

        if not self.options['silva'] :
            # create a reference phylogeny
            print "Aligning %s sequences with PAGAN%s ..." % (len(clust_keys), "" if len(clust_keys) < 50 else ", (this might take a while)")
            ref_alignment,ref_tree = Pagan().phylogenetic_alignment(phy_fasta)
        else :
            print "Aligning %s sequences with PAGAN against SILVA ..." % (len(clust_keys))
            s = self.options['silva']
            ref_alignment,ref_tree = Pagan().silva_phylogenetic_alignment(s + '.fasta', s + '.tree', phy_fasta)

        # write results
        b = BiomFile()

        # if a column is going to sum to 0, then it is not worth
        # including it
        samp = []
        for sample in self.samples :
            for key in phy_keys :
                if key in sample :
                    samp.append(sample)
                    break

        print "out of %d samples %d included in heatmap" % (len(self.samples), len(samp))

#        for i in self.samples :
#            print i
#
#        print "==="
#
#        for i in samp :
#            print i

        for index, sample in enumerate(samp) :
            #b.add_sample(("%d " % index) + sample.description(), sample.metadata)
            b.add_sample(sample.description(), sample.metadata)

        for key in clust_keys :
            b.add_otu(otu_names.get(key, "%s_unknown" % key))

        if self.options['silva'] :
            fq = FastqFile(phy_fasta + '.silva.pruned.fas')
            fq.open()

            # XXX this is not very robust
            for seq in fq :
                if ';' in seq.id :
                    b.add_otu(seq.id.split(';')[-1])

            fq.close()

        # XXX PRINT EVERYTHING
        for cindex in range(len(c.clusters)) :
            cluster = c.clusters[cindex]
            key = clust_keys[cindex]
            cluster_label = otu_names.get(key, "%s_unknown" % key)

            f = open(os.path.join(self.temp_directory, "cluster_%s.fasta" % (cluster_label)), 'w')
            for read in cluster :
                print >> f, self.seqdb.get(read).fasta().rstrip()
            f.close()
        # XXX /PRINT EVERYTHING

        for sindex in range(len(samp)) :
            sample = samp[sindex]
            for cindex in range(len(c.clusters)) :
                cluster = c.clusters[cindex]
                count = 0

                key = clust_keys[cindex]
                sample_label = "%s-%s_%s" % (sample.metadata['location'], sample.metadata['lemur'], sample.metadata['file'].split('-')[0])
                sample_label = "%s" % sample.metadata['id']
                cluster_label = otu_names.get(key, "%s_unknown" % key)

                f = open(os.path.join(self.temp_directory, "sample_%s_cluster_%s.fasta" % (sample_label, cluster_label)), 'w')
                for read in cluster :
                    if read in sample :
                        count += sample.seqcounts[read]
                        #print >> f, self.seqdb.get(read).fasta().rstrip()
                        current = self.seqdb.get(read)
                        print >> f, ">%s NumDuplicates=%d\n%s" % (current.id, sample.seqcounts[read], current.sequence)

                f.close()
                b.add_quantity(cindex, sindex, count)

        b.write_to(os.path.join(self.temp_directory, "reference_phyla.biom"))

