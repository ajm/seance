import sys
import os
import glob
import collections

from metagenomics.sample import NematodeSample
from metagenomics.filetypes import MetadataReader, DataFileError
from metagenomics.datatypes import SampleMetadata
from metagenomics.filters import *
from metagenomics.db import SequenceDB
from metagenomics.progress import Progress
from metagenomics.tools import Pagan, BlastN
from metagenomics.cluster import Cluster
from metagenomics.biom import BiomFile


class WorkFlow(object) :
    def __init__(self, options) :
        self.options = options
        self.data_directory = options['sffdir']
        self.temp_directory = options['outdir']
        self.metadata_file = options['metadata']
        self.seqdb = SequenceDB(compressed=options['compress'])

        self.samples = self.__create_samples()
        self.samples.sort()

    def __get_datafiles(self) :
        return glob.glob(os.path.join(self.data_directory, "*sff"))

    def __create_samples(self) :
        mdr = MetadataReader(self.metadata_file)
        mdr.process()

        tmp = []

        for sample in self.__get_datafiles() :
            fname = os.path.basename(sample)

            md = mdr.get(fname)

            if md == None :
                print >> sys.stderr, "Warning: skipping %s, metadata missing..." % fname
                continue

            s = NematodeSample(sample, self.temp_directory, self.options['mid-length'], self.seqdb, mdr.get(fname))
            tmp.append(s)

        return tmp

    def __build_filter(self) :
        mf = MultiFilter()

        if not self.options['dont-remove-nbases'] :
            mf.add(AmbiguousFilter())

        if self.options['denoise'] :
            return mf

        length = self.options['length'] + self.options['mid-length']
        if self.options['clip-primer'] :
            length += len(self.options['forward-primer'])

        if self.options['compress'] :
            mf.add(CompressedLengthFilter(length))
        else :
            mf.add(LengthFilter(length))

        if self.options['max-homopolymer'] > 0 :
            mf.add(HomopolymerFilter(self.options['max-homopolymer']))

        if self.options['minimum-quality'] != None :
            if self.options['window-length'] != None :
                mf.add(WindowedQualityFilter(self.options['minimum-quality'], self.options['window-length']))
            else :
                mf.add(AverageQualityFilter(self.options['minimum-quality']))

        return mf

    def __rebuild_database(self) :
        self.seqdb = {}

        p = Progress("Rebuild DB", len(self.samples))
        p.start()

        for sample in self.samples :
            try :
                sample.rebuild(self.seqdb)

            except DataFileError, dfe :
                print "\rError: %s\nHave all the samples been preprocessed?" % str(dfe)
                sys.exit(-1)

            sample.print_sample(extension=".rebuild")

        p.end()

    def preprocess(self) :
        mf = self.__build_filter()

        # read in all samples and perform basic quality filtering
        p = Progress("Reading", len(self.samples))
        p.start()

        for sample in self.samples :
            if self.options['denoise'] :
                sample.preprocess_denoise(mf, 
                        self.options['length'], 
                        self.options['forward-primer'], 
                        self.options['clip-primer'],
                        mid_errors=self.options['mid-errors'])

            elif self.options['compress'] :
                sample.preprocess_compress(mf, 
                        self.options['length'], 
                        mid_errors=self.options['mid-errors'])
            
            else :
                sample.preprocess(mf, 
                        self.options['length'], 
                        mid_errors=self.options['mid-errors'])
            
            p.increment()

        p.end()

        if self.options['compress'] :
            # get canonical sequences and truncate them to 'length'
            self.seqdb.finalise(self.options['length'])


        # chimera detection
        p = Progress("Chimeras", len(self.samples))
        p.start()
        
        for sample in self.samples :
            sample.detect_chimeras()
            p.increment()

        p.end()

        
        # print out everything in '.sample' files + '.database' file
        for sample in self.samples :
            sample.print_sample()
        
        self.seqdb.print_database(os.path.join(self.temp_directory, "database.fasta"))

        print >> sys.stderr, "\n" + str(self.seqdb)

    def __get_important_keys(self, duplicate_threshold, sample_threshold) :
        ref_count = collections.Counter()
        phy_keys = set()

        # first collect keys for all sequences that fit number of reads
        for sample in self.samples :
            for key in sample.seqcounts :
                if self.seqdb.get(key).duplicates >= duplicate_threshold :
                    ref_count[key] += 1

        # then for these keys see how many samples they occurred in
        for key,freq in ref_count.most_common() :
            if freq < sample_threshold :
                break

            phy_keys.add(key)

        # some samples contain reads you would not necessarily expect to see 
        # in any other samples (control samples) so add these separately, but
        # still respect duplicate_threshold
        for sample in self.samples :
            if sample.metadata['allow-singletons'] :
                for key in sample.seqcounts :
                    if self.seqdb.get(key).duplicates >= duplicate_threshold :
                        phy_keys.add(key)

        return list(phy_keys)
    
    def __write_fasta(self, keys, filename, names=None) :
        f = open(os.path.join(self.temp_directory, filename), 'w')

        if names is None :
            for key in keys :
                print >> f, self.seqdb.get(key).fasta()
        else :
            for key in keys :
                print >> f, ">%s" % names.get(key, "%s_unknown" % key)
                print >> f, self.seqdb.get(key).sequence

        f.close()

        return f.name

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

#    def otu(self) :
#        self.__rebuild_database()
#
#        c = Cluster(self.seqdb, self.options['otu-similarity'], self.options['otu-dup-threshold'])
#        c.create_clusters()
#
#        f = open(os.path.join(self.temp_directory, "database.clusters"), 'w')
#
#        for ci in range(len(c.clusters)) :
#            cluster = c.clusters[ci]
#            for key in cluster :
#                seq = self.seqdb.get(key)
#                print >> f, ">seq%d NumDuplicates=%d otu=%d" % (seq.id, seq.duplicates, ci)
#                print >> f, seq.sequence
#
#        f.close()

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

        # create a reference phylogeny
        print "Aligning %s sequences with PAGAN%s ..." % (len(phy_keys), "" if len(phy_keys) < 50 else ", (this might take a while)")
        ref_alignment,ref_tree = Pagan().phylogenetic_alignment(phy_fasta)


        # write results
        b = BiomFile()

        # if a column is going to sum to 0, then it is not worth
        # including it
        samp = []
        for sample in self.samples :
            for key in clust_keys :
                if key in sample :
                    samp.append(sample)
                    break

        for index, sample in enumerate(samp) :
            #b.add_sample(("%d " % index) + sample.description(), sample.metadata)
            b.add_sample(sample.description(), sample.metadata)

        for key in clust_keys :
            b.add_otu(otu_names.get(key, "%s_unknown" % key))

        for sindex in range(len(samp)) :
            sample = samp[sindex]
            for cindex in range(len(c.clusters)) :
                cluster = c.clusters[cindex]
                count = 0
                for read in cluster :
                    if read in sample :
                        count += sample.seqcounts[read]

                b.add_quantity(cindex, sindex, count)

        b.write_to(os.path.join(self.temp_directory, "reference_phyla.biom"))

