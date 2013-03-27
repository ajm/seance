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
        self.data_directory = options['datadir']
        self.temp_directory = options['tempdir']
        self.metadata_file = options['metadata']
        self.seqdb = SequenceDB()

        self.samples = self.__create_samples()
        self.samples.sort()

    def __get_datafiles(self) :
        return glob.glob(self.data_directory + os.sep + "*sff")

    def __create_samples(self) :
        mdr = MetadataReader(self.metadata_file)
        mdr.process()

        return map(lambda x : NematodeSample(x, self.temp_directory, self.options['mid-length'], self.seqdb, mdr.get(os.path.basename(x))), self.__get_datafiles())

    def __build_filter(self) :
        mf = MultiFilter()

        if not self.options['dont-remove-nbases'] :
            mf.add(AmbiguousFilter())

        mf.add(CompressedLengthFilter(self.options['compressed-length']))

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

            sample.print_sample_raw(extension=".rebuild2")

        p.end()

    def preprocess(self) :
        mf = self.__build_filter()

        p = Progress("Reading", len(self.samples))
        p.start()

        # read in all samples and perform basic quality filtering
        for sample in self.samples :
            sample.preprocess(mf, self.options['compressed-length'], mid_errors=self.options['mid-errors'])
            p.increment()

        p.end()

        # get canonical sequences
        self.seqdb.finalise()

        p = Progress("Chimeras", len(self.samples))
        p.start()
        
        # chimera detection + clustering based on multiple alignment
        for sample in self.samples :
            sample.detect_chimeras()
            sample.print_sample_raw()

            p.increment()

        p.end()

        # see everything in the database
        self.seqdb.print_database(self.temp_directory + os.sep + "database.fasta")

        print >> sys.stderr, "\n" + str(self.seqdb)

    def __get_important_keys(self, duplicate_threshold, sample_threshold) :
        ref_count = collections.Counter()
        phy_keys = []

        # first collect keys for all sequences that fit number of reads
        for sample in self.samples :
            for key in sample.seqcounts :
                if self.seqdb.get(key).duplicates >= duplicate_threshold :
                    ref_count[key] += 1

        # then for these keys see how many samples they occurred in
        for key,freq in ref_count.most_common() :
            if freq < sample_threshold :
                break

            phy_keys.append(key)

        return phy_keys
    
    def __write_fasta(self, keys, filename, names=None) :
        f = open(os.path.join(self.temp_directory, filename), 'w')

        for key in keys :
            if names is None :
                print >> f, ">%d" % key
            else :
                print >> f, ">%s" % names.get(key, "%d_unknown" % key)

            print >> f, self.seqdb.get(key).sequence

        f.close()

        return f.name

    def phylogeny(self) :
        self.__rebuild_database()

        phy_keys = self.__get_important_keys(self.options['phyla-read-threshold'], self.options['phyla-sample-threshold'])
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

    def otu(self) :
        self.__rebuild_database()

        c = Cluster(self.seqdb, self.options['otu-similarity'], self.options['otu-dup-threshold'])
        c.create_clusters()

        f = open(os.path.join(self.temp_directory, "database.clusters"), 'w')

        for ci in range(len(c.clusters)) :
            cluster = c.clusters[ci]
            for key in cluster :
                seq = self.seqdb.get(key)
                print >> f, ">seq%d NumDuplicates=%d otu=%d" % (seq.id, seq.duplicates, ci)
                print >> f, seq.sequence

        f.close()

    def otu_phylogeny(self) :
        self.__rebuild_database()

        phy_keys = self.__get_important_keys(self.options['phyla-read-threshold'], self.options['phyla-sample-threshold'])

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

        for sample in self.samples :
            b.add_sample(sample.sample_desc(), sample.metadata)

        for key in clust_keys :
            b.add_otu(otu_names.get(key, "%d_unknown" % key))

        for sindex in range(len(self.samples)) :
            sample = self.samples[sindex]
            for cindex in range(len(c.clusters)) :
                cluster = c.clusters[cindex]
                count = 0
                for read in cluster :
                    if read in sample :
                        count += sample.seqcounts[read]

                b.add_quantity(cindex, sindex, count)

        b.write_to(os.path.join(self.temp_directory, "reference_phyla.biom"))

