import sys
import os
import glob
import collections

from metagenomics.sample import NematodeSample
from metagenomics.filetypes import MetadataReader
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

    # XXX this rebuilds something akin to the error correcting database
    #     constructed during the 'preprocess' command, but it is suboptimal
    #     for anything else
#    def __rebuild_database(self) :
#        if len(self.seqdb) != 0 :
#            return
#
#        p = Progress("Rebuild DB", len(self.samples))
#        p.start()
#
#        for sample in self.samples :
#            sample.rebuild()
#            sample.print_sample_raw(extension=".rebuild")
#            
#            p.increment()
#
#        p.end()

    def __rebuild_database(self) :
        self.seqdb = {}

        p = Progress("Rebuild DB", len(self.samples))
        p.start()

        for sample in self.samples :
            sample.rebuild(self.seqdb)
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

        #print >> sys.stderr, "\n" + str(self.seqdb)

    def phylogeny(self) :
        self.__rebuild_database()

        # we want a 'reference' phylogeny against which to do phylogenetic
        # placement of everything
        #
        # user specifies what sequences are included in the reference phylogeny
        # by two parameters, number of reads + number of samples
        duplicate_threshold = self.options['phyla-read-threshold']
        sample_threshold = self.options['phyla-sample-threshold']
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

        # TODO BLAST all of these sequences and ensure that all of them have a high identity
        # with some known species 

        # write files out + align with PAGAN
        f = open(os.path.join(self.temp_directory, "reference_phyla.fasta"), 'w')

        for key in phy_keys :
            print >> f, ">seq%d" % key
            print >> f, self.seqdb.get(key).sequence

        f.close()


        #print "\n%d / %d unique sequences\n%.2f%% of total data\n" % (len(phy_keys), len(self.seqdb), 100 * sum(map(lambda x : self.seqdb[x].duplicates, phy_keys)) / float(sum(map(lambda x : x.duplicates, self.seqdb.values()))))
        #sys.exit(-1)


        # create a reference phylogeny
        print "Aligning %s sequences with PAGAN%s ..." % (len(phy_keys), "" if len(phy_keys) < 50 else ", (this might take a while)")
        ref_alignment,ref_tree = Pagan().phylogenetic_alignment(f.name)

        # perform phylogenetic placement of all the reads in a sample
        # probably best to do this in the Sample class
        for sample in self.samples :
            if len(sample) == 0 :
                continue

            print "\n\n\nPhylogenetic placement: %d sequences\n\n" % len(sample)
            queries = os.path.join(self.temp_directory, sample.sff.get_basename() + ".sample")
            placement = Pagan().phylogenetic_placement(ref_alignment, ref_tree, queries)

        # TODO do something with the result

    def db_phylogenetic_alignment(self, read_threshold) :
        keys = []

        # find stuff to align
        for sample in self.samples :
            for key in sample.seqcounts :
                if self.seqdb.get(key).duplicates >= read_threshold :
                    keys.append(key)

        # write out
        f = open(os.path.join(self.temp_directory, "reference_phyla.fasta"), 'w')

        for key in keys :
            print >> f, ">seq%d" % key
            print >> f, self.seqdb.get(key).sequence

        f.close()

        # align
        ref_alignment,ref_tree = Pagan().phylogenetic_alignment(f.name)
        
        return ref_alignment, ref_tree

    def otu(self) :
        self.__rebuild_database()

        # pair-wise alignments
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

        # old way using MSA
#        p = Progress("OTU", len(self.samples))
#        p.start()
#
#        for sample in self.samples :
#            sample.simple_cluster(self.options['otu-similarity'])
#            p.increment()
#
#        p.end()
    def otu_phylogeny(self) :
        self.__rebuild_database()

        duplicate_threshold = self.options['phyla-read-threshold']
        sample_threshold = self.options['phyla-sample-threshold']
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

        # TODO BLAST all of these sequences and ensure that all of them have a high identity
        # with some known species 

        # cluster everything
        c = Cluster(self.seqdb, self.options['otu-similarity'])
        c.create_clusters(keys=phy_keys)

        # write files out + align with PAGAN
        f = open(os.path.join(self.temp_directory, "reference_phyla.fasta"), 'w')

        for cindex in range(len(c.clusters)) :
            key = c.clusters[cindex][0] # representative sequence
            otu_name = "OTU_%d" % cindex
            #print >> f, ">seq%d" % key
            print >> f, ">%s" % str(cindex)
            print >> f, self.seqdb.get(key).sequence

        f.close()


        print "Running blastn to get OTU names..."
        otu_names = BlastN().get_names(f.name)


        f = open(os.path.join(self.temp_directory, "reference_phyla.fasta"), 'w')

        for cindex in range(len(c.clusters)) :
            key = c.clusters[cindex][0] # representative sequence
            otu_name = "OTU_%d" % cindex
            #print >> f, ">seq%d" % key
            print >> f, ">%s" % otu_names.get(str(cindex), "%d_unknown" % cindex)
            print >> f, self.seqdb.get(key).sequence

        f.close()


        

        #print "\n%d / %d unique sequences\n%.2f%% of total data\n" % (len(phy_keys), len(self.seqdb), 100 * sum(map(lambda x : self.seqdb[x].duplicates, phy_keys)) / float(sum(map(lambda x : x.duplicates, self.seqdb.values()))))
        #sys.exit(-1)


        # create a reference phylogeny
        print "Aligning %s sequences with PAGAN%s ..." % (len(phy_keys), "" if len(phy_keys) < 50 else ", (this might take a while)")
        ref_alignment,ref_tree = Pagan().phylogenetic_alignment(f.name)


        # write results
        b = BiomFile()

        for sindex in range(len(self.samples)) :
            b.add_sample(self.samples[sindex].sample_desc())

        for cindex in range(len(c.clusters)) :
            otu_name = "OTU_%d" % cindex
            b.add_otu(otu_names.get(str(cindex), "%d_unknown" % cindex))

        for sindex in range(len(self.samples)) :
            sample = self.samples[sindex]
            for cindex in range(len(c.clusters)) :
                cluster = c.clusters[cindex]
                count = 0
                for read in cluster :
                    if read in sample :
                        count += sample.seqcounts[read]

                b.add_quantity(cindex, sindex, count)

        b.write_to(os.path.join(self.temp_directory, "otus.biom"))

