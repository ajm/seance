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

        return map(lambda x : NematodeSample(x, self.temp_directory, self.options['mid-length'], self.seqdb, mdr.get(x)), self.__get_datafiles())

    def __build_filter(self) :
        mf = MultiFilter()

        if self.options['remove-nbases'] :
            mf.add(AmbiguousFilter())

        mf.add(CompressedLengthFilter(self.options['compressed-length']))

        if self.options['minimum-quality'] != None :
            if self.options['window-length'] != None :
                mf.add(WindowedQualityFilter(self.options['minimum-quality'], self.options['window-length']))
            else :
                mf.add(AverageQualityFilter(self.options['minimum-quality']))

        return mf

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
            if len(sample) != 0 :
                sample.detect_chimeras()
                sample.simple_cluster(0.97)
            else :
                sample.print_sample_raw()

            p.increment()

        p.end()

        # see everything in the database
        self.seqdb.print_database(self.temp_directory + os.sep + "database.fasta")

        print >> sys.stderr, "\n" + str(self.seqdb)

    def phylogeny(self) :
        if len(self.seqdb) == 0 :
            p = Progress("Reconstruct", len(self.samples))
            p.start()
            
            for sample in self.samples :
                sample.reconstruct()
                p.increment()
            
            p.end()


