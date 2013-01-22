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

        return map(lambda x : NematodeSample(x, self.temp_directory, self.seqdb, mdr.get(x)), self.__get_datafiles())

    def __build_filter(self, phase) :
        mf = MultiFilter()

        mf.add(MIDHomopolymer(True if phase == 1 else False))
        
        if self.options['remove-nbases'] :
            mf.add(AmbiguousFilter())

        mf.add(CompressedLengthFilter(self.options['compressed-length']))

        if self.options['minimum-quality'] != None :
            if self.options['window-length'] != None :
                mf.add(WindowedQualityFilter(self.options['minimum-quality'], self.options['window-length']))
            else :
                mf.add(AverageQualityFilter(self.options['minimum-quality']))

        return mf

    def run(self) :
        # two passes over all the samples
        # first : ignore reads where a homopolymer from the mid appears to bleed into the rest of the read
        # seconds : handle ignored reads
        for phase,ignore_first in [(1, False), (2, True)] :
            mf = self.__build_filter(phase)

            p = Progress("Reading samples (pass %d)" % phase, len(self.samples))
            p.start()

            for sample in self.samples :
                sample.preprocess(mf, 
                        self.options['compressed-length'], 
                        mid_errors=self.options['mid-errors'], 
                        ignore_first_homopolymer=ignore_first)
                p.increment()

            p.end()

        # get canonical sequences
        self.seqdb.finalise()

        # get ids from each sample
#        ids = collections.Counter()
#        for sample in self.samples :
#            for key in sample.seqcounts.keys() :
#                ids[key] += 1
#        
#        tmp = []
#        for key,count in ids.items() :
#            if count > 1 :
#                tmp.append(key)
#
        # print out samples
        for sample in self.samples :
            sample.print_sample()
#            sample.print_sample(tmp)


        print >> sys.stderr, "\n" + str(self.seqdb)
