import sys
import os
import glob

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
        # first pass: ignore mid homopolymer issues
        mf = self.__build_filter(1)

        p = Progress("Reading samples", len(self.samples))
        p.start()

        for sample in self.samples :
            sample.preprocess(mf, self.options['compressed-length'], ignore_first_homopolymer=False)
            p.increment()

        p.end()


        # second pass: handle mid homopolymer sequences
        mf = self.__build_filter(2)

        p = Progress("Reading samples again", len(self.samples))
        p.start()

        for sample in self.samples :
            sample.preprocess(mf, self.options['compressed-length'], ignore_first_homopolymer=True)
            p.increment()

        p.end()


        print >> sys.stderr, "\n" + str(self.seqdb)

        self.seqdb.finalise()
        
        for sample in self.samples :
            sample.print_sample()


