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

    def __build_filter(self) :
        mf = MultiFilter()

        mf.add(TrimLeftFilter(10))

        if self.options['remove-nbases'] :
            mf.add(AmbiguousFilter())

        if self.options['minlength'] != None :
            mf.add(LengthFilter(self.options['minlength']))

        if self.options['maxlength'] != None :
            mf.add(TrimRightFilter(self.options['maxlength']))

        if self.options['minquality'] != None :
            if self.options['winquality'] != None :
                mf.add(WindowedQualityFilter(self.options['minquality'], self.options['winquality']))
            else :
                mf.add(AverageQualityFilter(self.options['minquality']))

        return mf

    def run(self) :
        mf = self.__build_filter()

        p = Progress("Reading samples", len(self.samples))
        p.start()

        for sample in self.samples :
            sample.preprocess(mf)
            p.increment()

        p.end()

        print >> sys.stderr, "\n" + str(self.seqdb)

        self.seqdb.finalise()
        
        for sample in self.samples :
            sample.print_sample()

