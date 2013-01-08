import sys
import os
import glob

from metagenomics.sample import NematodeSample
from metagenomics.filetypes import MetaDataReader
from metagenomics.datatypes import SampleMetaData
from metagenomics.filters import LengthFilter, QualityFilter, MultiFilter
from metagenomics.db import SequenceDB

class WorkFlow(object) :
    def __init__(self, datadir, tmpdir, metadata_file) :
        self.data_directory = datadir
        self.temp_directory = tmpdir
        self.metadata_file = metadata_file
        self.seqdb = SequenceDB()

        self.samples = self.__create_samples()

    def __get_datafiles(self) :
        return glob.glob(self.data_directory + os.sep + "*sff")

    def __create_samples(self) :
        mdr = MetaDataReader(self.metadata_file)
        mdr.process()

        return map(lambda x : NematodeSample(x, self.temp_directory, self.seqdb, mdr.get(x)), self.__get_datafiles())

    def preprocess(self, length, qualitythresh, qualitywin) :
        mf = MultiFilter()
        mf.add(LengthFilter(length))
        mf.add(QualityFilter(qualitythresh, qualitywin))

        for s in self.samples :
            s.preprocess(mf)

    def run(self, length, qualitythresh, qualitywin) :
        self.preprocess(length, qualitythresh, qualitywin)

        print "\n" + str(self.seqdb)

        print "Median =", self.seqdb.median()

