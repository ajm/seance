import sys
import collections

from metagenomics.filetypes import SffFile, FastqFile
from metagenomics.tools import Sff2Fastq
from metagenomics.filters import Filter
from metagenomics.db import SequenceDB

class Sample(object) :
    def __init__(self, sff_fname, workingdir, seqdb) :
        self._sff = SffFile(sff_fname)
        self._workingdir = workingdir
        self._fastq = None
        self._db = seqdb
        self._seqcounts = collections.Counter()

    def process(self, filt) :
        self._fastq = Sff2Fastq().run(self._sff, self._workingdir)

        for seq in self._fastq.get_generator() :
            if filt.accept(seq) :
                self._seqcounts[self._db.put(seq)] += 1

        tmp = self._seqcounts.most_common()[0][1] if len(self._seqcounts) != 0 else 0
        print "%s:\tsample = %d\tdb = %s\t(most freq. = %d)" % \
              (self._sff.get_basename(), len(self), str(self._db), tmp)

    def merge(self, fromkey, tokey) :
        self._seqcounts[tokey] += self._seqcounts[fromkey]
        del self._seqcounts[fromkey]

    def __len__(self) :
        return sum(self._seqcounts.values())

    def __str__(self) :
        return "%s %s: len = %d" % \
               (self.__class__.__name__, self._sff.get_basename(), self.__len__())

class NematodeSample(Sample) :
    def __init__(self, sff_fname, workingdir, seqdb, metadata) :
        super(NematodeSample, self).__init__(sff_fname, workingdir, seqdb)
        self._metadata = metadata

