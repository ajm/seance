import sys
import operator
import bisect
import collections
import os

from functools import total_ordering

from seance.datatypes import Sequence, SequenceError, IUPAC
from seance.tools import Pagan, Aligner1D
from seance.system import System
from seance.progress import Progress


class SequenceDB(object) :
    def __init__(self, preprocessed=False) :
        if preprocessed :
            self._db = WrapperDict()
        else :
            self._db = SequenceDict()

    def put(self, seq) :
        return self._db.put(seq)

    def get(self, key) :
        return self._db.get(key)

    def num_reads(self) :
        return self._db.num_reads()

    def num_sequences(self) :
        return self._db.num_sequences()

    def print_database(self, fname) :
        self._db.print_database(fname)

    def __contains__(self, obj) :
        return obj in self._db

    def __len__(self) :
        return len(self._db)

    def __str__(self) :
        return str(self._db)

class SequenceDict(object) :
    def __init__(self) :
        self.db = {}
        self.translate = {}
        self.count = 0

    def generate_key(self, s) :
        self.translate[s] = len(self.db)
        return self.translate[s]

    def put(self, seq) :
        tkey = repr(seq)
        key = None

        if self.translate.has_key(tkey) :
            key = self.translate[tkey]
            self.db[key].merge(seq)
        else :
            key = self.generate_key(tkey)
            seq.id = str(key) # overwrite the original name
            self.db[key] = seq

        self.count += 1

        return key

    def get(self, key) :
        return self.db[key]

    def finalise(self, length) :
        pass

    def print_database(self, fname) :
        f = open(fname, 'w')

        for seqclust_key in sorted(self.db, reverse=True, key=lambda x : len(self.db[x])) :
            print >> f, ">seq%d NumDuplicates=%d" % (seqclust_key, self.db[seqclust_key].duplicates)
            print >> f, self.db[seqclust_key].sequence

        f.close()

    def has_sequence(self, cseq) :
        return self.translate.has_key(cseq)

    def __contains__(self, cseq) :
        return self.has_sequence(cseq)

    def num_reads(self) :
        return sum([ i.duplicates for i in self.db.values() ])

    def num_sequences(self) :
        return len(self)

    def __len__(self) :
        return len(self.db)

    def __str__(self) :
        return "%s: added = %d, unique = %d" % \
                (type(self).__name__, self.count, len(self))

class WrapperDict(dict) :
    def __init__(self) :
        super(WrapperDict, self).__init__()

    def put(self, seq) :
        key = seq.id
        count = seq.duplicates

        if self.has_key(key) :
            self.__getitem__(key).duplicates += count
        else :
            self.__setitem__(key, seq)

        return key

    def num_reads(self) :
        return sum([ i.duplicates for i in self.values() ])

    def num_sequences(self) :
        return len(self)

