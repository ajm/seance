import math
import bisect

from functools import total_ordering

class SequenceError(Exception) :
    pass

class IUPAC(object) :
    codes = "ACGTUMRWSYKVHDBN"

@total_ordering
class Sequence(object) :
    def __init__(self, seq, qual_str=None) :
        self.seq = seq
        self.qual = self.__generate_quals(qual_str)
        self.ltrim = 0
        self.rtrim = len(seq)
        self.lengths = [len(seq)]

        if len(self.seq) != len(self.qual) :
            raise SequenceError("lengths of sequence and qualities are not equal")

    def __generate_quals(self, qual_str) :
        if qual_str == None :
            return len(self.seq) * [126]
        
        return map(Sequence.convert_quality, qual_str)

    def __bounded(self, s) :
        return s[self.ltrim : self.rtrim]

    def sequence(self) :
        return self.__bounded(self.seq)

    def qualities(self) :
        return self.__bounded(self.qual)

    def truncate(self, length) :
        diff = len(self) - length
        if diff > 0 :
            self.rtrim -= diff

    def is_singular(self) :
        return len(self.lengths) == 1

    @staticmethod
    def convert_quality(q) :
        qu = ord(q) - 33
        
        if qu < 0 or qu > 126 :
            raise SequenceError("quality out of range (%d)" % qu)
        
        return qu

    def merge(self, seq2) :
        if not self.is_duplicate(seq2) :
            raise SequenceError('merging two sequences where one was not duplicate!')

        if len(seq2.seq) > len(self.seq) :
            self.seq = seq2.seq
            self.ltrim = seq2.ltrim
            self.rtrim = seq2.rtrim

        self.qual = map(lambda y : max(y), map(lambda *x : tuple(x), self.qual, seq2.qual))
        self.lengths += seq2.lengths

    def is_duplicate(self, seq2) :
        s1 = self.sequence()
        s2 = seq2.sequence()

        return s1.startswith(s2) or s2.startswith(s1)

    def __eq__(self, other) :
        return self.is_duplicate(other)

    def __lt__(self, other) :
        return (not self.is_duplicate(other)) and (self.__repr__() < repr(other))

    def __len__(self) :
        return self.rtrim - self.ltrim

    def __repr__(self) :
        return repr(self.sequence())

    def __str__(self) :
        return self.sequence()

@total_ordering
class SequenceCluster(object) :
    """ 
Each 'SequenceCluster' represents all strings that are duplicates of one another
when all runs of the same character are flated to a single instance of that character

e.g. AAABBBCCC -> ABC

All strings that are perfect prefixes of one another (without compressed) are
represented as single sequence objects. 

The goal is that the most prevalent sequence will be the final representative for each 
sequence cluster. If no clear representative sequence exists (by some margin) then we
can use PAGAN to align them together and infer the representative sequence.

    """
    def __init__(self, seq) :
        self._compressed_rep = self.__compress(seq.sequence())
        #self._sequences = [seq]
        self._sequences = SortedList([seq])

    def add(self, seq) :
        #for s in self._sequences :
        #    if s.is_duplicate(seq) :
        #        s.merge(seq)
        #        return

        #self._sequences.append(seq)
        self._sequences.insert(seq)

    def merge(self, other) :
        for s in other._sequences :
            self.add(s)

    def is_singular(self) :
        return (len(self._sequences) == 1) and self._sequences[0].is_singular()

    def __compress(self, seq) :
        tmp = seq[0]

        #for i in range(1, len(seq)) :
        #    if tmp[-1] != seq[i] :
        #        tmp += seq[i]

        for i in seq[1:] :
            if tmp[-1] != i :
                tmp += i

        return tmp

    def __lt__(self, other) :
        return repr(self) < repr(other)

    def __eq__(self, other) :
        return other._compressed_rep.startswith(self._compressed_rep) or \
               self._compressed_rep.startswith(other._compressed_rep)

    def __len__(self) :
        return len(self._sequences)

    def __repr__(self) :
        return repr(self._compressed_rep)

    def __str__(self) :
        return self._compressed_rep

class SortedList(object) :
    def __init__(self, dat=[]) :
        self._data = dat

    def insert(self, obj) :
        loc = bisect.bisect_left(self._data, obj) 

        if (len(self._data) == loc) or (not self._data[loc].is_duplicate(obj)) :
            self._data.insert(loc, obj)
        else :
            self._data[loc].merge(obj)

    def __getitem__(self, ind) :
        return self._data[ind]

    def __contains__(self, obj) :
        return self._data[bisect.bisect_left(self._data, obj)] == obj

    def __iter__(self) :
        return iter(self._data)

    def __len__(self) :
        return len(self._data)

    def __str__(self) :
        return str(self._data)

class SampleMetaData(object) :
    def __init__(self) :
        self.data = {}

    def put(self, key, value) :
        self.data[key] = value

    def get(self, key) :
        return self.data[key]

    def __str__(self) :
        s = "MetaData: "

        for k in self.data :
            s += ("%s = %s, " % (k, self.data[k]))

        return s

