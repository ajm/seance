import math
import functools

class SequenceError(Exception) :
    pass

class IUPAC(object) :
    codes = "ACGTUMRWSYKVHDBN"

@functools.total_ordering
class Sequence(object) :
    def __init__(self, seq, qual_str=None) :
        self.seq = seq
        self.qual = self.__generate_quals(qual_str)
        self.compressed = self.__compress_seq()
        self.ltrim = 0
        self.rtrim = len(seq)
        self.lengths = [len(seq)]

        if len(self.seq) != len(self.qual) :
            raise SequenceError("lengths of sequence and qualities are not equal")

    def __generate_quals(self, qual_str) :
        if qual_str == None :
            return len(self.seq) * [126]
        
        return map(Sequence.convert_quality, qual_str)

    def __compress_seq(self) :
        tmp = self.seq[0]
        
        for i in range(1, len(self.seq)) :
            if tmp[-1] != self.seq[i] :
                tmp += self.seq[i]
        
        return tmp

    def __eq__(self, other) :
        return self.is_duplicate(other)

    def __lt__(self, other) :
        return self.__repr__() < repr(other)

    def __len__(self) :
        return self.rtrim - self.ltrim

    def __bounded(self, s) :
        return s[self.ltrim : self.rtrim]

    def sequence(self) :
        return self.__bounded(self.seq)

    def qualities(self) :
        return self.__bounded(self.qual)

    def truncate(self, length) :
        diff = self.__len__() - length
        if diff > 0 :
            self.rtrim -= diff

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

        self.qual = map(lambda y : max(y), map(lambda *x : tuple(x), self.qual, seq2.qual))
        self.lengths += seq2.lengths

    # this obviously cannot find all duplicates, but makes the assumption
    # that sequences are sorted in descending order 
    def is_duplicate(self, seq2) :
        s1 = self.compressed # self.sequence()
        s2 = seq2.compressed # seq2.sequence()

        return s1.startswith(s2) or s2.startswith(s1)

    def __repr__(self) :
        #return repr(self.sequence())
        return repr(self.compressed)

    def __str__(self) :
        return self.sequence()

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

