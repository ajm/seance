import math

from functools import total_ordering

class SequenceError(Exception) :
    pass

class IUPAC(object) :
    codes = "ACGTUMRWSYKVHDBN"
    mapping = {
        'A'     : 'A',
        'C'     : 'C',
        'G'     : 'G',
        'T'     : 'T',
        'U'     : 'U',
        'AC'    : 'M',
        'AG'    : 'R',
        'AT'    : 'W',
        'CG'    : 'S',
        'CT'    : 'Y',
        'GT'    : 'K',
        'ACG'   : 'V',
        'ACT'   : 'H',
        'AGT'   : 'D',
        'CGT'   : 'B',
        'ACGT'  : 'N',
        ''     : '-'
    }

    @staticmethod
    def get(bases) :
        try :
            #return IUPAC.mapping[''.join(sorted(bases))]
            return IUPAC.mapping[''.join(sorted(filter(lambda x : x != '-', bases)))]

        except KeyError, ke :
            raise SequenceError("No IUPAC code for \"%s\"" % ''.join(sorted(bases)))

@total_ordering
class Sequence(object) :
    def __init__(self, seq, qual_str=None) :
        self.seq = seq
        self.qual_str = qual_str
        self.qual = self.__generate_quals(qual_str)
        self.ltrim = 0
        self.rtrim = len(seq)
        self.lengths = [len(seq)]

        if len(self.seq) != len(self.qual) :
            raise SequenceError("lengths of sequence and qualities are not equal (s=%d q=%d)" % (len(self.seq), len(self.qual)))

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

    def __contains__(self, ch) :
        return ch in self.sequence()

    def __getitem__(self, i) :
        return self.sequence()[i]

    def __eq__(self, other) :
        return self.is_duplicate(other)

    def __lt__(self, other) :
        return (not self.is_duplicate(other)) and (self.__repr__() < repr(other))

    def __len__(self) :
        return self.rtrim - self.ltrim

    def __repr__(self) :
        return repr(self.sequence())

    def __str__(self) :
        #return "@%s\n%s\n+\n%s" % ("seq", self.sequence(), self.__bounded(self.qual_str))
        tmp = '@' if self.qual_str else '>'
        s = "%s%s\n%s\n" % (tmp, "seq", self.sequence())
        if self.qual_str :
            s += ("+\n%s\n" % (self.__bounded(self.qual_str)))
        return s


class SampleMetadata(object) :
    def __init__(self) :
        self.data = {}

    def put(self, key, value) :
        self.data[key] = value

    def get(self, key) :
        return self.data[key]

    def __str__(self) :
        s = "Metadata: "

        for k in self.data :
            s += ("%s = %s, " % (k, self.data[k]))

        return s

