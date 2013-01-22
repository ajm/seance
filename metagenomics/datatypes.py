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
        ''      : '-'
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
        self.sequence = seq
        self.qual_str = qual_str

        if self.qual_str :
            if len(self.sequence) != len(self.qual_str) :
                raise SequenceError("lengths of sequence and qualities are not equal (s=%d q=%d)" % 
                    (len(self.sequence), len(self.qual_str)))

        self.qualities = self.__generate_quals(qual_str)
        self.compressed = self.__compress(self.sequence)
        self.duplicates = 1
        self.id = None

    def __compress(self, seq) :
        tmp = seq[0]

        for i in seq[1:] :
            if tmp[-1] != i :
                tmp += i

        return tmp

    def __compress2(self, seq, length) :
        ctmp = stmp = seq[0]

        for i in seq[1:] :
            stmp += i
            if ctmp[-1] != i :
                ctmp += i
            if len(ctmp) == length :
                break

        return ctmp, stmp

    def __generate_quals(self, qual_str) :
        if qual_str is None :
            return []

        return map(Sequence.quality_to_int, qual_str)

    def truncate(self, length) :
        self.sequence = self.sequence[:length]
        self.compressed = self.__compress(self.sequence)
        
        if self.qual_str :
            self.qual_str = self.qual_str[:length]
            self.qualities = self.qualities[:length]

    def ctruncate(self, length) :
        self.compressed, self.sequence = self.__compress2(self.sequence, length)
        
        if self.qual_str :
            self.qual_str = self.qual_str[:len(self.sequence)]
            self.qualities = self.qualities[:len(self.sequence)]

    def rtrim(self, char='-') :
        self.truncate(len(self.sequence.rstrip(char)))

    def ltrim(self, char='-') :
        self.sequence.lstrip(char)
        self.compressed = self.__compress(self.sequence)

        if self.qual_str :
            tmp = len(self.qual_str) - len(self.sequence)
            self.qual_str = self.qual_str[tmp:]
            self.qualities = self.qualities[tmp:]

    def remove_mid(self) :
        self.sequence = self.sequence[10:]
        self.compressed = self.__compress(self.sequence)

        if self.qual_str :
            self.qual_str = self.qual_str[10:]
            self.qualities = self.qualities[10:]

    # XXX necessary?
    def is_singular(self) :
        return self.duplicates == 1

    @staticmethod
    def quality_to_int(q) :
        qu = ord(q) - 33
        
        if qu < 0 or qu > 126 :
            raise SequenceError("quality out of range (%d)" % qu)
        
        return qu

    @staticmethod
    def int_to_quality(i) :
        return chr(i + 33)

    def merge(self, seq2) :
        if not self.is_duplicate(seq2) :
            raise SequenceError('cannot merge two non-duplicate sequences')

        if len(seq2) > len(self) :
            self.sequence = seq2.sequence

        if self.qual_str :
            self.qualities = map(lambda y : max(y), map(lambda *x : tuple(x), self.qualities, seq2.qualities))
            self.qual_str = ''.join(map(Sequence.int_to_quality, self.qualities))
        
        self.duplicates += seq2.duplicates

    def is_duplicate(self, seq2) :
        s1 = self.sequence
        s2 = seq2.sequence

        return s1.startswith(s2) or s2.startswith(s1)

    def __iadd__(self, other) :
        self.merge(other)

    def __iter__(self) :
        return self.sequence

    def __contains__(self, ch) :
        return ch in self.sequence

    def __getitem__(self, i) :
        return self.sequence[i]

    def __eq__(self, other) :
        return self.is_duplicate(other)

    def __lt__(self, other) :
        return repr(self) < repr(other)

    def __len__(self) :
        return len(self.sequence)

    def __repr__(self) :
        return repr(self.sequence)

    def __str__(self) :
        tmp = '@' if self.qual_str else '>'
        s = "%s%s\n%s\n" % (tmp, "seq", self.sequence)
        
        if self.qual_str :
            s += ("+\n%s\n" % self.qual_str)
        
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

