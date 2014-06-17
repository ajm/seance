import math
import datetime

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

    reverse_mapping = {
        'A' : 'A',
        'C' : 'C',
        'G' : 'G',
        'T' : 'T',
        'U' : 'U',
        'M' : 'AC',
        'R' : 'AG',
        'W' : 'AT',
        'S' : 'CG',
        'Y' : 'CT',
        'K' : 'GT',
        'V' : 'ACG',
        'H' : 'ACT',
        'D' : 'AGT',
        'B' : 'CGT',
        'N' : 'ACGT',
        '-' : ''
    }

    @staticmethod
    def close_enough(primer, sequence, diff) :
        if diff < 0 :
            return False

        if (len(primer) == 0) or (len(sequence) == 0) :
            return True

        m = IUPAC.equal(sequence[0], primer[0])

        return IUPAC.close_enough(primer[1:], sequence[1:], diff if m else diff-1) or \
               IUPAC.close_enough(primer[1:], sequence,     diff-1) or \
               IUPAC.close_enough(primer,     sequence[1:], diff-1)

    @staticmethod
    def seq_position(primer, sequence, diff) :
        for i in range(len(sequence) - len(primer)) :
            if IUPAC.close_enough(primer, sequence[i:], diff) :
                return i
        return -1
    
    @staticmethod
    def seq_position_reverse(primer, sequence, diff) :
        for i in range(len(sequence) - len(primer))[::-1] :
            if IUPAC.close_enough(primer, sequence[i:], diff) :
                return i
        return -1

#    @staticmethod
#    def equal(base, iupac) :
#        if   base == 'A' :
#            return iupac in "AMRWVHDN"
#        elif base == 'G' :
#            return iupac in "GRSKVDBN"
#        elif base == 'T' :
#            return iupac in "TWYKHDBN"
#        elif base == 'C' :
#            return iupac in "CMSYVHBN"

    # now both can contain iupac characters
    @staticmethod
    def equal(a, b) :
        a2 = IUPAC.reverse_mapping[a]
        b2 = IUPAC.reverse_mapping[b]

        for i in a2 :
            if i in b2 :
                return True
        
        return False

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
        self.duplicates = 1
        self.id = None

    def __generate_quals(self, qual_str) :
        if qual_str is None :
            return []

        return map(Sequence.quality_to_int, qual_str)

    def truncate(self, length) :
        self.sequence = self.sequence[:length]

        if self.qual_str :
            self.qual_str = self.qual_str[:length]
            self.qualities = self.qualities[:length]

    def rtrim(self, char='-') :
        self.truncate(len(self.sequence.rstrip(char)))

    def ltrim(self, char='-') :
        self.sequence = self.sequence.lstrip(char)

        if self.qual_str :
            tmp = len(self.qual_str) - len(self.sequence)
            self.qual_str = self.qual_str[tmp:]
            self.qualities = self.qualities[tmp:]

    def remove_mid(self, mid_length) :
        self.sequence = self.sequence[mid_length:]

        if self.qual_str :
            self.qual_str = self.qual_str[mid_length:]
            self.qualities = self.qualities[mid_length:]

    def ungap(self) :
        self.sequence = self.sequence.replace('-','')

    def back_translate(self) :
        self.sequence = self.sequence.replace('U','T')

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

    def fasta(self, duplabel=' NumDuplicates') :
        return ">%s%s=%d\n%s" % (self.id, duplabel, self.duplicates, self.sequence)

    def fastq(self, duplabel=' NumDuplicates') :
        if not self.qual_str :
            raise Exception('cannot print as fastq, no qualities present')

        return "@%s\n+\n%s" % (self.fasta(duplabel)[1:], self.qual_str)

    def __str__(self) :
        raise Exception

        # TODO sensible default

        tmp = '@' if self.qual_str else '>'
        s = "%s%s\n%s\n" % (tmp, self.id if self.id != None else "seq", self.sequence)
        
        if self.qual_str :
            s += ("+\n%s\n" % self.qual_str)
        
        return s

class SampleMetadata(object) :
    def __init__(self) :
        self.data = {}

    def defaults(self) :
        self.data = {
                'id'                : '',
                'file'              : '',
                'lemur'             : '',
                'date'              : datetime.date(1970, 1, 1),
                'location'          : '',
                'eggs'              : -1,
                'allow-singletons'  : False
                }

    def __setitem__(self, key, value) :
        self.data[key] = value

    def __getitem__(self, key) :
        return self.data[key]

    def items(self) :
        return self.data.items()

    def __str__(self) :
        s = "Metadata: "

        for k in self.data :
            s += ("%s = %s, " % (k, self.data[k]))

        return s

