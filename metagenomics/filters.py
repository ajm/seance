import sys
import abc
import operator

class FilterError(Exception) :
    pass

class Filter(object) :
    __metaclass__ = abc.ABCMeta
    
    def __init__(self) :
        pass

    @abc.abstractmethod
    def accept(self, seq) :
        pass

class NullFilter(object) :
    def accept(self) :
        return True

class MultiFilter(Filter) :
    def __init__(self) :
        self.filters = []

    def add(self, f) :
        self.filters.append(f)

    def accept(self, seq) :
        for f in self.filters :
            if not f.accept(seq) :
                return False
        return True

    def __len__(self) :
        return len(self.filters)

class LengthFilter(Filter) :
    def __init__(self, length) :
        if length < 0 :
            raise FilterError, "%s: length is negative (%d)" % (type(self).__name__, length)

        self.length = length

    def accept(self, seq) :
        return len(seq) > self.length

class CompressedLengthFilter(LengthFilter) :
    def __init__(self, length) :
        super(CompressedLengthFilter, self).__init__(length+11) 
        # XXX +11 for worst-case scenario: no runs of characters in the mid tag, 
        # but last character in the mid is wrongly extended
        # XXX this might be nonsense, esp. as mid is not necessarily 10...

    def accept(self, seq) :
        if len(seq.compressed) < self.length :
            return False

        #seq.ctruncate(self.length)
        return True

class AmbiguousFilter(Filter) :
    def accept(self, seq) :
        return 'N' not in seq

class AverageQualityFilter(Filter) :
    def __init__(self, qual) :
        self.qual = qual

    def accept(self, seq) :
        tmp = seq.qualities
        return (sum(tmp) / float(len(tmp))) >= self.qual

class WindowedQualityFilter(Filter) :
    def __init__(self, qual, winlen) :
        self.qual = qual
        self.winlen = winlen

    def accept(self, seq) :
        if len(seq) < self.winlen :
            return False

        qual = seq.qualities
        qualsum = sum(qual[:self.winlen])
        total = self.qual * self.winlen

        for i in range(len(seq) - self.winlen) :
            if i != 0 :
                qualsum = qualsum - qual[i-1] + qual[i + self.winlen - 1]

            if qualsum < total :
                return False

        return True

class MIDHomopolymer(Filter) :
    def __init__(self, reject=True) :
        self.op = operator.ne if reject else operator.eq

    def accept(self, seq) :
        return self.op(seq[9], seq[10])

