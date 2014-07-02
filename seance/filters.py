import sys
import abc
import operator
import logging

from seance.datatypes import IUPAC

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
        self.counts = []

    def add(self, f) :
        self.filters.append(f)
        self.counts.append(0)

    def accept(self, seq) :
        for index,f in enumerate(self.filters) :
            if not f.accept(seq) :
                self.counts[index] += seq.duplicates
                return False
        return True

    def reset(self) :
        for i in range(len(self.counts)) :
            self.counts[i] = 0

    def __len__(self) :
        return len(self.filters)

    def filter_counts(self) :
        return [(f.__class__.__name__ , self.counts[i]) for i,f in enumerate(self.filters)]

    def __str__(self) :
        s = ""
        for index,f in enumerate(self.filters) :
            s += "%s %d\n" % (f.__class__.__name__, self.counts[index])
        return s[:-1]

class LengthFilter(Filter) :
    def __init__(self, length) :
        if length < 0 :
            raise FilterError, "%s: length is negative (%d)" % (type(self).__name__, length)

        logging.getLogger('seance').info("created LengthFilter(length=%d)" % (length))
        self.length = length

    def accept(self, seq) :
        #if self.start != 0 :
        #    seq.remove_mid(self.start)

        if len(seq) >= self.length :
            seq.truncate(self.length)
            return True

        return False

class AmbiguousFilter(Filter) :
    def __init__(self) :
        logging.getLogger('seance').info("created AmbiguousFilter()")

    def accept(self, seq) :
        return 'N' not in seq

class MinimumQualityFilter(Filter) :
    def __init__(self, qual) :
        logging.getLogger('seance').info("created MinimumQualityFilter(qual=%d)" % (qual))
        self.qual = qual

    def accept(self, seq) :
        return min(seq.qualities) >= self.qual

class AverageQualityFilter(Filter) :
    def __init__(self, qual) :
        logging.getLogger('seance').info("created AverageQualityFilter(qual=%d)" % (qual))
        self.qual = qual

    def accept(self, seq) :
        tmp = seq.qualities
        return (sum(tmp) / float(len(tmp))) >= self.qual

class WindowedQualityFilter(Filter) :
    def __init__(self, qual, winlen) :
        logging.getLogger('seance').info("created WindowedQualityFilter(qual=%d, winlen=%d)" % (qual, winlen))
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

class HomopolymerFilter(Filter) :
    def __init__(self, maxlen) :
        logging.getLogger('seance').info("created HomopolymerFilter(maxlen=%d)" % (maxlen))
        self.maxlen = maxlen

    def accept(self, seq) :
        tmpchar = ""
        tmphp = 0

        for c in seq.sequence :
            if c != tmpchar :
                tmpchar = c
                tmphp = 0

            tmphp += 1

            if tmphp > self.maxlen :
                return False

        return True

class MidFilter(Filter) :
    def __init__(self, mid, err) :
        logging.getLogger('seance').info("created MidFilter(mid=%s, err=%d)" % (mid, err))
        self.mid = mid
        self.midlen = len(mid)
        self.err = err

    def _hamming(self, mid, seq) :
        return len(filter(lambda x: x[0] != x[1], zip(mid, seq)))

    def accept(self, seq) :
        seqmid = seq.sequence[:self.midlen]
        seq.remove_mid(self.midlen)
        return self._hamming(self.mid, seqmid) <= self.err

class PrimerFilter(Filter) :
    def __init__(self, primer, err, clip) :
        logging.getLogger('seance').info("created PrimerFilter(primer=%s, err=%d)" % (primer, err))
        self.primer = primer
        self.len = len(primer)
        self.err = err
        self.clip = clip

    def accept(self, seq) :
        seqprimer = seq.sequence[:self.len]

        ret = IUPAC.close_enough(self.primer, seq.sequence, self.err)

        if ret and self.clip :
            # primer part of sequence may be longer or
            # shorter, but it does not really matter
            # as terminal gaps are not included in our
            # definition of identity
            seq.remove_mid(self.len)
            
        return ret

