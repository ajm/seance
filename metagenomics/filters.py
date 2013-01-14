import sys
import abc

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

class TrimLeftFilter(Filter) :
    def __init__(self, trim) :
        self.trim = trim

    def accept(self, seq) :
        seq.ltrim = self.trim
        return True

class TrimRightFilter(Filter) :
    def __init__(self, maxlen) :
        self.maxlen = maxlen

    def accept(self, seq) :
        seq.truncate(self.maxlen)
        return True

class AmbiguousFilter(Filter) :
    def accept(self, seq) :
        return 'N' not in seq

class AverageQualityFilter(Filter) :
    def __init__(self, qual) :
        self._qual = qual

    def accept(self, seq) :
        tmp = seq.qualities()
        return (sum(tmp) / float(len(tmp))) >= self._qual

class WindowedQualityFilter(Filter) :
    def __init__(self, qual, winlen) :
        self.__qual = qual
        self.__winlen = winlen

    def accept(self, seq) :
        if len(seq) < self.__winlen :
            return False

        qual = seq.qualities()
        qualsum = sum(qual[:self.__winlen])
        total = self.__qual * self.__winlen

        for i in range(len(seq) - self.__winlen) :
            if i != 0 :
                qualsum = qualsum - qual[i-1] + qual[i + self.__winlen - 1]

            if qualsum < total :
                return False
            
        return True

