import sys
import abc

class Filter(object) :
    __metaclass__ = abc.ABCMeta
    
    def __init__(self) :
        pass

    @abc.abstractmethod
    def accept(self, seq) :
        return

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

class LengthFilter(Filter) :
    def __init__(self, length) :
        self.length = length

    def accept(self, seq) :
        return len(seq) > self.length

class QualityFilter(Filter) :
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

class DuplicateFilter(Filter) :
    def __init__(self) :
        self.__prev = None

    def accept(self, seq) :
        if self.__prev == None :
            self.__prev = seq
            return True

        if self.__prev.is_duplicate(seq) :
            self.__prev.merge(seq)
            return False
        else :
            self.__prev = seq
            return True

