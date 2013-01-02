import sys
import bisect

# TODO actually write an efficient version of this
#      at the moment I just want the interface to 
#      be clean and don't care about speed
class SequenceDB(object) :
    def __init__(self) :
        self._db = []
        self._count = 0

    # put a sequence and return a key
    def put(self, seq) :
        self._count += 1
        try :
            tmp = self._db.index(seq)
            self._db[tmp].merge(seq)
            return tmp

        except ValueError, ve :
            self._db.append(seq)
            return len(self._db) - 1

    # return the sequence object associated with the sequence key
    def get(self, key) :
        return self._db[key]

    def singulars(self) :
        return map(lambda x : len(x.lengths), self._db).count(1)

    def __len__(self) :
        return len(self._db)

    def __str__(self) :
        return "%d (total added = %d, singulars = %d)" % (len(self._db), self._count, self.singulars())

# this is junk
# i thought a sorted list would do, but then the insert
# time is O(n)
class SortedList(object) :
    def __init__(self) :
        self._data = []

    def insert(self, obj) :
        loc = bisect.bisect_left(self._data, obj) 

        if (len(self._data) == loc) or (self._data[loc] != obj) :
            self._data.insert(loc, obj)
        else :
            self._data[loc].merge(obj)

        return self._data[loc].sequence_number(), self._data[loc]

    def __contains__(self, obj) :
        return self._data[bisect.bisect_left(self._data, obj)] == obj

