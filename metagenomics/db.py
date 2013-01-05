import sys
import bisect

from metagenomics.datatypes import Sequence, SequenceCluster

class SequenceDB(object) :
    def __init__(self) :
        self._db = SequenceTree()

    def put(self, seq) :
        return self._db.put(seq)

    def get(self, key) :
        return self._db.get(key)

    def __len__(self) :
        return len(self._db)

    def __str__(self) :
        return str(self._db)

class SequenceList(object) :
    def __init__(self) :
        self._db = []
        self._count = 0

    def put(self, seq) :
        self._count += 1
        try :
            tmp = self._db.index(seq)
            self._db[tmp].merge(seq)
            return tmp

        except ValueError, ve :
            self._db.append(seq)
            return len(self._db) - 1

    def get(self, key) :
        return self._db[key]

    def singulars(self) :
        return map(lambda x : len(x.lengths), self._db).count(1)

    def __len__(self) :
        return len(self._db)

    def __str__(self) :
        return "%s: added = %d, unique = %d, singulars = %d" % (type(self).__name__, self._count, len(self._db), self.singulars())

class SequenceTree(object) :
    _count = 0
    _db = {}

    def __init__(self) :
        self.left = None
        self.right = None
        self.data = None

        self.count = 0
        self.token = SequenceTree._count
        SequenceTree._db[self.token] = self
        SequenceTree._count += 1

    def put(self, seq) :
        return self.add_cluster(SequenceCluster(seq))

    def get(self, key) :
        return SequenceTree._db[key].data

    def add_cluster(self, clust) :
        self.count += 1

        if self.data is None :
            self.data = clust
            return self.token

        elif self.data == clust :
            self.data.merge(clust)
            return self.token

        elif self.data < clust :
            if self.left is None :
                self.left = SequenceTree()
            
            return self.left.add_cluster(clust)

        else :
            if self.right is None :
                self.right = SequenceTree()
            
            return self.right.add_cluster(clust)

    def singulars(self) :
        tmp = 0
        if self.data and self.data.is_singular() :
            tmp += 1
            
        if self.left :
            tmp += self.left.singulars()
            
        if self.right :
            tmp += self.right.singulars()

        return tmp

    def max_cluster(self) :
        tmp = [0]

        if self.left :
            tmp.append(self.left.max_cluster())

        if self.right :
            tmp.append(self.right.max_cluster())

        if self.data :
            tmp.append(len(self.data))

        return max(tmp)

    def __len__(self) :
        return SequenceTree._count

    def __str__(self) :
        return "%s: added = %d, unique = %d, singulars = %d, max. cluster = %d" % (type(self).__name__, self.count, SequenceTree._count, self.singulars(), self.max_cluster())

