import sys
import operator
import bisect
import collections

from functools import total_ordering

from metagenomics.datatypes import Sequence, SequenceError, IUPAC
from metagenomics.tools import Pagan, Aligner1D
from metagenomics.system import System
from metagenomics.progress import Progress

class SequenceDB(object) :
    def __init__(self) :
        self._db = SequenceTree()

    def put(self, seq) :
        return self._db.put(seq)

    def get(self, key) :
        return self._db.get(key)

    def debug(self) :
        self._db.debug()

    def finalise(self) :
        self._db.finalise()

    def __len__(self) :
        return len(self._db)

    def __str__(self) :
        return str(self._db)

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

    def debug(self) :
        for i in SequenceTree._db.values() :
            f = open("cluster_%d.fasta" % i.token, 'w')
            print >> f, str(i.data)
            f.close()

    def finalise(self) :
        p = Progress("Alignment", len(SequenceTree._db.values()))
        p.start()

        for sc in SequenceTree._db.values() :
            sc.data.generate_canonical_sequence()
            sc.data.write_fasta(".canon")
            p.increment()

        p.end()

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
        return "%s: added = %d, unique = %d, singulars = %d, max. cluster = %d" % \
                (type(self).__name__, self.count, SequenceTree._count, self.singulars(), self.max_cluster())

@total_ordering
class SequenceCluster(object) :
    """ 
Each 'SequenceCluster' represents all strings that are duplicates of one another
when all runs of the same character are flated to a single instance of that character

e.g. AAABBBCCC -> ABC

All strings that are perfect prefixes of one another (without compression) are
represented as single sequence objects. 

The goal is that the most prevalent sequence will be the final representative for each 
sequence cluster. If no clear representative sequence exists (by some margin) then we
can use PAGAN to align them together and infer the representative sequence.

    """
    def __init__(self, seq) :
        self._compressed_rep = self.__compress(seq.sequence())
        self._sequences = SortedList([seq])
        self._canonical_sequence = None

    def merge(self, other) :
        for seq in other._sequences :
            self._sequences.insert(seq)

    def is_singular(self) :
        return (len(self._sequences) == 1) and self._sequences[0].is_singular()

    def __compress(self, seq) :
        tmp = seq[0]

        for i in seq[1:] :
            if tmp[-1] != i :
                tmp += i

        return tmp

    def write_fasta(self, ext=".fasta") :
        fname = System.tempfilename(ext)
        f = open(fname, 'w')
        print >> f, str(self)
        f.close()

        return fname

    def generate_canonical_sequence(self) :
        # there is only one sequence anyway, 
        # so no alignment necessary
        if len(self._sequences) == 1 :
            self._canonical_sequence = self._sequences[0]
            return

        #aligned = Pagan().get_454_alignment(self.write_fasta())
        aligned = Aligner1D().get_alignment(self.write_fasta())

        chars = {}

        for seq in aligned :
            seq.trim() # without this trim canonical sequence is median length of sequences in cluster
            for i in range(len(seq)) :
                if not chars.has_key(i) :
                    chars[i] = collections.Counter()
                
                chars[i][seq[i]] += seq.duplicates

        aligned.close()

        tmp = ""
        for i in range(max(chars.keys())) :
            # tmp2 is list of tuples of the form (character, frequency)
            tmp2 = chars[i].most_common()

            max_freq = tmp2[0][1]
            bases = []

            for char,freq in tmp2 :
                if freq == max_freq :
                    bases.append(char)
                    continue
                else :
                    break

            try :
                iupac = IUPAC.get(bases)

            except SequenceError, se:
                print bases, tmp2
                print aligned.get_filename()
                raise se

            if iupac != '-' :
                tmp += iupac

        self._canonical_sequence = Sequence(tmp)

        #print str(self._canonical_sequence)
        #print aligned.get_filename()

    def __lt__(self, other) :
        return repr(self) < repr(other)

    def __is_prefix(self, seq1, seq2) :
        return seq1.startswith(seq2) or seq2.startswith(seq1)

    # comparing the compressed representative strings and ignoring the 
    # first character is in case there was a homopolymer error in the 
    # MID that is not the same as the first character in the true sequence
    def __eq__(self, other) :
        #return self.__is_prefix(self._compressed_rep, other._compressed_rep) or \
        #       self.__is_prefix(self._compressed_rep[1:], other._compressed_rep) or \
        #       self.__is_prefix(self._compressed_rep, other._compressed_rep[1:])
        return self.__is_prefix(self._compressed_rep, other._compressed_rep)

    # TODO : make this more natural
    def __len__(self) :
        return sum(map(lambda x : len(x.lengths), self._sequences))

    def __repr__(self) :
        return repr(self._compressed_rep)

    def __str__(self) :
        #return self._compressed_rep
        tmp = ""

        for i in range(len(self._sequences)) :
            count = len(self._sequences[i].lengths)

            tmp += (">seq%d NumDuplicates=%d\n" % (i, count))
            tmp += self._sequences[i].sequence()
            tmp += "\n"

        if self._canonical_sequence != None :
            tmp += (">canonical NumDuplicates=%d\n" % len(self))
            tmp += self._canonical_sequence.sequence()

        return tmp

class SortedList(object) :
    def __init__(self, dat=[]) :
        self._data = dat

    def insert(self, obj) :
        loc = bisect.bisect_left(self._data, obj) 

        if (len(self._data) == loc) or (not self._data[loc].is_duplicate(obj)) :
            self._data.insert(loc, obj)
        else :
            self._data[loc].merge(obj)

    def __getitem__(self, ind) :
        return self._data[ind]

    def __contains__(self, obj) :
        return self._data[bisect.bisect_left(self._data, obj)] == obj

    def __iter__(self) :
        return iter(self._data)

    def __len__(self) :
        return len(self._data)

    def __str__(self) :
        return str(self._data)


