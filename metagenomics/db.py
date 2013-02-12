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
        self._db = SequenceDict()

    def put(self, seq) :
        return self._db.put(seq)

    def get(self, key) :
        return self._db.get(key)

    def debug(self) :
        self._db.debug()

    def finalise(self) :
        self._db.finalise()

    def print_database(self, fname) :
        self._db.print_database(fname)

    def __contains__(self, obj) :
        return obj in self._db

    def __len__(self) :
        return len(self._db)

    def __str__(self) :
        return str(self._db)

class SequenceDict(object):
    def __init__(self) :
        self.db = {}
        self.translate = {}
        self.count = 0

    def generate_key(self, s) :
        self.translate[s] = len(self.db)
        return self.translate[s]

    def put(self, seq) :
        c = SequenceCluster(seq)
        tkey = repr(c)
        key = None

        #if self.db.has_key(key) :
        #    self.db[key].merge(c)
        #else :
        #    self.db[key] = c

        if self.translate.has_key(tkey) :
            key = self.translate[tkey]
            self.db[key].merge(c)
        else :
            key = self.generate_key(tkey)
            self.db[key] = c

        self.count += 1

        return key

    def get(self, key) :
        return self.db[key]

    def finalise(self) :
        p = Progress("Alignment", len(self.db))
        p.start()

        for sc in self.db.values() :
            sc.generate_canonical_sequence()
            sc.write_fasta(".canon")
            p.increment()

        p.end()

    def singulars(self) :
        tmp = 0

        for k in self.db :
            if self.db[k].is_singular() :
                tmp += 1

        return tmp

    def max_cluster(self) :
        tmp = -1

        for k in self.db :
            if len(self.db[k]) > tmp :
                tmp = len(self.db[k])

        return tmp

    def print_database(self, fname) :
        f = open(fname, 'w')

        for seqclust_key in sorted(self.db, reverse=True, key=lambda x : len(self.db[x])) :
            print >> f, ">seq%d NumDuplicates=%d" % (seqclust_key, len(self.db[seqclust_key]))
            print >> f, self.db[seqclust_key].canonical.sequence

        f.close()

    def has_sequence(self, cseq) :
        return self.translate.has_key(cseq)

    def __contains__(self, cseq) :
        return self.has_sequence(cseq)

    def __len__(self) :
        return len(self.db)

    def __str__(self) :
        return "%s: added = %d, unique = %d, singulars = %d, max. cluster = %d" % \
                (type(self).__name__, self.count, len(self), self.singulars(), self.max_cluster())

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
        self.compressed = seq.compressed
        self.sequences = SortedList([seq])
        self.canonical = None

    def merge(self, other) :
        for seq in other.sequences :
            self.sequences.insert(seq)

    def is_singular(self) :
        return (len(self.sequences) == 1) and self.sequences[0].is_singular()

    def write_fasta(self, ext=".fasta") :
        fname = System.tempfilename(ext)
        f = open(fname, 'w')
        print >> f, str(self)
        f.close()

        return fname

    def generate_canonical_sequence(self) :
        # there is only one sequence anyway, 
        # so no alignment necessary
        if len(self.sequences) == 1 :
            self.canonical = self.sequences[0]
            return

        # 1. align sequences
        #aligned = Pagan().get_454_alignment(self.write_fasta())
        aligned = Aligner1D().get_alignment(self.write_fasta())
        aligned.open()

        chars = {}

        # 2. parse aligned file
        for seq in aligned :
            # pagan adds the consensus sequence to the end, but
            # we will do it ourselves using a weighted majority
            # vote and allow IUPAC codes
            if seq.id == ">consensus" :
                continue

            # trim trailing gaps
            # or else canonical sequence will be median length
            # of sequences in cluster
            seq.rtrim('-')

            for i in range(len(seq)) :
                if not chars.has_key(i) :
                    chars[i] = collections.Counter()
                
                chars[i][seq[i]] += seq.duplicates

        aligned.close()

        # 3. construct canonical sequence
        tmp = ""
        for i in range(max(chars.keys())+1) :
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

            # most frequent character might be ambiguous
            # so give IUPAC code instead (not possible with Aligner1D)
            try :
                iupac = IUPAC.get(bases)

            except SequenceError, se:
                print bases, tmp2
                print aligned.get_filename()
                raise se

            if iupac != '-' :
                tmp += iupac

        self.canonical = Sequence(tmp)

    def __iadd__(self, other) :
        self.merge(other)

    def __lt__(self, other) :
        return repr(self) < repr(other)

    def __eq__(self, other) :
        return self.compressed == other.compressed

    def __len__(self) :
        return sum(map(lambda x : x.duplicates, self.sequences))

    def __repr__(self) :
        return self.compressed

    def __str__(self) :
        tmp = ""

        for i in range(len(self.sequences)) :
            tmp += (">seq%d NumDuplicates=%d\n" % (i, self.sequences[i].duplicates))
            tmp += self.sequences[i].sequence
            tmp += "\n"

        if self.canonical != None :
            tmp += (">canonical NumDuplicates=%d\n" % len(self))
            tmp += self.canonical.sequence

        return tmp

class SortedList(object) :
    def __init__(self, data=[]) :
        self.data = data

    def insert(self, obj) :
        loc = bisect.bisect_left(self.data, obj) 

        if (len(self.data) == loc) or (not self.data[loc].is_duplicate(obj)) :
            self.data.insert(loc, obj)
        else :
            self.data[loc].merge(obj)

    def __getitem__(self, ind) :
        return self.data[ind]

    def __contains__(self, obj) :
        return self.data[bisect.bisect_left(self.data, obj)] == obj

    def __iter__(self) :
        return iter(self.data)

    def __len__(self) :
        return len(self.data)

    def __str__(self) :
        return str(self.data)

