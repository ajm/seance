import sys
import os
import collections

from metagenomics.system import System
from metagenomics.tools import Pagan
from metagenomics.progress import Progress


class Cluster(object) :
    def __init__(self, db, similiarity_threshold, copynumber_threshold=0) :
        self.db = db
        self.similarity_threshold = similiarity_threshold
        self.copynumber_threshold = copynumber_threshold
        self.clusters = []
        #self.profiles = []

    def to_kmers(self, seq, k) :
        s = collections.Counter()
        for i in range(len(seq)-k+1) :
            s[seq[i:i+k]] += 1
        return s

    def rough_similarity2(self, seq1, seq2, k) :
        s1 = self.to_kmers(seq1, k)
        s2 = self.to_kmers(seq2, k)

        s2.subtract(s1)

        count = 0
        for i in s2.values() :
            if i > 0 :
                count += i

        return 1.0 - (count / float(len(seq2)-k))

    def rough_similarity(self, prof, seq, k) :
        s = self.to_kmers(seq, k)

        s.subtract(prof)

        count = 0
        for i in s.values() :
            if i > 0 :
                count += i

        return 1.0 - (count / float(len(seq)-k))

    def alignment_similarity(self, seq1, seq2) :
        # write out
        f = open(os.path.join(System.tempdir(), 'tmp'), 'w')

        print >> f, seq1.fasta()
        print >> f, seq2.fasta()

        f.close()


        # align
        aligned = []
        fq = Pagan().get_454_alignment(f.name)
        fq.open()

        for seq in fq :
            if seq.id == ">consensus" :
                continue

            aligned.append(seq.sequence)

        fq.close()


        # if things are really dissimilar they do not align
        # so just give up here for this cluster
        if len(aligned) != 2 :
            return 0.0

        # calculate distance
        #same = len(filter(lambda x: x[0] == x[1], zip(aligned[0], aligned[1])))
        leng = float(min(len(aligned[0]), len(aligned[1])))

        last_gap = False
        diff = 0
        for c1,c2 in zip(aligned[0], aligned[1]) :
            gap = '-' in (c1,c2)

            if last_gap and gap :
                continue

            if c1 != c2 :
                diff += 1

            last_gap = gap

        #return same / leng
        return (leng - diff) / leng

    def create_clusters(self, keys=None) :
        seqcount = collections.Counter()

        if keys == None :
            keys = self.db

        for key in keys :
            seqcount[key] = self.db.get(key).duplicates

        p = Progress("OTU", len(seqcount))
        p.start()

        for key,freq in seqcount.most_common() :
            if freq < self.copynumber_threshold :
                break

            clustered = False

            seq = self.db.get(key)

            for cindex in range(len(self.clusters)) :
                #ksim = self.rough_similarity(self.profiles[cindex], seq, 4)
                #if ksim < 0.8 :
                #    continue

                c = self.clusters[cindex]
                similarity = self.alignment_similarity(self.db.get(c[0]), self.db.get(key))

                if similarity >= self.similarity_threshold :
                    c.append(key)
                    clustered = True
                    break

            if not clustered :
                self.clusters.append([key])
                #self.profiles.append(self.to_kmers(seq, 4))

            # the rich get richer
            self.clusters.sort(key=len, reverse=True)

            p.increment()

        p.end()

        print "Number of clusters = %d" % len(self.clusters)

