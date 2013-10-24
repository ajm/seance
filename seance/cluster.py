import sys
import os
import collections
import logging

from seance.system import System
from seance.tools import Pagan
from seance.progress import Progress


class Cluster(object) :
    def __init__(self, db, similiarity_threshold, verbose) :
        self.db = db
        self.similarity_threshold = similiarity_threshold
        self.clusters = []
        self.verbose = verbose
        self.log = logging.getLogger('seance')

    def all(self) :
        tmp = []

        for c in self.clusters :
            tmp += c

        return tmp

    def centroids(self) :
        return [ c[0] for c in self.clusters ]

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

        # delete tmp files
        os.remove(f.name)
        os.remove(fq.get_filename())

        # if things are really dissimilar they do not align
        # so just give up here for this cluster
        if len(aligned) != 2 :
            return 0.0

        # calculate distance
        #same = len(filter(lambda x: x[0] == x[1], zip(aligned[0], aligned[1])))
        leng = float(min(len(aligned[0]), len(aligned[1])))

        last_gap = True 
        diff = 0
        for c1,c2 in zip(aligned[0], aligned[1]) :
            gap = '-' in (c1,c2)

            if last_gap and gap :
                continue

            if c1 != c2 :
                diff += 1

            last_gap = gap

        if last_gap :
            diff -= 1

        #return same / leng
        return (leng - diff) / leng

    def create_clusters(self, keys=None) :
        seqcount = collections.Counter()

        if keys == None :
            keys = self.db

        for key in keys :
            seqcount[key] = self.db.get(key).duplicates

        p = Progress("Clustering", len(seqcount), self.verbose)
        p.start()

        for key,freq in seqcount.most_common() :
            clustered = False

            seq = self.db.get(key)

            for c in self.clusters :
                if self.alignment_similarity(self.db.get(c[0]), seq) >= self.similarity_threshold :
                    c.append(key)
                    clustered = True
                    break

            if not clustered :
                self.clusters.append([key])

            self.clusters.sort(key=len, reverse=True)
            p.increment()

        p.end()

        self.log.info("number of clusters = %d" % len(self.clusters))

