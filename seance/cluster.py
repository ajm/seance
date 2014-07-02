import sys
import os
import collections
import logging

from seance.system import System
from seance.tools import Pagan
from seance.progress import Progress
from seance.system import System


class Cluster(object) :
    def __init__(self, db, similiarity_threshold, verbose) :
        self.db = db
        self.similarity_threshold = similiarity_threshold
        self.clusters = []
        self.verbose = verbose
        self.log = logging.getLogger('seance')

    def __len__(self) :
        return len(self.clusters)

    def all(self) :
        tmp = []

        for c in self.clusters :
            tmp += c

        return tmp

    def centroids(self) :
        return [ c[0] for c in self.clusters ]
    
    def merge(self, names) :
        self.log.info("merging clusters based on labels")
        
        tmp = {} # cluster --> species

        for k,v in names.iteritems() :
            data = v.split('_')
            tmp[k] = '_'.join(data[1:-1])

        unknowns = 0
        merged_clusters = {} # species --> cluster

        # clusters are in descending order of size
        for c in self.clusters :
            try :
                species = tmp[c[0]]
            except KeyError :
                species = "unknown_%d" % unknowns
                unknowns += 1

            if species in merged_clusters :
                merged_clusters[species] += c
            else :
                merged_clusters[species] = c

        self.log.info("%d clusters merged into %d clusters based on blast hits" % (len(self.clusters), len(merged_clusters)))

        self.clusters = sorted(merged_clusters.values(), key=len, reverse=True)

    def distance(self, aligned, homopolymer_correction) :
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

        return (leng - diff) / leng

    def distance2(self, aligned, homopolymer_correction) :
        leng = float(min(len(aligned[0]), len(aligned[1])))

        # terminal homopolymer can cause problems
        #
        # pagan does 
        # XYYYY
        # X-YYY
        #
        # instead of 
        # XYYYYZ
        # XYYY-Z
        if homopolymer_correction :
            last = aligned[0][-1]

            a0 = aligned[0].rstrip(last)
            a1 = aligned[1].rstrip(last)

            if (len(a0) < len(a1)) and (a1[-1] == '-') :
                aligned[1] = a1[:-1] + (last * (len(a1) - len(a0)))

            elif (len(a1) < len(a0)) and (a0[-1] == '-') :
                aligned[0] = a0[:-1] + (last * (len(a0) - len(a1)))


        last_match = '-'
        last_gap = True 
        diff = 0

        for c1,c2 in zip(aligned[0], aligned[1]) :
            gap = '-' in (c1,c2)

            # only count gaps as a single difference
            # so continue
            if last_gap and gap :
                continue

            # however, if the gap is trailing a homopolymer
            # we may want to also ignore it
            if gap :
                # only one can have a gap
                # if last character was a match and 
                # the other current character is the same
                #   => then we are in a homopolymer
                #   => if homopolymer_correction, continue
                c = c1 if c1 != '-' else c2
                if c == last_match and homopolymer_correction :
                    continue

            if c1 != c2 :
                diff += 1

            last_match = c1 if c1 == c2 else '-'
            last_gap = gap

        if last_gap :
            diff -= 1

        return (leng - diff) / leng

    def alignment_similarity(self, seq1, seq2, homopolymer_correction) :
        # write out
        f = open(System.tempfilename(ext='cluster'), 'w')

        print >> f, seq1.fasta()
        print >> f, seq2.fasta()

        f.close()

        # align
        aligned = []
        if homopolymer_correction :
            fq = Pagan().get_454_alignment(f.name)
        else :
            fq = Pagan().get_alignment(f.name)
        
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

        return self.distance2(aligned, homopolymer_correction)

    def create_clusters(self, keys=None, homopolymer_correction=True) :
        seqcount = collections.Counter()

        if keys == None :
            keys = self.db

        for key in keys :
            seqcount[key] = self.db.get(key).duplicates

        p = Progress("Clustering", len(seqcount))
        p.start()

        for key,freq in seqcount.most_common() :
            clustered = False

            seq = self.db.get(key)

            for c in self.clusters :
                if self.alignment_similarity(self.db.get(c[0]), seq, homopolymer_correction) >= self.similarity_threshold :
                    c.append(key)
                    clustered = True
                    break

            if not clustered :
                self.clusters.append([key])

            self.clusters.sort(key=len, reverse=True)
            p.increment()

        p.end()

        #self.log.info("number of clusters = %d" % len(self.clusters))

    def create_clusters2(self, keys, homopolymer_correction=True, singletons=[], sample_threshold=1) :
        self.create_clusters(keys=keys.keys(), homopolymer_correction=homopolymer_correction)

        new_clusters = []

        for c in self.clusters :
            tmp = set()
            add_cluster = False

            for k in c :
                # short circuit in case we have any control samples
                if k in singletons :
                    add_cluster = True
                    break

                tmp.update(keys[k])

            if add_cluster or (len(tmp) >= sample_threshold) :
                new_clusters.append(c)

        self.log.info("removed %d clusters due to sample threshold" % (len(self.clusters) - len(new_clusters)))
        self.log.info("created %d clusters" % len(new_clusters))

        # i think we always want to see this
        print "removed %d clusters due to sample threshold" % (len(self.clusters) - len(new_clusters))
        print "created %d clusters" % len(new_clusters)

        self.clusters = new_clusters

