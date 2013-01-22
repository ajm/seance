import sys
import collections
import os
import copy

from metagenomics.filetypes import SffFile, FastqFile
from metagenomics.tools import Sff2Fastq, GetMID
from metagenomics.filters import Filter
from metagenomics.db import SequenceDB

class Sample(object) :
    def __init__(self, sff_fname, workingdir, seqdb) :
        self.sff = SffFile(sff_fname)
        self.workingdir = workingdir
        self.fastq = Sff2Fastq().run(self.sff, self.workingdir)
        self.db = seqdb
        self.seqcounts = collections.Counter()

    def preprocess(self, filt, compressed_length, mid_errors=1, ignore_first_homopolymer=False) :

        mid = GetMID().run(self.fastq.get_filename())

        self.fastq.open()

        for seq in self.fastq :
            if filt.accept(seq) :
                # ensure there are not too many errors in the mid
                if self.__hamming_distance(mid, seq[:10]) > mid_errors :
                    continue

                # everything else assumes the mid does not exist
                seq.remove_mid()

                # if the last character of the mid is the same as the first character
                # of the read itself, there is a problem as the mid may end in a homopolymer
                # error, therefore we do two passes
                # during the second pass, the first homopolymer can be trimmed iff the resulting
                # sequence is already in the database, if this sequence is genuine (ie: not starting
                # with a homopolymer from the mid, then any more will be clustered with it anyway)
                if ignore_first_homopolymer :
                    #print "ifh %d %s" % (len(seq.compressed[1:compressed_length+1]), str(set(map(len, self.db._db.translate.keys()))))

                    if seq.compressed[1:compressed_length+1] in self.db :
                        #print "pre ", seq.sequence
                        seq.ltrim(seq[0])
                        #print "post", seq.sequence
                        #print ""

                # the database can only handle sequences where the compressed representation of
                # that sequence are of identical length
                seq.ctruncate(compressed_length)

                #self.seqcounts[self.db.put(seq)] += 1

                token = self.db.put(seq)
                self.seqcounts[token] += 1

        self.fastq.close()

    def print_sample(self, keys=None) :
        f = open(self.workingdir + os.sep + self.sff.get_basename() + ".sample", 'w')

        for key,freq in self.seqcounts.most_common() :
            if (keys is None) or (key in keys) :
                print >> f, ">seq%d NumDuplicates=%d" % (key, freq)
                print >> f, self.db.get(key).canonical.sequence

        f.close()

    def __hamming_distance(self, mid, seq) :
        return len(filter(lambda x: x[0] != x[1], zip(mid, seq)))

    def __merge(self, fromkey, tokey) :
        self.seqcounts[tokey] += self.seqcounts[fromkey]
        del self.seqcounts[fromkey]

    def __most_numerous_freq(self) :
        return self.seqcounts.most_common()[0][1] if len(self.seqcounts) != 0 else 0

    def __len__(self) :
        return sum(self.seqcounts.values())

    def __str__(self) :
        return "%s: %s (#seq=%d unique=%d mostfreq = %d)" % \
               (type(self).__name__, self.sff.get_basename(), len(self), len(self.seqcounts), self.__most_numerous_freq())

class NematodeSample(Sample) :
    def __init__(self, sff_fname, workingdir, seqdb, metadata) :
        super(NematodeSample, self).__init__(sff_fname, workingdir, seqdb)
        self.metadata = metadata

