import sys
import collections
import os
import copy
import math

from functools import total_ordering

from metagenomics.filetypes import SffFile, FastqFile
from metagenomics.tools import Sff2Fastq, GetMID, Pagan, Uchime, PyroNoise
from metagenomics.filters import Filter
from metagenomics.db import SequenceDB

class Sample(object) :
    def __init__(self, sff_fname, workingdir, mid_length, seqdb) :
        self.sff = SffFile(sff_fname)
        self.workingdir = workingdir
        self.fastq = None #Sff2Fastq().run(self.sff, self.workingdir)
        self.mid_length = mid_length
        self.db = seqdb
        self.seqcounts = collections.Counter()
        self.chimeras = []

    def __contains__(self, seqkey) :
        return seqkey in self.seqcounts

    def preprocess(self, filt, length, mid_errors=0) :
        self.fastq = Sff2Fastq().run(self.sff, self.workingdir)
        mid = GetMID(self.mid_length).run(self.fastq.get_filename())

        self.fastq.open()

        for seq in self.fastq :
            if filt.accept(seq) :
                if self.__hamming_distance(mid, seq[:self.mid_length]) > mid_errors :
                    continue

                seq.remove_mid(self.mid_length)
                seq.truncate(length)
                self.seqcounts[self.db.put(seq)] += 1

        self.fastq.close()
    
    def preprocess_compress(self, filt, compressed_length, mid_errors=0) :
        self.fastq = Sff2Fastq().run(self.sff, self.workingdir)
        mid = GetMID(self.mid_length).run(self.fastq.get_filename())

        self.fastq.open()

        for seq in self.fastq :
            if filt.accept(seq) :
                # ensure there are not too many errors in the mid
                if self.__hamming_distance(mid, seq[:self.mid_length]) > mid_errors :
                    continue

                # everything else assumes the mid does not exist
                seq.remove_mid(self.mid_length)

                # the database can only handle sequences where the compressed representation of
                # that sequence are of identical length
                seq.ctruncate(compressed_length)

                self.seqcounts[self.db.put(seq)] += 1

        self.fastq.close()
    
    def preprocess_denoise(self, filt, length, primer, mid_errors=0) :
        self.fastq = Sff2Fastq().run(self.sff, self.workingdir)
        mid = GetMID(self.mid_length).run(self.fastq.get_filename())

        self.fastq = PyroNoise().run(self.sff.get_filename(), primer, mid)
        self.fastq.open()

        for seq in self.fastq :
            if filt.accept(seq) :
                self.seqcounts[self.db.put(seq)] += 1

        self.fastq.close()

    def rebuild(self, newdb) :
        self.db = newdb
        self.seqcounts.clear()

        self.fastq = FastqFile(os.path.join(self.workingdir, self.sff.get_basename() + ".sample"))
        self.fastq.open()

        for seq in self.fastq :
            count = seq.duplicates
            key = seq.id

            if key not in self.db :
                self.db[key] = seq
            else :
                self.db[key].duplicates += count

            self.seqcounts[key] = count

        self.fastq.close()

    def detect_chimeras(self) :
        if len(self) == 0 :
            return
        
        self.chimeras = Uchime().run(self.print_sample(duplicate_label="/ab"))

    def print_sample(self, whitelist=None, duplicate_label=" NumDuplicates", extension=".sample") :
        f = open(os.path.join(self.workingdir, self.sff.get_basename() + extension), 'w')

        for key,freq in self.seqcounts.most_common() :
            if ((whitelist is None) or (key in whitelist)) and (key not in self.chimeras) :
                #print >> sys.stderr, key, freq
                print >> f, ">%s%s=%d" % (key, duplicate_label, freq)
                print >> f, self.db.get(key).sequence

        f.close()

        return f.name

    def __hamming_distance(self, mid, seq) :
        return len(filter(lambda x: x[0] != x[1], zip(mid, seq)))

    def __average_length(self) :
        total_length = 0
        num_seq = 0

        for key,freq in self.seqcounts.iteritems() :
            total_length += (freq * len(self.db.get(key).sequence))
            num_seq += freq

        return total_length / num_seq if num_seq != 0 else 0

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

@total_ordering
class NematodeSample(Sample) :
    def __init__(self, sff_fname, workingdir, mid_length, seqdb, metadata) :
        super(NematodeSample, self).__init__(sff_fname, workingdir, mid_length, seqdb)
        self.metadata = metadata

    def sample_desc(self) :
        d = self.metadata.get('date')
        date = '/'.join(map(str, [d.day, d.month, d.year]))
        return ' '.join([self.metadata.get('location'), self.metadata.get('lemur'), date])

    def __eq__(self, other) :
        return (self.metadata.get('location'), self.metadata.get('lemur'), self.metadata.get('date')) == \
                (other.metadata.get('location'), other.metadata.get('lemur'), other.metadata.get('date'))

    def __lt__(self, other) :
        return (self.metadata.get('location'), self.metadata.get('lemur'), self.metadata.get('date')) < \
                (other.metadata.get('location'), other.metadata.get('lemur'), other.metadata.get('date'))

