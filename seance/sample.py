import sys
import collections
import os
import copy
import math
import logging

from functools import total_ordering

from seance.filetypes import FastqFile
from seance.tools import Uchime
from seance.filters import Filter
from seance.db import SequenceDB


class Sample(object) :
    def __init__(self, fastq, seqdb, filters=None, chimeras=False) :
        self.log = logging.getLogger('seance')
        self.fastq = fastq
        self.filters = filters
        self.db = seqdb

        self.seqcounts = collections.Counter()
        self.chimeras = []

        if self.filters != None :
            self.__filter_load()

            if chimeras :
                self.__detect_chimeras()
        else :
            self.__raw_load()

    def __contains__(self, seqkey) :
        return seqkey in self.seqcounts

    # True if at least one of the keys is present
    def contains(self, keys) :
        for k in keys :
            if k in self.seqcounts :
                return True
        return False

    def __filter_load(self) :
        self.fastq.open()

        for seq in self.fastq :
            if self.filters.accept(seq) :
                #print "accept", seq.id, seq.duplicates, seq.sequence[:60]
                self.seqcounts[self.db.put(seq)] += seq.duplicates
            #else :
            #    print "reject", seq.id, seq.duplicates, seq.sequence[:60]

        self.fastq.close()

        self.log.info("filter results\n" + str(self.filters))
        self.log.info("accepted %d sequences" % (sum(self.seqcounts.values())))

    def __raw_load(self) :
        self.fastq.open()

        for seq in self.fastq :
            count = seq.duplicates
            key = seq.id

            if key not in self.db :
                self.db.put(seq)
            else :
                self.db.get(key).duplicates += count

            # each sequences only appears once per sample on reload
            self.seqcounts[key] = count

        self.fastq.close()

        self.log.info("loaded %d reads (%d unique sequences)" % \
                (sum(self.seqcounts.values()), len(self.seqcounts)))

    def __detect_chimeras(self) :
        if len(self) == 0 :
            return

        self.chimeras = Uchime().run(self.print_sample(duplicate_label="/ab"))
        self.log.info("%d chimeric sequences" % len(self.chimeras))
        self.log.info("%d sequences in sample (minus chimeras)" % len(self))

    def print_sample(self, duplicate_label=" NumDuplicates", extension=".sample") :
        f = open(self.fastq.get_filename() + extension, 'w')

        for key,freq in self.seqcounts.most_common() :
            if key not in self.chimeras :
                print >> f, ">%s%s=%d" % (key, duplicate_label, freq)
                print >> f, self.db.get(key).sequence

        f.close()

        return f.name
    
    def __len__(self) :
        return sum([freq for key,freq in self.seqcounts.most_common() if key not in self.chimeras])

    def __str__(self) :
        return "%s %s (%d sequences, %d unique)" % \
               (type(self).__name__, self.fastq.get_filename(), len(self), len(self.seqcounts))

@total_ordering
class MetadataSample(Sample) :
    def __init__(self, fastq, seqdb, metadata) :
        super(MetadataSample, self).__init__(fastq, seqdb)
        self.metadata = metadata

    def description(self) :
        d = self.metadata['date']
        date = '/'.join(map(str, [d.day, d.month, d.year]))
        return ' '.join([self.metadata['id'], self.metadata['location'], self.metadata['lemur'], date])

    def __eq__(self, other) :
        return (self.metadata['location'], self.metadata['lemur'], self.metadata['date']) == \
                (other.metadata['location'], other.metadata['lemur'], other.metadata['date'])

    def __lt__(self, other) :
        return (self.metadata['location'], self.metadata['lemur'], self.metadata['date']) < \
                (other.metadata['location'], other.metadata['lemur'], other.metadata['date'])

