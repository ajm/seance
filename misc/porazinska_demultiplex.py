from sys import stderr, exit, argv
from os import mkdir
from os.path import splitext, join, exists, isdir
from collections import Counter
from Bio import SeqIO


fname = 'SRR522901.fastq'
froot,fext = splitext(fname)

outdir = 'soil_data'
primer = 'TACAAAGGGCAGGGACGTAAT'

mid2file = {
    #'AGTATCTGGCT' : 'litter1',
    #'AGTCGACGGCT' : 'litter2',
    #'ACGTGCAGACT' : 'litter3',
    #'ATAGCGGGACT' : 'litter4',
    #'ACAGACGGACT' : 'canopy1',
    #'ATCGTGGGTAT' : 'canopy2',
    #'AGCACGGGTAT' : 'canopy3',
    #'AGTGCTGGTAT' : 'canopy4',
    'ACTAGCTGTAT' : 'soil1',
    'ATACTCTGTAT' : 'soil2',
    'ACGTGAGGGAT' : 'soil3',
    'AGACTGAGGAT' : 'soil4'
}

counts = Counter()

def close_enough(seq1, seq2, diff) :
    if diff < 0 :
        return False

    if (len(seq1) == 0) or (len(seq2) == 0) :
        return True

    return close_enough(seq1[1:], seq2[1:], diff if seq1[0] == seq2[0] else diff-1) or \
           close_enough(seq1[1:], seq2,     diff-1) or \
           close_enough(seq1,     seq2[1:], diff-1)

def contains_primer(read, primer_seq) :
    for i in range(10, 20) : # (12,17) : #start looking around the end of the 15th base (4 char tag + 11 char mid)
        if close_enough(read[i:], primer_seq, 4) :
            return True
    return False

if not exists(fname) :
    print >> stderr, "Could not find '%s'!\n" % fname
    exit(1)

# check directory
if not exists(outdir) :
    mkdir(outdir)

if exists(outdir) and not isdir(outdir) :
    print >> stderr, "'%s' exists, but is not a directory" % outdir
    exit(1)

# make output files
for i in mid2file :
    mid2file[i] = open(join(outdir, froot + "_" + mid2file[i] + '.fastq'), 'w')

#rejects = open('rejected.fastq', 'w')

# process sequences
for s in SeqIO.parse(fname, 'fastq') :
    counts['all'] += 1

    # some barcodes seem to have been reused with the other primer
    # so only use sequences matching the correct primer
    if not contains_primer(str(s.seq), primer) :
        continue

    counts['correct_primer'] += 1

    mid = str(s.seq[4:15])
    if mid in mid2file :
        print >> mid2file[mid], s.format('fastq').rstrip()
        counts[mid2file[mid].name] += 1
    else :
        tmp = []
        for i in mid2file :
            if close_enough(mid, i, 2) :
                tmp.append(i)
            
        if len(tmp) == 0 :
            #print >> rejects, s.format('fastq').rstrip()
            counts['rejected'] += 1
        elif len(tmp) == 1 :
            print >> mid2file[tmp[0]], s.format('fastq').rstrip()
            counts[mid2file[tmp[0]].name] += 1
        else :
            counts['ambiguous'] += 1

# close
for i in mid2file :
    mid2file[i].close()

#rejects.close()

# print
for i in sorted([ i.name for i in mid2file.values() ]) :
    print "%16s%8d" % (i, counts[i])


