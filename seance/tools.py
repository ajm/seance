import sys
import abc
import os
import commands
import re
import urllib2
import logging
import glob
import shutil
import collections

from os.path import abspath, join, dirname
from seance.filetypes import SffFile, FastqFile
from seance.datatypes import IUPAC

class ExternalProgramError(Exception) :
    pass

class ExternalProgramNotInstalledError(Exception) :
    pass

class ExternalProgram(object) :
    __metaclass__ = abc.ABCMeta

    def __init__(self, pname) :
        self.programname = pname
        self.log = logging.getLogger('seance')
    
    @staticmethod
    def exists(programname) :
        return ExternalProgram.get_path(programname) != None
    
    @staticmethod
    def get_path(programname) :
        for p in os.environ['PATH'].split(os.pathsep) :
            progpath = join(p, programname)
            if os.path.isfile(progpath) :
                # there may be another executable with the correct
                # permissions lower down in the path, but the shell
                # would not find it, so just return here...
                if os.access(progpath, os.X_OK) :
                    return progpath
                else :
                    return None
        
        return None
    
    def system(self, command) :
        self.log.debug(command)
        retcode = os.system(command) 
        if retcode != 0 :
            raise ExternalProgramError("'%s' return code %d" % (command, retcode))

class Sff2Fastq(ExternalProgram) :
    def __init__(self) :
        super(Sff2Fastq, self).__init__('sff2fastq')
        self.command = "sff2fastq -o %s %s 2> /dev/null"

    def run(self, sff, outdir) :
        if not isinstance(sff, SffFile) :
            raise ExternalProgramError("argument is not an SffFile")

        fastq_fname = outdir + os.sep + sff.get_basename() + ".fastq"

        self.log.info("running sff2fastq (%s, %s)" % (sff.get_filename(), fastq_fname))

        try :
            self.system(self.command % (fastq_fname, sff.get_filename()))
        
        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)
        
        return FastqFile(fastq_fname)

class GetMID(object) :
    def __init__(self, length) :
        self.length = length
        self.command = "grep -B1 \"^+\" %s | grep -v \"^[+-]\" | awk '{ print substr($0, 0, " + str(self.length) + ") }' | sort | uniq -c | sort -g | tail -1 | awk '{ print $2 }'"
        self.log = logging.getLogger('seance')

    def run(self, fastq_name) :
        status,output = commands.getstatusoutput(self.command % fastq_name)

        if status != 0 :
            raise ExternalProgramError("%s: %s", type(self).__name__, output)

        output = output.strip()

        # if the file is empty, there is nothing to read, so
        # returning an empty string will not matter...
        if output == '':
            self.log.info("mid empty")
            return output

        if re.match("[GATC]{%d}" % self.length, output) == None :
            self.log.error("%s does not look like a MID" % (output))
            sys.exit(1)
            #raise ExternalProgramError("%s: %s does not look like a MID" % (type(self).__name__, output))

        self.log.info("mid = %s" % output)

        return output

class GetMID2(object) :
    def __init__(self, length) :
        self.length = length
        self.log = logging.getLogger('seance')

    def run(self, fastq) :
        count = collections.Counter()
        fastq.open()
        
        for s in fastq :
            count[s[:self.length]] += s.duplicates

        fastq.close()

        if len(count) == 0 :
            return ""

        mid,midcount = count.most_common()[0]
        self.log.info("mid = %s" % mid)
        return mid

class Pagan(ExternalProgram) :
    def __init__(self) :
        super(Pagan, self).__init__('pagan')

    def get_alignment(self, fasta_fname) :
        command = "pagan --use-consensus --use-duplicate-weights --pileup-alignment --queryfile %s --outfile %s &> /dev/null"
        out_fname = fasta_fname + ".out"

        try :
            self.system(command % (fasta_fname, out_fname))

        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)

        # pagan always adds '.fas'
        return FastqFile(out_fname + ".fas")

    def get_454_alignment(self, fasta_fname) :
        command = "pagan --use-consensus --use-duplicate-weights --homopolymer --pileup-alignment --queryfile %s --outfile %s &> /dev/null"
        out_fname = fasta_fname + ".out"

        try :
            self.system(command % (fasta_fname, out_fname))

        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)
        
        # pagan always adds '.fas'
        return FastqFile(out_fname + ".fas")

    def phylogenetic_alignment(self, fasta_fname) :
        command = "pagan --seqfile %s --outfile %s --raxml-tree --xml" # &> /dev/null"
        out_fname = fasta_fname + ".out"

        try :
            self.system(command % (fasta_fname, out_fname))

        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)

        return out_fname + ".fas", out_fname + ".tre", out_fname + ".xml"

    def phylogenetic_placement(self, ref_alignment, ref_tree, queries) :
        command = "pagan --ref-seqfile %s --ref-treefile %s --queryfile %s --outfile %s --fast-placement \
                                --test-every-node --exhaustive-placement --one-placement-only --homopolymer --xml" #--output-nhx-tree"
        #command = "pagan --ref-seqfile %s --ref-treefile %s --queryfile %s --outfile %s \
        #                         --test-every-node --exhaustive-placement --one-placement-only --homopolymer --xml"
        out_fname = queries + ".placement"

        try :
            self.system(command % (ref_alignment, ref_tree, queries, out_fname))

        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)

        return out_fname + ".fas"

    def silva_phylogenetic_alignment(self, ref_alignment, ref_tree, queries) :
        #command = "pagan --ref-seqfile %s --ref-treefile %s --queryfile %s --outfile %s " + \
        #   "--fast-placement --use-anchors --prune-extended-alignment --test-every-terminal-node " + \
        #   "--xml --trim-extended-alignment --score-only-ungapped"
        #command = "pagan --ref-seqfile %s --ref-treefile %s --queryfile %s --outfile %s " + \
        #    "--use-anchors --use-exonerate-local --test-every-terminal-node --query-distance 0.01 " + \
        #    "--one-placement-only --output-nhx-tree --xml --prune-extended-alignment --trim-extended-alignment " + \
        #    "--prune-keep-number 0"
        command = "pagan --ref-seqfile %s \
                         --ref-treefile %s \
                         --queryfile %s \
                         --outfile %s \
                         --terminal-nodes \
                         --one-placement-only \
                         --output-nhx-tree \
                         --xml \
                         --trim-extended-alignment \
                         --prune-keep-number 0 \
                         --prune-extended-alignment"

        #out_fname = os.path.splitext(queries)[0] + ".silva"
        out_fname = queries + ".silva"

        try :
            self.system(command % (ref_alignment, ref_tree, queries, out_fname))

        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)

        for f in [ queries + '.silva.' + i for i in ['fas','nhx_tree','xml'] ] :
            os.remove(f)

        return out_fname + ".pruned.fas", out_fname + ".pruned.tre", out_fname + ".pruned.xml"

class Aligner1D(object) :
    def __init__(self) :
        self.__reset()

    def __reset(self) :
        self.sequences = []
        self.charfreqs = []
        self.maxlens = {}

    def __test_max(self, pos, freq) :
        try :
            if self.maxlens[pos] < freq :
                self.maxlens[pos] = freq

        except KeyError, ke :
            self.maxlens[pos] = freq

    def __convert(self, seq) :
        tmp = []

        char = seq[0]
        freq = 1

        for i in seq[1:] :
            if i == char :
                freq += 1
            else :
                self.__test_max(len(tmp), freq)
                tmp.append((char, freq))

                char = i
                freq = 1

        self.__test_max(len(tmp), freq)
        tmp.append((char, freq))

        return tmp

    def get_alignment(self, fname) :
        self.__reset()

        fq = FastqFile(fname)

        fq.open()

        for seq in fq :
            self.sequences.append(seq)
            self.charfreqs.append(self.__convert(seq))

        fq.close()

        s = ""
        for i in range(len(self.charfreqs)) :
            s += (">seq%s\n" % self.sequences[i].id)
            tmp = self.charfreqs[i]
            for pos in range(len(self.maxlens)) :
                try :
                    char,freq = tmp[pos]
                except IndexError, ie :
                    char,freq = 'X', 0    

                s += ((char * freq) + ('-' * (self.maxlens[pos] - freq)))
            s += "\n"

        # write out again to file to keep same API
        # TODO: streamline, this is mostly for testing so I don't need to
        # change anything else
        outf = open(fname + ".out", 'w')
        print >> outf, s
        outf.close()

        return FastqFile(outf.name)

class Uchime(ExternalProgram) :
    def __init__(self) :
        super(Uchime, self).__init__('uchime')
        self.command = "uchime --input %s --uchimeout %s &> /dev/null"

    def __parse(self, fname) :
        tmp = []
        f = open(fname)

        for line in f :
            line = line.strip()

            if line == "" :
                continue

            data = line.split()
            
            if data[16] == 'Y' :
                # data[1] = "seqXXX/ab=YYY"
                #seqid = int(data[1].split('/')[0][3:])
                seqid = int(data[1].split('/')[0])
                tmp.append(seqid)

        f.close()

        return tmp

    def run(self, fname) :
        out_fname = fname + ".uchime"

        try :
            self.system(self.command % (fname, out_fname))

        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)

        return self.__parse(out_fname)

class BlastN(ExternalProgram) :
    def __init__(self) :
        super(BlastN, self).__init__('blastn')
        self.command = "blastn -query %s -db nr -remote -num_alignments 10 -outfmt 10"
        self.url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta"

    def __get_complete_desc(self, name) :
        try :
            f = urllib2.urlopen(self.url % name, None, 1)
            #return '_'.join(f.readline().split('|')[-1].strip().split()[:2]).replace('.', '') # ;-P
            line = f.readline()
            tmp = '_'.join(line.split('|')[-1].strip().split()[:2])
            for badchar in " \n\t:,)(;][" : # these cause problems for RaxML
                tmp = tmp.replace(badchar, '')

            if tmp == '' :
                self.log.warn("querying ncbi eutils for %s failed (empty file), retrying..." % (name))
                return self.__get_complete_desc(name)

            return tmp
        
        except urllib2.HTTPError, he :
            self.log.warn("querying ncbi eutils for %s failed (%s), retrying..." % (name, str(he)))
            return self.__get_complete_desc(name)

        except urllib2.URLError, ue :
            self.log.warn("querying ncbi eutils for %s failed (%s), retrying..." % (name, str(ue)))
            return self.__get_complete_desc(name)

    def get_names(self, fasta_fname) :
        s,o = commands.getstatusoutput(self.command % fasta_fname)
        
        if s != 0 :
            self.log.error("blastn returned %d" % s)
            sys.exit(1)

        names = {}

        for line in o.split('\n') :
            fields = line.split(',')
            try :
                name = int(fields[0])

            except ValueError, ve :
                self.log.warn("blast returned an id I do not understand (%s), skipping..." % (fields[0]))
                continue

            # only accept the first one
            if name in names :
                continue
            
            try :
                desc = fields[1].split('|')

            except IndexError :
                self.log.warn("could not split line from blastn: %s" % str(fields))
                continue

            if re.match(".+\.\d+", desc[3]) :
                names[name] = "%s_%s_%s" % (fields[0], self.__get_complete_desc(desc[3]), fields[2])

        return names

class PyroDist(ExternalProgram) :
    def __init__(self) :
        super(PyroDist, self).__init__('PyroDist')
        self.command = "PyroDist -in %s -out %s -rin %s &> /dev/null"

    def __lookup_file(self) :
        f = join(dirname(ExternalProgram.get_path(self.programname)), "LookUp.dat")

        if not os.path.exists(f) :
            raise ExternalProgramError("%s - cannot find LookUp.dat (it needs to be in the same dir as %s)" % (self.programname, self.programname))

        return f

    def run(self, datfile, outfile) :
        try :
            self.system(self.command % (datfile, outfile, self.__lookup_file()))

        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)

        return "%s.fdist" % outfile
    
class FCluster(ExternalProgram) :
    def __init__(self) :
        super(FCluster, self).__init__('FCluster')
        self.command = "FCluster -in %s -out %s &> /dev/null"

    def run(self, fdistfile, outfile) :
        try :
            self.system(self.command % (fdistfile, outfile))

        except ExternalProgramError, epe :
            self.log.error(str(epe))            
            sys.exit(1)

        return "%s.list" % outfile

class PyroNoise(ExternalProgram) :
    def __init__(self) :
        super(PyroNoise, self).__init__('PyroNoise')
        self.command = "PyroNoise -din %s -out %s -lin %s -rin %s -s 60.0 -c 0.01 &> /dev/null"

    def __lookup_file(self) :
        f = join(dirname(ExternalProgram.get_path(self.programname)), "LookUp.dat")

        if not os.path.exists(f) :
            raise ExternalProgramError("%s - cannot find LookUp.dat (it needs to be in the same dir as %s)" % (self.programname, self.programname))

        return f

    def run(self, datfile, listfile, outfile) :
        try :
            self.system(self.command % (datfile, outfile, listfile, self.__lookup_file()))

        except ExternalProgramError, epe :
            self.log.error(str(epe))            
            sys.exit(1)

        return "%s_cd.fa" % outfile

class AmpliconNoise(ExternalProgram) :
    def __init__(self) :
        super(AmpliconNoise, self).__init__('AmpliconNoise')

    def close_enough_old(self, primer, sequence, errors) :
        err_count = 0

        for i,j in zip(sequence, primer) :
            if not IUPAC.equal(i, j) :
                err_count += 1

        #if len(a) > 10 :
        #    print a, b, err_count

        return err_count <= errors

    def close_enough(self, primer, sequence, diff) :
        if diff < 0 :
            return False

        if (len(primer) == 0) or (len(sequence) == 0) :
            return True

        m = IUPAC.equal(sequence[0], primer[0])

        return self.close_enough(primer[1:], sequence[1:], diff if m else diff-1) or \
               self.close_enough(primer[1:], sequence, diff-1) or \
               self.close_enough(primer, sequence[1:], diff-1)

    def extract(self, sff, outdir, primer, barcode, barcode_errors, max_homopolymer) :
        try :
            from Bio import SeqIO
        except ImportError :
            print >> sys.stderr, "BioPython not installed (only required for working with SFF files)"
            sys.exit(1)

        barcode_len = len(barcode)
        primer_len = len(primer)

        raw_seq_total = 0

        names = []
        flows = []
        flowlens = []

        for record in SeqIO.parse(sff.get_filename(), "sff") :
            raw_seq_total += 1
            good_bases = record.seq[record.annotations["clip_qual_left"] : record.annotations["clip_qual_right"]]
            barcode_seq = good_bases[:barcode_len]
            primer_seq = good_bases[barcode_len : barcode_len + primer_len]

            new_length = 0

            for i in range(0, len(record.annotations["flow_values"]), 4) : 
                signal = 0
                noise = 0

                for j in range(4) :
                    f = float(record.annotations["flow_values"][i + j]) / 100.0

                    if int(f + 0.5) > max_homopolymer :
                        break

                    if f > 0.5 :
                        signal += 1
                        if f < 0.7 :
                            noise += 1

                if noise > 0 or signal == 0 :
                    break

                new_length += 1

            new_length *= 4

            if new_length > 450 :
                new_length = 450

            if new_length >= 360 and \
                    IUPAC.close_enough(barcode, barcode_seq, barcode_errors) and \
                    IUPAC.close_enough(primer, primer_seq, 2) :
                flows.append(record.annotations["flow_values"])
                flowlens.append(new_length)
                names.append(record.id)



        if len(flows) == 0 :
            self.log.info("kept 0/%d sequences" % raw_seq_total)
            return 0, None

        # output pyronoise input file
        # see http://userweb.eng.gla.ac.uk/christopher.quince/Software/PyroNoise.html
        f = open(join(outdir, "flows.dat"), 'w')

        print >> f, "%d %d" % (len(flows), max([ len(i) for i in flows ]))
        for i in range(len(flows)) :
            print >> f, " ".join([ names[i], str(flowlens[i]) ] + [ "%.2f" % (float(i) / 100.0) for i in flows[i] ])

        f.close()

        self.log.info("kept %d/%d sequences" % (len(flows), raw_seq_total))
        return len(flows), f.name

    def run(self, sff, outdir, forward_primer, barcode, barcode_errors, max_homopolymer) :
        if not isinstance(sff, SffFile) :
            raise ExternalProgramError("argument is not an SffFile")

        output_name = abspath(join(outdir, sff.get_basename() + '.fasta'))

        numseq,fname = self.extract(sff, outdir, forward_primer, barcode, barcode_errors, max_homopolymer)

        # just so the rest of the pipeline can be run and there be a record
        # of the sample containing zero sequences
        if numseq == 0 :
            open(output_name, 'w').close()
            return FastqFile(output_name)

        # well... this is a mess
        # PyroDist does not like being given 1 sequence
        if numseq == 1 :
            try :
                from Bio import SeqIO
            except ImportError: 
                print >> sys.stderr, "BioPython not installed (only required for working with SFF files)"
                sys.exit(1)

            fout = open(output_name, 'w')

            for r in SeqIO.parse(sff.get_filename(), 'sff-trim') :
                print >> fout, ">seq0 NumDuplicates=1\n%s" % r.seq
            
            fout.close()
            
            return FastqFile(output_name)

        outfile = join(outdir, "flows")
        
        distfile = PyroDist().run(fname, outfile)
        listfile = FCluster().run(distfile, outfile)
        fafile   = PyroNoise().run(fname, listfile, outfile)

        # read fa file
        # add NumDulicates fields
        # output with correct file name
        fout = open(output_name, 'w')
        f = FastqFile(outfile + "_cd.fa")

#        seq2qual = {}
#        q = open(outfile + "_cd.qual")
#        seqname = None
#        for line in q :
#            if line.startswith('>') :
#                seqname = line.rstrip()[1:]
#            else :
#                seq2qual[seqname] = ''.join([ Sequence.int_to_quality(int(i)) for i in line.split() ])
#        q.close()

        f.open()

        count = 0
        for seq in f :
            dups = int(seq.id.split('_')[-1])
            print >> fout, ">seq%d NumDuplicates=%d\n%s" % (count, dups, seq.sequence)
            count += 1

        f.close()
        fout.close()

        # delete intermediate files
        shutil.rmtree(outfile)
        for fname in glob.glob(outfile + '*') :
            os.remove(fname)

        return FastqFile(output_name)

class PyroNoiseMothur(ExternalProgram) :
    def __init__(self) :
        super(PyroNoiseMothur, self).__init__('PyroNoise')
        self.command = "mothur \"#sff.multiple(file=tmp.txt, maxhomop=%d, pdiffs=0, bdiffs=%d, minflows=360)\" &> /dev/null"

    def delete_intermediates(self, prefix) :
        suffix = ['.fasta', '.flow', '.flow.files', '.qual', '.scrap.flow', '.summary', '.trim.flow']
        files = ['tmp.oligos', 'tmp.txt', 'not found.shhh.fasta', 'not found.shhh.names'] + glob.glob("mothur.*.logfile")
        for f in [ prefix + i for i in suffix ] + files :
            try :
                os.remove(f)
            except OSError :
                pass

    def run(self, sff, outdir, forward_primer, barcode, barcode_errors, max_homopolymers) :
        # 1. ensure SFF file name does not contain hyphens
        # 2. write singular.txt
        #   eg: Tg_2_25062012_1.sff singular.oligos
        # 3. write singular.oligos
        #   eg: forward AGRGGTGAAATYCGTGGAC
        #       barcode TACAG Tg_2_25062012_1
        # 4. run mothur "#sff.multiple(file=singular.txt, maxhomop=8, pdiffs=2, bdiffs=1)"
        # 5. output : singular.singular.fasta   - rename
        #             singular.singular.groups  - kill
        #             singular.singular.names   - kill
        
        if not isinstance(sff, SffFile) :
            raise ExternalProgramError("argument is not an SffFile")

        cwd = os.getcwd()
        
        old_sff_name = abspath(sff.get_filename())
        new_sff_name = sff.get_basename().replace('-', '_')
        output_name = abspath(join(outdir, sff.get_basename() + '.fasta'))
        sff_prefix,sff_ext = os.path.splitext(new_sff_name)

        os.chdir(outdir)

        #print os.getcwd()
        #print old_sff_name
        #print new_sff_name        

        try :
            os.symlink(old_sff_name, new_sff_name)
        except OSError :
            pass
            #print >> sys.stderr, "could not create symlink"

        # create batch file
        f = open('tmp.txt', 'w')
        print >> f, "%s tmp.oligos" % new_sff_name
        f.close()

        # create oligos file
        f = open('tmp.oligos', 'w')
        print >> f, "forward %s" % forward_primer
        print >> f, "barcode %s %s" % (barcode, new_sff_name)
        f.close()

        # run mothur
        if os.system(self.command % (max_homopolymers, barcode_errors)) != 0 :
            # mothur segfaults on empty files (or is it files where nothing
            # passes qc?), just create an empty file and return
            open(output_name, 'w').close()
            os.remove(new_sff_name)
            self.delete_intermediates(sff_prefix)
            os.chdir(cwd)
            return FastqFile(output_name)

        # rename to what we want to call it
        os.rename('tmp.tmp.fasta', output_name)

        # 6.
        # merge output files to get number of duplicates 
        # in fasta file
        counts = {}
        f = open(new_sff_name[:-3] + 'shhh.trim.summary')
        f.readline()
        for line in f :
            seqname,start,end,nbases,ambigs,polymer,numseqs = line.strip().split()
            counts[seqname] = int(numseqs)
        f.close()

        f = FastqFile(new_sff_name[:-3] + 'shhh.trim.fasta')
        out_fasta = open(output_name, 'w')
        f.open()
        for seq in f :
            print >> out_fasta, "%s NumDuplicates=%d\n%s\n" % (seq.id, counts[seq.id[1:]], seq.sequence)
        f.close()
        out_fasta.close()

        self.delete_intermediates(sff_prefix)
        os.chdir(cwd)

        return FastqFile(output_name)

class Cutadapt(ExternalProgram) :
    def __init__(self) :
        super(Cutadapt, self).__init__('cutadapt')

    def run(self, fastq, outdir, forwardprimer, reverseprimer=None) :
        command = "cutadapt"

        if (forwardprimer == None) and (reverseprimer == None) :
            return fastq

        if forwardprimer != None :
            command += (" -b %s" % forwardprimer)

        if reverseprimer != None :
            command += (" -a %s" % reverseprimer)

        #command += " --trimmed-only -e 0.0"

        outfile = "%s.adaptortrim" % join(outdir, fastq.get_basename())
        command += (" %s > %s 2> /dev/null" % (fastq.get_filename(), outfile))

        try :
            self.system(command)

        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)

        return FastqFile(outfile)

if __name__ == '__main__' :
    log = logging.getLogger('seance')
    log.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter('%(levelname)s %(message)s'))

    log.addHandler(ch)

    # ---

    #p = PyroNoiseMothur()
    p = AmpliconNoise() 
   
    sff_filename = "/Users/ajm/research/coding/mothra/test_files/Tg_12-25062012-1.sff"
    outdir = "."
    primer = "AGRGGTGAAATYCGTGGAC"
    barcode = "CATGC"
    barcode_errors = 0
    max_homopolymers = 8

    p.run(SffFile(sff_filename), outdir, primer, barcode, barcode_errors, max_homopolymers)

