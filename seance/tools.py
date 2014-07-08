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
import socket

from os.path import abspath, join, dirname
from seance.filetypes import SffFile, FastqFile
from seance.datatypes import IUPAC
from seance.progress import Progress


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
    
    def system(self, command, debug=True) :
        if debug :
            self.log.debug(command)
        retcode = os.system(command) 
        if retcode != 0 :
            raise ExternalProgramError("'%s' return code %d" % (' '.join(command.split()), retcode))

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
        self.log.debug("mid = %s (%s)" % (mid, fastq.get_filename()))

        if re.match("[GATC]{%d}" % self.length, mid) == None :
            self.log.error("%s does not look like a MID" % (mid))
            sys.exit(1)

        return mid

class Pagan(ExternalProgram) :
    def __init__(self) :
        super(Pagan, self).__init__('pagan')

    def get_alignment(self, fasta_fname) :
        command = "pagan --use-consensus --use-duplicate-weights --pileup-alignment --queryfile %s --outfile %s &> /dev/null"
        out_fname = fasta_fname + ".out"

        try :
            self.system(command % (fasta_fname, out_fname), debug=False)

        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)

        # pagan always adds '.fas'
        return FastqFile(out_fname + ".fas")

    def get_454_alignment(self, fasta_fname) :
        command = "pagan --use-consensus --use-duplicate-weights --homopolymer --pileup-alignment --queryfile %s --outfile %s &> /dev/null"
        out_fname = fasta_fname + ".out"

        try :
            self.system(command % (fasta_fname, out_fname), debug=False)

        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)
        
        # pagan always adds '.fas'
        return FastqFile(out_fname + ".fas")

    def phylogenetic_alignment(self, fasta_fname, tree_fname=None) :
        out_fname = fasta_fname + ".out"

        if tree_fname :
            command = "pagan --seqfile %s --treefile %s --outfile %s --xml" % (fasta_fname, tree_fname, out_fname)
        else :
            command = "pagan --seqfile %s --raxml-tree --outfile %s --xml" % (fasta_fname, out_fname)

        try :
            self.system(command)

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
                         --prune-extended-alignment \
                         --prune-keep-closest"

        #out_fname = os.path.splitext(queries)[0] + ".silva"
        out_fname = queries + ".silva"

        try :
            self.system(command % (ref_alignment, ref_tree, queries, out_fname))

        except ExternalProgramError, epe :
            self.log.error(str(epe))
            sys.exit(1)

        for f in [ queries + '.silva.' + i for i in ['fas','nhx_tree','xml'] ] :
            os.remove(f)

        return out_fname + ".pruned_closest.fas", out_fname + ".pruned_closest.tre", out_fname + ".pruned_closest.xml"

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
    def __init__(self, verbose) :
        super(BlastN, self).__init__('blastn')
        self.command = "blastn -query %s -db nr -remote -task megablast -outfmt 10 -perc_identity 90 -entrez_query 'all[filter] NOT (environmental samples[organism] OR metagenomes[orgn])'"
        self.url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=xml"
        self.regex = {  "blast"     : "Org\-ref_taxname",
                        "taxonomy"  : "OrgName_lineage" }
        self.verbose = verbose

    def __merge_taxonomy(self, names) :
        tmp = collections.defaultdict(list)
        for name in names :
            for i,v in enumerate(name.split(';')) :
                if v == 'error' :
                    return v

                tmp[i].append(v)
        
        s = ""
        for i in sorted(tmp.keys()) :
            if len(set(tmp[i])) == 1 :
                s += (";%s" % tmp[i][0])
            else :
                break
                #s += (";{%s}" % ','.join(sorted(tmp[i])))
        
        if len(s) == 0 :
            return "cannot label (matches multiple domains!)"

        return s[1:]

    def __get_desc_common(self, names, method, progress) :
        #return self.__merge_taxonomy([ self.__get_desc(i[0], method) for i in names ])

        tmp = []

        for i in names :
            tmp.append( self.__get_desc_wrapper(i[0], method, 10) )
            progress.increment()

        return self.__merge_taxonomy(tmp)

    def __get_desc_concat(self, names, method) :
        return ','.join([ self.__get_desc(i[0], method) for i in names ])

    def __get_desc_wrapper(self, name, method, max_attempts, return_on_failure='error') :
        attempts = 0

        while True :
            try :
                return self.__get_desc(name, method)
            except :
                attempts += 1
                
            if attempts == max_attempts :
                print "" # there is a progress bar
                self.log.error("tried to query '%s' using NCBI eutils (exceeded %d attempts)" % (name, max_attempts))
                return return_on_failure

    def __get_desc(self, name, method) :
        try :
            f = urllib2.urlopen(self.url % name, None, 5)
            tmp = None

            for line in f :
                m = re.match("\W*<%s>(.*)</%s>" % (self.regex[method], self.regex[method]), line) 
                if m :
                    tmp = m.group(1)
                    break

            f.close()

            if tmp is None :
                self.log.error("Error calling NCBI eutils with '%s'" % name)
                tmp = "none"

            if method == 'taxonomy' :
                tmp = tmp.replace(" ", "")
            
            #self.log.info(tmp)
            return tmp

        except urllib2.HTTPError, he :
            self.log.debug("querying ncbi eutils for %s failed (%s), retrying..." % (name, str(he)))
            raise he

        except urllib2.URLError, ue :
            self.log.debug("querying ncbi eutils for %s failed (%s), retrying..." % (name, str(ue)))
            raise ue

        except socket.timeout, to :
            self.log.debug("querying ncbi eutils for %s failed (%s), retrying..." % (name, str(to)))
            raise to

    def get_names(self, fasta_fname, method) :
        #self.log.info("getting OTU names (this may take a while)...")

        if method not in self.regex.keys() :
            self.log.error("'%s' is not a valid labelling method" % method)
            sys.exit(1)

        s,o = commands.getstatusoutput(self.command % fasta_fname)
        
        if s != 0 :
            self.log.error("blastn returned %d" % s)
            sys.exit(1)

        names = collections.defaultdict(list)
        scores = {}

        for line in o.split('\n') :
            fields = line.split(',')
            try :
                name = int(fields[0])

            except ValueError :
                name = fields[0]

            try :
                score = float(fields[-1])

            except ValueError, ve :
                self.log.warn("problem with blast result (%s), skipping..." % (str(ve)))
                self.log.debug(line)
                continue

            if name not in scores :
                scores[name] = score

            # only keep the highest scoring hits
            if scores[name] != score :
                continue

            try :
                desc = fields[1].split('|')

                if re.match(".+\.\d+", desc[3]) :
                    names[name].append((desc[3], fields[2]))
                    #names[name] = "%s_%s_%s" % (fields[0], self.__get_desc(desc[3], method), fields[2])

            except IndexError :
                self.log.warn("could not split line from blastn: %s" % str(fields))
                continue

        p = Progress("OTU naming (%s)" % method, sum([len(i) for i in names.values()]) if method == 'taxonomy' else len(names))
        p.start()

        # now generate labels
        if method == 'blast' :
            for name in names :
                tmp = names[name][0]
                names[name] = "%s_%s" % (self.__get_desc_wrapper(tmp[0], method, 10), tmp[1])
                p.increment()
                #names[name] = "%s_%s" % (self.__get_desc_concat(names[name], method), tmp[1])
        elif method == 'taxonomy' :
            for name in names :
                names[name] = "%s" % (self.__get_desc_common(names[name], method, p))

        p.end()

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

    def extract(self, sff, outdir, primer, primer_errors, barcode, barcode_errors, max_homopolymer) :
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
                    IUPAC.close_enough(primer, primer_seq, primer_errors) :
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

