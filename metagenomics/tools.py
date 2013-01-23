import sys
import abc
import os
import commands
import re

from metagenomics.filetypes import SffFile, FastqFile

class ExternalProgramError(Exception) :
    pass

class ExternalProgramNotInstalledError(Exception) :
    pass

class ExternalProgram(object) :
    __metaclass__ = abc.ABCMeta

    def __init__(self, pname) :
        self.programname = pname
    
    @staticmethod
    def exists(programname) :
        for p in os.environ['PATH'].split(os.pathsep) :
            progpath = os.path.join(p, programname)
            if os.path.isfile(progpath) :
                # there may be another executable with the correct
                # permissions lower down in the path, but the shell
                # would not find it, so just return here...
                return os.access(progpath, os.X_OK) 
        
        return False
    
    def system(self, command) :
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

        try :
            self.system(self.command % (fastq_fname, sff.get_filename()))
        
        except ExternalProgramError, epe :
            print >> sys.stderr, "Error: " + str(epe)
            sys.exit(-1)
        
        return FastqFile(fastq_fname)

class GetMID(object) :
    def __init__(self) :
        self.command = "grep -A1 \"^@\" %s | grep -v \"^[@-]\" | awk '{ print substr($0, 0, 10) }' | sort | uniq -c | sort -g | tail -1 | awk '{ print $2 }'"

    def run(self, fastq_name) :
        status,output = commands.getstatusoutput(self.command % fastq_name)

        if status != 0 :
            raise ExternalProgramError("%s: %s", type(self).__name__, output)

        output = output.strip()
        if re.match("[GATC]{10}", output) == None :
            raise ExternalProgramError("%s: %s does not look like a MID", type(self).__name__, output)

        return output

class Pagan(ExternalProgram) :
    def __init__(self) :
        super(Pagan, self).__init__('pagan')
        #self.command = "pagan --use-consensus --consensus-minimum=3 --use-duplicate-weigths --454 --queryfile %s --outfile %s &> /dev/null"
        #self.command = "pagan --use-consensus --use-duplicate-weigths --homopolymer --pileup-alignment --use-prefix-anchors --no-terminal-edges --queryfile %s --outfile %s &> /dev/null"
        self.command = "pagan --use-consensus --use-duplicate-weigths --homopolymer --pileup-alignment --queryfile %s --outfile %s &> /dev/null"

    def get_454_alignment(self, fasta_fname) :
        out_fname = fasta_fname + ".out"

        try :
            self.system(self.command % (fasta_fname, out_fname))

        except ExternalProgramError, epe :
            print >> sys.stderr, "Error: " + str(epe)
            sys.exit(-1)
        
        # pagan always adds '.fas'
        return FastqFile(out_fname + ".fas")

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
            s += ("%s\n" % self.sequences[i].id)
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

