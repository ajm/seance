import sys
import abc
import os

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

class Pagan(ExternalProgram) :
    def __init__(self) :
        super(Pagan, self).__init__('pagan')
        self.command = "pagan --use-consensus --consensus-minimum=3 --use-duplicate-weigths --454 --queryfile %s --outfile %s &> /dev/null"

    def get_454_alignment(self, fasta_fname) :
        out_fname = fasta_fname + ".out"

        try :
            self.system(self.command % (fasta_fname, out_fname))

        except ExternalProgramError, epe :
            print >> sys.stderr, "Error: " + str(epe)
            sys.exit(-1)
        
        # pagan always adds '.fas'
        return FastqFile(out_fname + ".fas")
