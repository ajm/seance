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
    
    @abc.abstractmethod
    def run(self, *args) :
        pass

    def system(self, command) :
        retcode = os.system(command) 
        if retcode != 0 :
            raise ExternalProgramError("'%s' return code %d" % (command, retcode))

class Sff2Fastq(ExternalProgram) :
    def __init__(self) :
        ExternalProgram.__init__(self, 'sff2fastq')
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

class TagCleaner(ExternalProgram) :
    def __init__(self) :
        ExternalProgram.__init__(self, 'tagCleaner')
        self.command = ""

    def run(self, *args) :
        raise NotImplemented


class Pagan(ExternalProgram) :
    def __init__(self) :
        super(Pagan, self).__init__('pagan')

    def run(self, *args) :
        raise NotImplemented

