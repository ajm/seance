import sys
import os
import tempfile
import multiprocessing

from metagenomics.tools import ExternalProgram

class System(object) :
    __tmpdir__ = None

    def __init__(self) :
        pass

    @staticmethod
    def tempdir(tmpdir=None) :
        if tmpdir != None :
            System.__tmpdir__ = tmpdir
            return

        if System.__tmpdir__ != None :
            return System.__tmpdir__
        elif os.environ.has_key('TMPDIR') :
            return os.environ['TMPDIR']
        else :
            return '/tmp'

    @staticmethod
    def tempfilename(ext="") :
        return tempfile.mktemp(dir=System.tempdir()) + ext

    def check_local_installation(self, required_programs) :
        bad = False

        # ensure required programs are installed locally
        for prog in required_programs :
            if not ExternalProgram.exists(prog) :
                print >> sys.stderr, "Error: '%s' is not installed." % prog
                bad = True

        if bad :
            print >> sys.stderr, "Error: Exiting due to required external programs not being installed..."
            sys.exit(-1)

    def check_file(self, filename) :
        if not os.path.exists(filename) :
            print >> sys.stderr, "Error: '%s' does not exist." % filename
            return False

        if not os.path.isfile(filename) :
            print >> sys.stderr, "Error: '%s' exists, but is not a file!" % filename
            return False

        return True

    def check_directory(self, dirname, create=False) :
        if not os.path.exists(dirname) :
            if create :
                print >> sys.stderr, "Info: '%s' does not exist, creating..." % dirname
                try :
                    os.makedirs(dirname)
                except OSError, ose :
                    print >> sys.stderr, str(ose)
                    return False
            else :
                print >> sys.stderr, "Error: '%s' does not exist." % dirname
                return False

        elif not os.path.isdir(dirname) :
            print >> sys.stderr, "Error: '%s' exists, but is not a directory!" % dirname
            return False

        else :
            return True

    def check_directories(self, dirs) :
        return False not in map(lambda x : self.check_directory(x[0], x[1]), dirs)

    def run_commands_multithreaded(self, commands) :
        p = multiprocessing.Pool(multiprocessing.cpu_count())

        retcodes = p.map(os.system, commands)

        bad = False
        for i in range(len(retcodes)) :
            if retcodes[i] != 0 :
                print >> sys,stderr, "Error: '%s' failed" % commands[i]
                bad = True

        if bad :
            sys.exit(-1)


