import sys
import os
import tempfile
import multiprocessing
import logging

from seance.tools import ExternalProgram


class System(object) :
    __tmpdir__ = None

    def __init__(self) :
        self.log = logging.getLogger('seance')

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
        fd,fpath = tempfile.mkstemp(dir=System.tempdir(), suffix=ext)
        os.close(fd)
        return fpath

    @staticmethod
    def is_installed(p) :
        return ExternalProgram.get_path(p) != None

    def check_local_installation(self, required_programs) :
        bad = False

        # ensure required programs are installed locally
        for prog in required_programs :
            path = ExternalProgram.get_path(prog)

            if path != None :
                self.log.info("found %s @ %s" % (prog, path))
            else :
                self.log.error("%s is not installed" % prog)
                bad = True

        if bad :
            self.log.error("Exiting due to required external programs not being installed...")
            sys.exit(1)

    def check_file(self, filename) :
        if not os.path.exists(filename) :
            self.log.error("%s does not exist" % filename)
            return False

        if not os.path.isfile(filename) :
            self.log.error("%s exists, but is not a file" % filename)
            return False

        return True

    def check_directory(self, dirname, create=False) :
        if not os.path.exists(dirname) :
            if create :
                self.log.info("%s does not exist, creating..." % dirname)
                try :
                    os.makedirs(dirname)
                    return True

                except OSError, ose :
                    self.log.error(str(ose))
                    return False
            else :
                self.log.error("%s does not exist" % dirname)
                return False

        elif not os.path.isdir(dirname) :
            self.log.error("%s exists, but is not a directory" % dirname)
            return False

        else :
            return True

    def check_files(self, filenames) :
        return False not in map(lambda x : self.check_file(x), filenames)

    def check_directories(self, dirs) :
        return False not in map(lambda x : self.check_directory(x[0], x[1]), dirs)

    def run_commands_multithreaded(self, commands) :
        p = multiprocessing.Pool(multiprocessing.cpu_count())

        retcodes = p.map(os.system, commands)

        bad = False
        for i in range(len(retcodes)) :
            if retcodes[i] != 0 :
                self.log.error("%s failed" % commands[i])
                bad = True

        if bad :
            sys.exit(1)

