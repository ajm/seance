import sys
import os
import getopt
import multiprocessing
import signal

from metagenomics.workflow import WorkFlow
from metagenomics.tools import ExternalProgram

def err(s) :
    print >> sys.stderr, "Error: " + s

def out(s) :
    print "Info:  " + s

def get_default_options() :
    return {
            "datadir" : None,
            "tempdir" : None,
            "metadata" : None,
            "verbose" : True
           }

def get_required_programs() :
    return [
            "sff2fastq"
            ]

def clean_up() :
    pass

def handler_sigterm(signal, frame) :
    clean_up()
    sys.exit(0)

def check_local_installation() :
    bad = False

    # ensure required programs are installed locally
    for prog in get_required_programs() :
        if not ExternalProgram.exists(prog) :
            err("'%s' is not installed." % prog)
            bad = True

    if bad :
        err("Exiting due to required external programs not being installed...")
        sys.exit(-1)

def init() :
    # install a signal handler to die cleanly
    signal.signal(signal.SIGINT, handler_sigterm)

    # check that necessary software is installed
    check_local_installation()

def usage(progname) :
    print >> sys.stderr, "Usage: %s [OPTIONS]" % progname
    
    options = get_default_options()

    for k in options :
        print >> sys.stderr, "\t-%c\t--%s\t(default = %s)" % (k[0], k, options[k])

    print >> sys.stderr, ""

def parse_args() :
    options = get_default_options()

    try :
        opts,args = getopt.getopt(
                        sys.argv[1:],
                        "d:t:m:h",
                        ["datadir=", "tempdir=", "metadata=", "help"]
                    )

    except getopt.GetoptError, err :
        print >> sys.stderr, str(err)
        usage(sys.argv[0])
        sys.exit(-1)

    for o,a in opts :
        if o in ('-h', '--help') :
            usage(sys.argv[0])
            sys.exit(0)

        elif o in ('-d', '--datadir') :
            options['datadir'] = a

        elif o in ('-t', '--tempdir') :
            options['tempdir'] = a

        elif o in ('-m', '--metadata') :
            options['metadata'] = a

        else :
            assert False, "unhandled option %s" % o

    return options

def check_dirs(options) :
    # check data dir
    tmp = options['datadir']
    if not os.path.exists(tmp) :
        err("specified data directory '%s' does not exist." % tmp)
        return False
    else :
        if not os.path.isdir(tmp) :
            err("specified data directory '%s' is not a directory." % tmp)
            return False

    # check tmp data
    tmp = options['tempdir']
    if not os.path.exists(tmp) :
        out("temp directory '%s' does not exist, creating..." % tmp)
        try :
            os.makedirs(tmp)

        except OSError, ose :
            err(str(ose))
            return False
    else :
        if not os.path.isdir(tmp) :
            err("specified temp directory '%s' exists, but it not a directory!" % tmp)
            return False

    # check metadata file
    tmp = options['metadata']
    if not os.path.isfile(tmp) :
        err("metadata file does not exist.")
        return False

    return True

def mandatory_options_set(options) :
    mandatory = ["datadir", "tempdir", "metadata"]
    ret = True
    
    for m in mandatory :
        if options[m] == None :
            err("%s must be set." % m)
            ret = False
    
    return ret

def run_commands_multithreaded(commands) :
    p = multiprocessing.Pool(multiprocessing.cpu_count())
    
    retcodes = p.map(os.system, commands)

    bad = False
    for i in range(len(retcodes)) :
        if retcodes[i] != 0 :
            err("'%s' failed" % commands[i])
            bad = True

    if bad :
        sys.exit(-1)

def main() :
    init()

    options = parse_args()
    
    if not mandatory_options_set(options) or not check_dirs(options) :
        sys.exit(-1)

    wf = WorkFlow(options['datadir'], options['tempdir'], options['metadata'])
    wf.run()

if __name__ == '__main__' :
    main()

