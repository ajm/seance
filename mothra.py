import sys
import os
import getopt
import multiprocessing
import signal

from metagenomics.workflow import WorkFlow
from metagenomics.system import System


def get_default_options() :
    return {
            "datadir"       : None,
            "tempdir"       : None,
            "metadata"      : None,
            "verbose"       : True,
            "minlength"     : 100,
            "minquality"    : 15,
            "winquality"    : 10
           }

def get_mandatory_options() :
    return ["datadir", "tempdir", "metadata"]

def get_required_programs() :
    return ["sff2fastq"]

def clean_up() :
    pass

def handler_sigterm(signal, frame) :
    clean_up()
    sys.exit(0)

def usage(progname) :
    print >> sys.stderr, "Usage: %s [OPTIONS]" % progname
    
    options = get_default_options()

    for k in sorted(options) :
        print >> sys.stderr, "\t-%c\t--%s\t(default = %s)" % (k[0], k, options[k])

    print >> sys.stderr, ""

def expect_int(parameter, argument) :
    try :
        return int(argument)

    except ValueError, ve :
        print >> sys.stderr, "Problem parsing argument for %s: %s\n" % (parameter, str(ve))
        usage(sys.argv[0])
        sys.exit(-1)

def parse_args() :
    options = get_default_options()

    try :
        opts,args = getopt.getopt(
                        sys.argv[1:],
                        "d:t:m:hl:q:w:",
                        ["datadir=", "tempdir=", "metadata=", "help", "minlength=", "minquality=", "winquality="]
                    )

    except getopt.GetoptError, err :
        print >> sys.stderr, str(err) + "\n"
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

        elif o in ('-l', '--minlength') :
            options['minlength'] = expect_int("minlength", a)

        elif o in ('-q', '--minquality') :
            options['minquality'] = expect_int("minquality", a)

        elif o in ('-w', '--winquality') :
            options['winquality'] = expect_int("winquality", a)

        else :
            assert False, "unhandled option %s" % o

    return options

def mandatory_options_set(options) :
    ret = True
    
    for m in get_mandatory_options() :
        if options[m] == None :
            print >> sys.stderr, "Error: %s must be set." % m
            ret = False
    
    return ret

def main() :
    signal.signal(signal.SIGINT, handler_sigterm)

    system = System()
    system.check_local_installation(get_required_programs())

    options = parse_args()

    if not mandatory_options_set(options) :
        sys.exit(-1)

    if False in [system.check_directories([(options['datadir'], False), (options['tempdir'], True)]), \
                 system.check_file(options['metadata'])] :
        sys.exit(-1)

    wf = WorkFlow(options['datadir'], options['tempdir'], options['metadata'])
    wf.run(options['minlength'], options['minquality'], options['winquality'])

if __name__ == '__main__' :
    main()

