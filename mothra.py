import sys
import os
import getopt
import signal

from metagenomics.workflow import WorkFlow
from metagenomics.system import System


def get_default_options() :
    return {
            'datadir'           : None,
            'tempdir'           : None,
            'metadata'          : None,
            
            'compressed-length' : 300, 
            'minimum-quality'   : 20,
            'window-length'     : None,
            'dont-remove-nbases': False,
            'mid-errors'        : 0,
            'mid-length'        : 5,

            'phyla-read-threshold'   : 10,
            'phyla-sample-threshold' : 2,

            'verbose'           : False
           }

def get_mandatory_options() :
    return ['datadir', 'tempdir', 'metadata']

def get_required_programs() :
    return ['sff2fastq', 'pagan', 'raxml', 'exonerate', 'uchime']

def clean_up() :
    pass

def handler_sigterm(signal, frame) :
    clean_up()
    sys.exit(0)

def usage() :
    options = get_default_options()

    print >> sys.stderr, """Usage: %s command [OPTIONS]
    Mandatory:
        -d      --datadir   (default = %s)
        -t      --tempdir   (default = %s)
        -m      --metadata  (default = %s)

    Preprocess options:
        -q      --minimum-quality     (default = %s)
        -w      --window-length       (default = %s)
        -n      --dont-remove-nbases  (default = %s)
        -l      --compressed-length   (default = %s)
        -e      --mid-errors          (default = %s)
        -g      --mid-length          (default = %s)

    Phylogeny options:
        -a      --phyla-read-threshold      (default = %s)
        -b      --phyla-sample-threshold    (default = %s)

    Misc options:
        -v      --verbose             (default = %s)
        -h      --help
""" % (sys.argv[0], 
       str(options['datadir']),                 str(options['tempdir']), 
       str(options['metadata']),                str(options['minimum-quality']), 
       str(options['window-length']),           str(options['dont-remove-nbases']), 
       str(options['compressed-length']),       str(options['mid-errors']), 
       str(options['mid-length']),              str(options['phyla-read-threshold']), 
       str(options['phyla-sample-threshold']),  str(options['verbose']))

def expect_int(parameter, argument) :
    try :
        return int(argument)

    except ValueError, ve :
        print >> sys.stderr, "Problem parsing argument for %s: %s\n" % (parameter, str(ve))
        usage()
        sys.exit(-1)

def parse_args(args) :
    options = get_default_options()

    try :
        opts,args = getopt.getopt(
                        args,
                        "d:t:m:hnq:w:l:e:g:a:b:",
                        [   "help", 
                            "verbose", 
                            "datadir=", 
                            "tempdir=", 
                            "metadata=", 
                            "minimum-quality=", 
                            "window-length=", 
                            "remove-nbases", 
                            "compressed-length=",
                            "mid-errors=",
                            "mid-length=",
                            "phyla-read-threshold=",
                            "phyla-sample-threshold="
                        ]
                    )

    except getopt.GetoptError, err :
        print >> sys.stderr, str(err) + "\n"
        usage()
        sys.exit(-1)

    for o,a in opts :
        if o in ('-h', '--help') :
            usage()
            sys.exit(0)

        elif o in ('-d', '--datadir') :
            options['datadir'] = a

        elif o in ('-t', '--tempdir') :
            options['tempdir'] = a

        elif o in ('-m', '--metadata') :
            options['metadata'] = a

        elif o in ('-q', '--minimum-quality') :
            options['minquality'] = expect_int("minquality", a)

        elif o in ('-w', '--window-length') :
            options['winquality'] = expect_int("winquality", a)

        elif o in ('-n', '--dont-remove-nbases') :
            options['dont-remove-nbases'] = True

        elif o in ('-c', '--compress-length') :
            options['compress-length'] = expect_int("compressed-length", a)
        
        elif o in ('-e', '--mid-errors') :
            options['mid-errors'] = expect_int("mid-errors", a)

        elif o in ('-g', '--mid-length') :
            options['mid-length'] = expect_int("mid-length", a)

        elif o in ('-a', '--phyla-read-threshold') :
            options['phyla-read-threshold'] = expect_int("phyla-read-threshold", a)

        elif o in ('-b', '--phyla-sample-threshold') :
            options['phyla-sample-threshold'] = expect_int("phyla-sample-threshold", a)

        elif o in ('-v', '--verbose') :
            options['verbose'] = True

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

    if len(sys.argv) < 2 :
        usage()
        return -1

    command = sys.argv[1]
    if command not in ['all', 'summary', 'preprocess', 'phylogeny'] :
        print >> sys.stderr, "Error: unknown command '%s'" % command
        usage()
        return -1

    options = parse_args(sys.argv[2:])

    if not mandatory_options_set(options) :
        return -1

    if False in [system.check_directories([(options['datadir'], False), (options['tempdir'], True)]), \
                 system.check_file(options['metadata'])] :
        return -1

    System.tempdir(options['tempdir'])


    wf = WorkFlow(options)
    if command == 'preprocess' or command == 'all' :
        wf.preprocess()

    if command == 'phylogeny' or command == 'all' :
        wf.phylogeny()

    return 0

if __name__ == '__main__' :
    sys.exit(main())

