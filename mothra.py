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

            'otu-similarity'    : 0.97,
            'otu-dup-threshold' : 100,

            'verbose'           : False
           }

def get_mandatory_options() :
    return ['datadir', 'tempdir', 'metadata']

def get_required_programs() :
    return ['sff2fastq', 'pagan', 'raxml', 'exonerate', 'uchime']

def get_commands() :
    return ['all', 'preprocess', 'phylogeny', 'otu']

def clean_up() :
    pass

def handler_sigterm(signal, frame) :
    clean_up()
    sys.exit(0)

def bold(s) :
    return "\033[1m%s\033[0m" % s

def bold_all(l) :
    return map(bold, l)

def quote(s) :
    return "'%s'" % s

def quote_all(l) :
    return map(quote, l)

def list_sentence(l) :
    if len(l) < 2 :
        return "".join(l)
    return "%s and %s" % (', '.join(l[:-1]), l[-1])

def usage() :
    options = get_default_options()

    print >> sys.stderr, """Usage: %s command [OPTIONS]

Legal commands are %s (see below for options).
%s assumes that the following programs are installed: %s.
    
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

    OTU options:
        -o      --otu-similarity      (default = %s)
        -x      --otu-dup-threshold   (default = %s)

    Misc options:
        -v      --verbose   (default = %s)
        -h      --help

""" % (sys.argv[0],                             list_sentence(quote_all(bold_all(get_commands()))), 
       sys.argv[0],                             list_sentence(bold_all(get_required_programs())),
       str(options['datadir']),                 str(options['tempdir']), 
       str(options['metadata']),                str(options['minimum-quality']), 
       str(options['window-length']),           str(options['dont-remove-nbases']), 
       str(options['compressed-length']),       str(options['mid-errors']), 
       str(options['mid-length']),              str(options['phyla-read-threshold']), 
       str(options['phyla-sample-threshold']),  str(options['otu-similarity']), 
       str(options['otu-dup-threshold']),       str(options['verbose']))

def expect_cast(parameter, argument, func) :
    try :
        return func(argument)

    except ValueError, ve :
        print >> sys.stderr, "Problem parsing argument for %s: %s\n" % (parameter, str(ve))
        usage()
        sys.exit(-1)

def expect_int(parameter, argument) :
    return expect_cast(parameter, argument, int)

def expect_float(parameter, argument) :
    return expect_cast(parameter, argument, float)

def parse_args(args) :
    options = get_default_options()

    try :
        opts,args = getopt.getopt(
                        args,
                        "d:t:m:hnq:w:l:e:g:a:b:o:x:",
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
                            "phyla-sample-threshold=",
                            "otu-similarity=",
                            "otu-dup-threshold="
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

        elif o in ('-o', '--otu-similarity') :
            options['otu-similarity'] = expect_float("otu-similarity", a)

        elif o in ('-x', '--otu-dup-threshold') :
            options['otu-dup-threshold'] = expect_int("otu-dup-threshold", a)

        elif o in ('-v', '--verbose') :
            options['verbose'] = True

        else :
            assert False, "unhandled option %s" % o


    check_options(options)

    return options

def mandatory_options_set(options) :
    ret = True
    
    for m in get_mandatory_options() :
        if options[m] == None :
            print >> sys.stderr, "Error: %s must be set." % m
            ret = False
    
    return ret

def check_options(options) :
    system = System()

    if not mandatory_options_set(options) :
        sys.exit(-1)

    if False in [system.check_directories([(options['datadir'], False), (options['tempdir'], True)]), \
                 system.check_file(options['metadata'])] :
        sys.exit(-1)

    if options['phyla-sample-threshold'] <= 0 :
        print >> sys.stderr, "Error: phyla-sample-threshold must be > 0 (read %d)" % options['phyla-sample-threshold']
        sys.exit(-1)

    if options['phyla-read-threshold'] <= 0 :
        print >> sys.stderr, "Error: phyla-read-threshold must be > 0 (read %d)" % options['phyla-read-threshold']
        sys.exit(-1)

    if options['otu-similarity'] < 0.0 or options['otu-similarity'] > 1.0 :
        print >> sys.stderr, "Error: otu-similarity must be between 0.0 and 1.0 (read %.2f)" % options['otu-similarity']
        sys.exit(-1)

    if options['otu-dup-threshold'] <= 0 :
        print >> sys.stderr, "Error: otu-dup-threshold must be > 0 (read %d)" % options['otu-dup-threshold']
        sys.exit(-1)

def main() :
    signal.signal(signal.SIGINT, handler_sigterm)

    if len(sys.argv) < 2 :
        usage()
        return -1

    command = sys.argv[1]
    if command not in get_commands() :
        print >> sys.stderr, "Error: unknown command '%s'" % command
        usage()
        return -1

    options = parse_args(sys.argv[2:])

    system = System()
    system.check_local_installation(get_required_programs())
    System.tempdir(options['tempdir']) # some objects need this set

    wf = WorkFlow(options)
    if command in ('preprocess', 'all') :
        wf.preprocess()

    if command in ('otu', 'all') :
        wf.otu_phylogeny()

    if command in ('phylogeny', 'all') :
        wf.phylogeny()

    return 0

if __name__ == '__main__' :
    sys.exit(main())

