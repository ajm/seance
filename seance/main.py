import sys
import os
import getopt
import logging

from os.path import splitext, join

from seance.workflow import WorkFlow
from seance.system import System


verbose_level = logging.DEBUG # logging.INFO
default_dir = './out_dir'
default_prefix = join(default_dir, 'seance')

def get_default_options() :
    global default_dir

    return {
            'input-files'       : [],
            'outdir'            : default_dir,
            'metadata'          : None,

            'denoise'           : False,
            'forwardprimer'     : None,
            'reverseprimer'     : None,
            'clipprimers'       : False,
            'miderrors'         : 0,
            'midlength'         : 5,

            'length'            : 300, 
            'quality'           : 25,
            'windowlength'      : None,
            'removeambiguous'   : False,
            'maxhomopolymer'    : 8,
            'chimeras'          : False,

            'duplicate-threshold'       : 2,
            'sample-threshold'          : 2,
            'contamination-threshold'   : 0,
            'otu-similarity'            : 0.97,
            'blast-centroids'           : False,
            'no-homopolymer-correction' : False,

            'silva-prefix'      : None,

            'output-prefix'     : default_prefix,
            'cluster-fasta'     : default_prefix + '.clusters.fasta',
            'cluster-biom'      : default_prefix + '.clusters.biom',
            'phylogeny-fasta'   : default_prefix + '.phylogeny.fasta',
            'phylogeny-tree'    : default_prefix + '.phylogeny.tree',

            'verbose'           : False
           }

def get_all_programs() :
    return ['sff2fastq', 'pagan', 'raxml', 'exonerate', 'uchime', 'blastn', 'mothur', 'cutadapt']

def get_required_programs(command, options) :
    tmp = []

    if command == 'preprocess' :
        tmp.append('uchime')

        if '.sff' in [splitext(i)[1] for i in options['input-files']] :
            tmp.append('sff2fastq')

        if options['denoise'] :
            tmp.append('mothur')

        if options['clipprimers'] :
            tmp.append('cutadapt')

    elif command == 'cluster' :
        tmp.append('pagan')

        if options['blast-centroids'] :
            tmp.append('blastn')

    elif command == 'phylogeny' :
        tmp.append('pagan')
        tmp.append('exonerate')

        if options['silva-prefix'] == None :
            tmp.append('raxml')

    elif command == 'heatmap' :
        pass

    return tmp

def get_commands() :
    return ['preprocess', 'cluster', 'phylogeny', 'heatmap']

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

    print >> sys.stderr, """Usage: %s command [OPTIONS] <sff files>

Legal commands are %s (see below for options).
%s assumes that the following programs are \ninstalled: %s.

    Preprocess options:
        -o DIR      --outdir=DIR            (default = %s)

        -d          --denoise               (default = %s)
        -p SEQ      --forwardprimer=SEQ     (default = %s)
        -r SEQ      --reverseprimer=SEQ     (default = %s)
        -k          --clipprimers           (default = %s)
        -e NUM      --miderrors=NUM         (default = %s)
        -g NUM      --midlength=NUM         (default = %s)

        -l NUM      --length=NUM            (default = %s)
        -q NUM      --quality=NUM           (default = %s)
        -w NUM      --windowlength=NUM      (default = %s)
        -x NUM      --maxhomopolymer=NUM    (default = %s)
        -n          --removeambiguous       (default = %s)
                    --chimeras              (default = %s)

    Cluster options:
        -o DIR      --outdir=DIR            (default = %s)
        -m FILE     --metadata=FILE         (default = %s)
        -a NUM      --duplicates=NUM        (default = %s)
        -b NUM      --samples=NUM           (default = %s)
        -c NUM      --contamination=NUM     (default = %s)
        -t REAL     --similarity=REAL       (default = %s)
                    --blastcentroids        (default = %s)
                    --nohomopolymer         (default = %s)
                    --output=FILEPREFIX     (default = %s)

    Phylogeny options:
        -o DIR      --outdir=DIR            (default = %s)
        -m FILE     --metadata=FILE         (default = %s)
                    --clusters=FILE         (default = %s)
        -s FILE     --silva=FILEPREFIX      (default = %s)
                    --output=FILEPREFIX     (default = %s)

    Heatmap options:
                    --biom=FILE             (default = %s)
                    --tree=FILE             (default = %s)

    Misc options:
        -v          --verbose               (default = %s)
        -h          --help
""" % (
        sys.argv[0],
        list_sentence(quote_all(bold_all(get_commands()))),
        sys.argv[0],
        list_sentence(bold_all(get_all_programs())),
        str(options['outdir']),
        str(options['denoise']),                 
        str(options['forwardprimer']),
        str(options['reverseprimer']),
        str(options['clipprimers']),
        str(options['miderrors']),
        str(options['midlength']),
        str(options['length']),
        str(options['quality']), 
        str(options['windowlength']),
        str(options['maxhomopolymer']),
        str(options['removeambiguous']), 
        str(options['chimeras']),
        str(options['outdir']),
        str(options['metadata']),
        str(options['duplicate-threshold']),
        str(options['sample-threshold']), 
        str(options['contamination-threshold']),
        str(options['otu-similarity']),
        str(options['blast-centroids']),
        str(options['no-homopolymer-correction']),
        str(options['output-prefix']),
        str(options['outdir']),
        str(options['metadata']),
        str(options['cluster-fasta']),
        str(options['silva-prefix']),
        str(options['output-prefix']),
        str(options['cluster-biom']),
        str(options['phylogeny-tree']),
        str(options['verbose'])
      )

def setup_logging(verbose) :
    global verbose_level

    log = logging.getLogger('seance')
    log.setLevel(logging.DEBUG)

    fh = logging.FileHandler('seance.log')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(levelname)s %(message)s'))

    ch = logging.StreamHandler()
    ch.setLevel(verbose_level if verbose else logging.WARNING)
    ch.setFormatter(logging.Formatter('%(levelname)s %(message)s'))

    log.addHandler(fh)
    log.addHandler(ch)

    return log

def expect_cast(parameter, argument, func) :
    try :
        return func(argument)

    except ValueError, ve :
        print >> sys.stderr, "Problem parsing argument for %s: %s\n" % (parameter, str(ve))
        usage()
        sys.exit(1)

def expect_int(parameter, argument) :
    return expect_cast(parameter, argument, int)

def expect_float(parameter, argument) :
    return expect_cast(parameter, argument, float)

def expect_iupac(parameter, argument) :
    tmp = argument.upper()
    for i in tmp :
        if i not in "TAGCRYSWKMBDHVN" :
            print >> sys.stderr, "Problem parsing argument for %s: contains illegal character '%s'\n" % (parameter, i)
            sys.exit(1)

    return tmp

def parse_args(args) :
    options = get_default_options()

    # if user sets --denoise, then we don't really need to do our own
    # length and quality check, but if they explicitly set them, do it
    tmp = {}
    tmp['quality'] = options['quality']
    tmp['length'] = options['length']

    options['quality'] = None
    options['length'] = None

    try :
        opts,args = getopt.getopt(
                        args,
                        "o:dp:r:ke:g:l:q:w:x:nm:a:b:c:t:s:vh",
                        [   "outdir=", 
                            "denoise", 
                            "forwardprimer=", 
                            "reverseprimer=", 
                            "clipprimers", 
                            "miderrors=", 
                            "midlength=", 
                            "length=",
                            "quality=",
                            "windowlength=",
                            "maxhomopolymer=",
                            "removeambiguous",
                            "metadata=",
                            "duplicates=",
                            "samples=",
                            "contamination=",
                            "similarity=",
                            "silva=",
                            "verbose",
                            "help",
                            "blastcentroids",
                            "chimeras",
                            "nohomopolymer"
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

        elif o in ('-v', '--verbose') :
            options['verbose'] = True

        elif o in ('-o', '--outdir') :
            options['outdir'] = a

        elif o in ('-m', '--metadata') :
            options['metadata'] = a

        elif o in ('-d', '--denoise') :
            options['denoise'] = True

        elif o in ('-p', '--forwardprimer') :
            options['forwardprimer'] = expect_iupac("forwardprimer", a)

        elif o in ('-r', '--reverseprimer') :
            options['reverseprimer'] = expect_iupac("reverseprimer", a)

        elif o in ('-k', '--clipprimers') :
            options['clipprimers'] = True

        elif o in ('-e', '--miderrors') :
            options['miderrors'] = expect_int("miderrors", a)

        elif o in ('-g', '--mid-length') :
            options['midlength'] = expect_int("midlength", a)

        elif o in ('-l', '--length') :
            options['length'] = expect_int("length", a)
        
        elif o in ('-q', '--quality') :
            options['quality'] = expect_int("quality", a)

        elif o in ('-w', '--windowlength') :
            options['windowlength'] = expect_int("windowlength", a)

        elif o in ('-x', '--maxhomopolymer') :
            options['maxhomopolymer'] = expect_int("maxhomopolymer", a)

        elif o in ('-n', '--removeambiguous') :
            options['removeambiguous'] = True

        elif o in ('-a', '--duplicates') :
            options['duplicate-threshold'] = expect_int("duplicates", a)

        elif o in ('-b', '--samples') :
            options['sample-threshold'] = expect_int("samples", a)

        elif o in ('-c', '--contamination') :
            options['contamination-threshold'] = expect_int("contamination", a)

        elif o in ('-t', '--similarity') :
            options['otu-similarity'] = expect_float("similarity", a)

        elif o in ('-s', '--silva') :
            options['silva'] = a

        elif o in ('--blastcentroids') :
            options['blast-centroids'] = True

        elif o in ('--chimeras') :
            options['chimeras'] = True

        elif o in ('--nohomopolymer') :
            options['no-homopolymer-correction'] = True

        else :
            assert False, "unhandled option %s" % o


    options['input-files'] = args

    # reset default quality and length if denoise is not set
    # and user has not explicitly set them
    if not options['denoise'] :
        if not options['quality'] :
            options['quality'] = tmp['quality']
        if not options['length'] :
            options['length'] = tmp['length']


    return options

def check_options(command, options) :
    system = System()

    if not system.check_directory(options['outdir'], create=True) :
        sys.exit(1)

    if command == 'preprocess' :
        if not system.check_files(options['input-files']) :
            sys.exit(1)

        if options['denoise'] and (options['forwardprimer'] is None) :
            print >> sys.stderr, "Error: forward primer must be specified with denoise flag"
            sys.exit(1)

    elif command == 'cluster' :
        if options['metadata'] is None :
            print >> sys.stderr, "Error: you must specify a metadata file"
            sys.exit(1)

        if not system.check_file(options['metadata']) :
            sys.exit(1)

        if options['sample-threshold'] <= 0 :
            print >> sys.stderr, "Error: sample-threshold must be > 0 (read %d)" % options['sample-threshold']
            sys.exit(1)

        if options['duplicate-threshold'] <= 0 :
            print >> sys.stderr, "Error: read-threshold must be > 0 (read %d)" % options['duplicate-threshold']
            sys.exit(1)

        if options['otu-similarity'] < 0.0 or options['otu-similarity'] > 1.0 :
            print >> sys.stderr, "Error: similarity must be between 0.0 and 1.0 (read %.2f)" % options['otu-similarity']
            sys.exit(1)

    elif command == 'phylogeny' :
        if not system.check_file(options['cluster-fasta']) :
            sys.exit(1)

    elif command == 'heatmap' :
        if not system.check_files([options['cluster-biom'], options['phylogeny-tree']]) :
            sys.exit(1)

def main() :
    if (len(sys.argv) < 2) or (sys.argv[1] in ('-h', '--help', 'help')) :
        usage()
        return 1

    command = sys.argv[1]
    if command not in get_commands() :
        print >> sys.stderr, "Error: unknown command '%s'" % command
        usage()
        return 1

    options = parse_args(sys.argv[2:])
    setup_logging(options['verbose'])
    check_options(command, options)

    system = System()
    system.check_local_installation(get_required_programs(command, options))
    System.tempdir(options['outdir']) # some objects need this set

    wf = WorkFlow(options)

    if command == 'preprocess' :
        return wf.preprocess()
    
    elif command == 'cluster' :
        return wf.cluster()
    
    elif command == 'phylogeny' :
        return wf.phylogeny()
    
    elif command == 'heatmap' :
        return wf.heatmap()

    return 1

if __name__ == '__main__' :
    try :
        sys.exit(main())
    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by user"

