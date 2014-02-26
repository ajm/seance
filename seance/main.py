import os
import getopt
import logging

from sys import argv, exit, stderr
from os.path import splitext, join

from seance.workflow import WorkFlow
from seance.system import System


__program__ = "seance"
verbose_level = logging.DEBUG # logging.INFO
default_dir = './out'
default_prefix = 'seance'

def get_default_options(fillin=False) :
    global default_dir

    tmp = {
            'input-files'       : [],
            'outdir'            : default_dir,
            'metadata'          : None,

            'denoise'           : False,
            'forwardprimer'     : None,
            'reverseprimer'     : None,
            'clipprimers'       : False,
            'miderrors'         : 0,
            'midlength'         : 5,

            'length'            : 250, 
            'quality-method'    : 'min', # 'average', 'window', 'none'
            'quality'           : 25,
            'windowlength'      : 10,
            'removeambiguous'   : True,
            'maxhomopolymer'    : 8,
            'chimeras'          : False,

            'total-duplicate-threshold'     : 1,
            'sample-threshold'              : 2,
            'duplicate-threshold'           : 2,
            'otu-similarity'                : 0.99,
            'blast-centroids'               : False,
            'merge-blast-hits'              : False,
            'no-homopolymer-correction'     : False,

            'output-prefix'     : default_prefix,
            'cluster-fasta'     : None,
            'cluster-biom'      : None,
            'phylogeny-fasta'   : None,
            'phylogeny-tree'    : None,
            'phylogeny-xml'     : None,

            'silva-fasta'       : None,
            'silva-tree'        : None,

            'heatmap-no-tree'   : False,
            'wasabi-url'        : 'http://wasabi2.biocenter.helsinki.fi:8000',
            'wasabi-user'       : None,

            'verbose'           : False
           }

    if fillin :
        apply_output_prefix(tmp)

    return tmp

def apply_output_prefix(d, command='all') :
    tmp = d['output-prefix']
    tmp = join(d['outdir'], tmp)

    if not d['cluster-fasta'] :
        d['cluster-fasta']   = tmp + '.cluster.fasta'

    if not d['cluster-biom'] :
        d['cluster-biom']    = tmp + '.cluster.biom'
    
    if not d['phylogeny-fasta'] :
        d['phylogeny-fasta'] = tmp + '.phylogeny.fasta'

    if not d['phylogeny-tree'] :
        d['phylogeny-tree']  = tmp + '.phylogeny.tree'

    if not d['phylogeny-xml'] :
        d['phylogeny-xml'] = tmp + '.phylogeny.xml'

def get_all_programs() :
    return ['sff2fastq', 'pagan', 'raxml', 'bppphysamp', 'exonerate', 'uchime', 'blastn', 'ampliconnoise']

def get_required_programs(command, options) :
    tmp = []

    if command == 'preprocess' :
        if options['chimeras'] :
            tmp.append('uchime')

        if '.sff' in [splitext(i)[1] for i in options['input-files']] :
            tmp.append('sff2fastq')

        if options['denoise'] :
            tmp.append('PyroDist')
            tmp.append('FCluster')
            tmp.append('PyroNoise')

        #if options['clipprimers'] :
        #    tmp.append('cutadapt')

    elif command == 'cluster' :
        tmp.append('pagan')

        if options['blast-centroids'] :
            tmp.append('blastn')

    elif command == 'phylogeny' :
        tmp.append('pagan')
        tmp.append('exonerate')

        if not options['silva-fasta'] :
            tmp.append('raxml')

    return tmp

def get_commands() :
    return ['preprocess', 'cluster', 'phylogeny', 'heatmap', 'wasabi']

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

def usage(command='all') :
    global __program__

    options = get_default_options(True)
    outdir = options['outdir']

    print >> stderr, """Usage: %s <command> [OPTIONS] <sff files>
Legal commands are %s (see below for options).
%s assumes that the following programs are installed: %s.\n""" % \
               (bold(__program__),
                list_sentence(quote_all(bold_all(get_commands()))), 
                bold(__program__),
                list_sentence(bold_all(get_all_programs())))
    
    print >> stderr, """    Common options:
        -o DIR      --outdir=DIR            (default = %s)
        -v          --verbose\n""" % \
                (str(options['outdir']))

    if command in ('preprocess','all') :
        print >> stderr, """    Preprocess options:
        -p SEQ      --forwardprimer=SEQ     (default = %s)
        -r SEQ      --reverseprimer=SEQ     (default = %s)
        -k          --clipprimers           (default = %s)

        -e NUM      --miderrors=NUM         (default = %s)
        -g NUM      --midlength=NUM         (default = %s)

        -l NUM      --length=NUM            (default = %s)
        -x NUM      --maxhomopolymer=NUM    (default = %s)
        -n          --keepambiguous         (default = %s)

                    --qualitymethod=X       (default = %s, options = (none, min, average, window))
        -q NUM      --quality=NUM           (default = %s)
        -w NUM      --windowlength=NUM      (default = %s)

        -d          --denoise               (default = %s)
                    --chimeras              (default = %s)\n""" % \
               (str(options['forwardprimer']),
                str(options['reverseprimer']),
                str(options['clipprimers']),
                str(options['miderrors']),
                str(options['midlength']),
                str(options['length']),
                str(options['maxhomopolymer']),
                str(not options['removeambiguous']),
                options['quality-method'],
                str(options['quality']), 
                str(options['windowlength']),
                str(options['denoise']),
                str(options['chimeras']))

    if command in ('cluster','all') :
        print >> stderr, """    Cluster options:
        -m FILE     --metadata=FILE         (default = %s)

        -a NUM      --totalduplicates=NUM   (default = %s)
        -b NUM      --samples=NUM           (default = %s)
        -c NUM      --duplicates=NUM        (default = %s)
        
        -t REAL     --similarity=REAL       (default = %s)
        
                    --blastcentroids        (default = %s)
                    --mergeblasthits        (default = %s)
                    --nohomopolymer         (default = %s)
                    --output=FILEPREFIX     (default = %s{.cluster.fasta, 
                                                          .cluster.biom})\n""" % \
               (options['metadata'],
                str(options['total-duplicate-threshold']),
                str(options['sample-threshold']), 
                str(options['duplicate-threshold']),
                str(options['otu-similarity']),
                str(options['blast-centroids']),
                str(options['merge-blast-hits']),
                str(options['no-homopolymer-correction']),
                options['output-prefix'])

    if command in ('phylogeny','all') :
        print >> stderr, """    Phylogeny options:
                    --clusters=FILE         (default = %s)
        -s FILE     --silva=FILEPREFIX      (default = %s, expects FILEPREFIX{.fasta, .tree})
                    --output=FILEPREFIX     (default = %s{.phylogeny.fasta, 
                                                          .phylogeny.tree, 
                                                          .phylogeny.xml})\n""" % \
               (options['cluster-fasta'],
                options['silva-fasta'],
                options['output-prefix'])

    if command in ('heatmap','all') :
        print >> stderr, """    Heatmap options:
                    --biom=FILE             (default = %s)
                    --tree=FILE             (default = %s)
                    --notree
                    --output=FILE           (default = %s.pdf)\n""" % \
               (options['cluster-biom'],
                options['phylogeny-tree'],
                options['output-prefix'])

    if command in ('wasabi','all') :
        print >> stderr, """    Wasabi options:
                    --xml=FILE              (default = %s)
                    --user=USER
                    --url=URL               (default = %s)\n""" % \
               (options['phylogeny-xml'],
                options['wasabi-user'],
                options['wasabi-url'])

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
        print >> stderr, "Problem parsing argument for %s: %s\n" % (parameter, str(ve))
        usage()
        exit(1)

def expect_int(parameter, argument) :
    return expect_cast(parameter, argument, int)

def expect_float(parameter, argument) :
    return expect_cast(parameter, argument, float)

def expect_iupac(parameter, argument) :
    tmp = argument.upper()
    for i in tmp :
        if i not in "TAGCRYSWKMBDHVN" :
            print >> stderr, "Problem parsing argument for %s: contains illegal character '%s'\n" % (parameter, i)
            exit(1)

    return tmp

def parse_args(command, args) :
    options = get_default_options()

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
                            "qualitymethod=",
                            "quality=",
                            "windowlength=",
                            "maxhomopolymer=",
                            "keepambiguous",
                            "metadata=",
                            "totalduplicates=",
                            "samples=",
                            "duplicates=",
                            "similarity=",
                            "silva=",
                            "verbose",
                            "help",
                            "blastcentroids",
                            "mergeblasthits",
                            "chimeras",
                            "nohomopolymer",
                            "output=",
                            "biom=",
                            "tree=",
                            "notree",
                            "clusters=",
                            "xml=",
                            "url=",
                            "user="
                        ]
                    )

    except getopt.GetoptError, err :
        print >> stderr, str(err) + "\n"
        usage()
        exit(1)

    for o,a in opts :
        if o in ('-h', '--help') :
            usage(command)
            exit(0)

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
        
        elif o in ('--qualitymethod') :
            methods = ['none', 'min', 'average', 'window']
            if a in methods :
                options['quality-method'] = a
            else :
                print >> stderr, "ERROR %s is not a valid quality method (valid options: %s)" % \
                        (bold(a), list_sentence(bold_all(methods)))
                exit(1)

        elif o in ('-q', '--quality') :
            options['quality'] = expect_int("quality", a)

        elif o in ('-w', '--windowlength') :
            options['windowlength'] = expect_int("windowlength", a)

        elif o in ('-x', '--maxhomopolymer') :
            options['maxhomopolymer'] = expect_int("maxhomopolymer", a)

        elif o in ('-n', '--keepambiguous') :
            options['removeambiguous'] = False

        elif o in ('-a', '--totalduplicates') :
            options['total-duplicate-threshold'] = expect_int("total-duplicates", a)

        elif o in ('-b', '--samples') :
            options['sample-threshold'] = expect_int("samples", a)

        elif o in ('-c', '--duplicates') :
            options['duplicate-threshold'] = expect_int("duplicates", a)

        elif o in ('-t', '--similarity') :
            options['otu-similarity'] = expect_float("similarity", a)

        elif o in ('-s', '--silva') :
            options['silva-fasta'] = a + '.fasta'
            options['silva-tree'] = a + '.tree'

        elif o in ('--blastcentroids') :
            options['blast-centroids'] = True

        elif o in ('--mergeblasthits') :
            options['merge-blast-hits'] = True

        elif o in ('--chimeras') :
            options['chimeras'] = True

        elif o in ('--nohomopolymer') :
            options['no-homopolymer-correction'] = True

        elif o in ('--clusters') :
            options['cluster-fasta'] = a

        elif o in ('--biom') :
            options['cluster-biom'] = a

        elif o in ('--tree') :
            options['phylogeny-tree'] = a

        elif o in ('--output') :
            options['output-prefix'] = a

        elif o in ('--xml') :
            options['phylogeny-xml'] = a

        elif o in ('--url') :
            options['wasabi-url'] = a

        elif o in ('--user') :
            options['wasabi-user'] = a

        elif o in ('--notree') :
            options['heatmap-no-tree'] = True

        else :
            assert False, "unhandled option %s" % o


    options['input-files'] = args

    return options

def check_options(command, options) :
    system = System()

    apply_output_prefix(options, command)

    if not system.check_directory(options['outdir'], create=True) :
        exit(1)

    if command == 'preprocess' :
        if not system.check_files(options['input-files']) :
            exit(1)

    elif command == 'cluster' :
        if options['metadata'] is None :
            print >> stderr, "Error: you must specify a metadata file"
            exit(1)

        if not system.check_file(options['metadata']) :
            exit(1)

        for i in ('duplicate-threshold', 'total-duplicate-threshold', 'sample-threshold') :
            if options[i] <= 0 :
                print >> stderr, "Error: %s must be > 0 (read %d)" % (i, options[i])
                exit(1)

        # i think pagan's sensitivity limit is ~80%
        if options['otu-similarity'] < 0.8 or options['otu-similarity'] > 1.0 :
            print >> stderr, "Error: similarity must be between 0.8 and 1.0 (read %.2f)" % options['otu-similarity']
            exit(1)

    elif command == 'phylogeny' :
        if not system.check_file(options['cluster-fasta']) :
            exit(1)

        if options['silva-fasta'] :
            if not system.check_files([options['silva-fasta'], options['silva-tree']]) :
                exit(1)

    elif command == 'heatmap' :
        if not system.check_files([options['cluster-biom']]) :
            exit(1)

    elif command == 'wasabi' :
        if not options['wasabi-user'] :
            print >> stderr, "Error: you must specify your wasabi username!\n"
            exit(1)

        if not system.check_files([options['phylogeny-xml']]) :
            exit(1)

def main() :
    if (len(argv) < 2) or (argv[1] in ('-h', '--help', 'help')) :
        usage()
        return 1

    command = argv[1]
    if command not in get_commands() :
        print >> stderr, "Error: unknown command '%s'" % command
        usage()
        return 1

    options = parse_args(command, argv[2:])
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

    elif command == 'wasabi' :
        return wf.wasabi()

    return 1

if __name__ == '__main__' :
    try :
        exit(main())
    except KeyboardInterrupt :
        print >> stderr, "Killed by user"

