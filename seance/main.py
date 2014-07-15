#Seance: a bioinformatics pipeline for reference-based phylogenetic analysis
#Copyright (C) 2014  Alan Medlar (amedlar@gmail.com)
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import getopt
import logging

from sys import argv, exit, stderr
from os.path import splitext, join

from seance.workflow import WorkFlow
from seance.system import System


__program__ = "seance"
verbose_level = logging.INFO
default_dir = './out'
default_prefix = 'seance'

def get_default_options(command, fillin=False) :
    global default_dir

    tmp = {
            'input-files'       : [],
            'outdir'            : default_dir,
            'metadata'          : None,

            'denoise'           : False,
            'forwardprimer'     : None,
            'reverseprimer'     : None,
            'clipprimers'       : False,
            'miderrors'         : 1,
            'midlength'         : 5,
            'primererrors'      : 2,

            'length'            : 250, 
            'quality-method'    : 'min', # 'average', 'window', 'none'
            'quality'           : 20,
            'windowlength'      : 50,
            'removeambiguous'   : True,
            'maxhomopolymer'    : 8,
            'chimeras'          : False,

            'total-duplicate-threshold'     : 1,
            'sample-threshold'              : 1,
            'duplicate-threshold'           : 2,
            'otu-similarity'                : 0.99,
            'merge-blast-hits'              : False,
            'no-homopolymer-correction'     : False,

            'summary-file'      : None,

            'labels'            : None,
            'label-missing'     : False,

            'delimiter'         : '\t',
            'prefix'            : default_prefix,
            'cluster-fasta'     : None,
            'cluster-biom'      : None,
            'phylogeny-fasta'   : None,
            'phylogeny-tree'    : None,
            'phylogeny-xml'     : None,

            'silva-fasta'       : None,
            'silva-tree'        : None,
            'denovo'            : False,

            'subset'            : None,

            'heatmap-no-tree'   : False,
            'heatmap-pdf'       : 'seance.pdf',
            'heatmap-out-tree'  : None,
            'heatmap-flip-tree' : False,
            'heatmap-ladderise' : False,
            'heatmap-scale'     : 0.05,
            'heatmap-tree-height' : 20,
            'heatmap-label-clip'  : [],
            'heatmap-label-tokens' : -1,

            'wasabi-url'        : 'http://wasabi2.biocenter.helsinki.fi:8000',
            'wasabi-user'       : None,

            'verbose'           : False
           }

    if fillin :
        apply_prefix(tmp)

    if command == 'label' :
        tmp['label'] = 'taxonomy'

    return tmp

def apply_prefix(d, command='all') :
    tmp = d['prefix']
    tmp = join(d['outdir'], tmp)

    d['summary-file'] = join(d['outdir'], 'summary.csv')

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
    return ['sff2fastq', 'pagan', 'raxml', 'bppphysamp', 'exonerate', 'uchime', 'blastn', 'PyroDist', 'FCluster', 'PyroNoise']

def test_system(command=None, options=None, exit_on_failure=False, output=False) :
    binaries = { 
        'preprocess' : {
            '*'         : ['sff2fastq'],
            'chimeras'  : ['uchime'],
            'denoise'   : ['PyroDist', 'FCluster', 'PyroNoise'],
        },
        'cluster' : {
            '*'         : ['pagan'],
            'labels'    : ['blastn']
        },
        'label' : {
            '*'         : ['blastn']
        },
        'summary' : {},
        'showcounts' : {},
        'showlabels' : {},
        'phylogeny' : {
            '*'         : ['pagan', 'exonerate', 'bppphysamp'],
            'denovo'    : ['raxml']
        },
        'heatmap' : {},
        'wasabi' : {},
        'test' : {}
    }

    def print_out(s) :
        if output :
            print s
        else :
            if 'not found' in s :
                print >> stderr, s.lstrip()

    fail = False

    for b in [command] if command else binaries :
        if binaries[b] :
            print_out("checking system for %s command dependancies :" % bold(b))
            for o in binaries[b] :
                if o == '*' or not options or options[o] :
                    for p in binaries[b][o] :
                        installed = System.is_installed(p)
                        print_out("    %s %s%s" % (p, "" if o == "*" else "(needed for --%s) " % o, 
                            bold_green("found.") if installed else bold_red("not found!")))
                    
                    if not installed :
                        fail = True    
            print_out("")

    # for some reason using pip to install cairo always fails
    # so check for it here instead
#    print "checking for python modules :"
#    try :
#        import cairo
#        print "    pycairo (needed by 'heatmap' command) " + bold_green("found.")
#    except ImportError :
#        print "    pycairo (needed by 'heatmap' command) " + bold_red("not found!")
#        fail = True
#    print ""

    if exit_on_failure and fail :
        exit(1)

def get_commands() :
    return ['test', 'preprocess', 'summary', 'cluster', 'label', 'showcounts', 'showlabels', 'phylogeny', 'heatmap', 'wasabi']

def bold_green(s) :
    return "\033[32m%s\033[0m" % s

def bold_red(s) :
    return "\033[31m%s\033[0m" % s

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

    options = get_default_options('all', True)
    outdir = options['outdir']

    print >> stderr, """Usage: %s <command> [OPTIONS] <sff files>
Legal commands are %s (see below for options).
%s assumes that the following programs are installed: %s.\n""" % \
               (bold(__program__),
                list_sentence(quote_all(bold_all(get_commands()))), 
                bold(__program__),
                list_sentence(bold_all(get_all_programs())))
    
    print >> stderr, """    Common options:
        -o DIR          --outdir=DIR            (default = %s)
        -p FILEPREFIX   --prefix=FILEPREFIX     (default = %s, overrided by biom,tree,clusters,xml)
        -v              --verbose\n""" % \
                (options['outdir'], options['prefix'])

    if command in ('preprocess','all') :
        print >> stderr, """    Preprocess options:
        -f SEQ          --forwardprimer=SEQ     (default = %s)
        -r SEQ          --reverseprimer=SEQ     (default = %s)
        -k              --clipprimers           (default = %s)
                        --primererrors=NUM      (default = %s)

        -e NUM          --miderrors=NUM         (default = %s)
        -g NUM          --midlength=NUM         (default = %s)

        -l NUM          --length=NUM            (default = %s)
        -x NUM          --maxhomopolymer=NUM    (default = %s)
        -n              --keepambiguous         (default = %s)

                        --qualmethod=X          (default = %s, options = (none, min, average, window))
        -q NUM          --quality=NUM           (default = %s)
        -w NUM          --windowlength=NUM      (default = %s)

        -d              --denoise               (default = %s)
                        --chimeras              (default = %s)\n""" % \
               (str(options['forwardprimer']),
                str(options['reverseprimer']),
                str(options['clipprimers']),
                str(options['primererrors']),
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
        -m FILE         --metadata=FILE         (default = %s)

        -a NUM          --totalduplicates=NUM   (default = %s)
        -b NUM          --samples=NUM           (default = %s)
        -c NUM          --duplicates=NUM        (default = %s)
        
        -t REAL         --similarity=REAL       (default = %s)
        
                        --labels=X              (default = none, options = (none, blast, taxonomy))
                        --mergeclusters         (default = %s)
                        --nohomopolymer         (default = %s)\n""" % \
               (options['metadata'],
                str(options['total-duplicate-threshold']),
                str(options['sample-threshold']), 
                str(options['duplicate-threshold']),
                str(options['otu-similarity']),
                str(options['merge-blast-hits']),
                str(options['no-homopolymer-correction']))

    if command in ('label', 'all') :
        print >> stderr, """    Label options:
                        --missing               (only fetch missing labels)
                        --labels=X              (default = taxonomy, options = (none, blast, taxonomy))\n"""

    if command in ('showcounts', 'all') :
        print >> stderr, """    Count options:
                        --delimiter=STR         (default = '\\t')\n"""

    if command in ('phylogeny','all') :
        print >> stderr, """    Phylogeny options:
                        --refalignment=FILE     (expects fasta)
                        --reftree=FILE          (expects newick)
                        --denovo                (don't perform phylogenetic placement, just use RaXML)
                        --subset=STR            (subset clusters for phylogenetic placement, e.g. --subset 'nematoda')
                        --clusters=FILE         (default = %s)\n""" % \
               (options['cluster-fasta'])

    if command in ('heatmap','all') :
        print >> stderr, """    Heatmap options:
                        --biom=FILE             (default = %s)
                        --tree=FILE             (default = %s)
                        --notree
                        --subset=STR            (subset samples for inclusion into heatmap)
                        --outtree=FILE          (output tree used in separate newick file)
                        --fliptree              (flip tree upside down)
                        --ladderise             (ladderise the tree)
                        --labelclip=X           (comma separated values, e.g.: nematoda,arthropoda)
                        --labeltokens=NUM       (need last NUM tokens in taxonomical labelling)
                        --scale=FLOAT           (set length of tree scale, default = %.2f)
                        --height=INT            (set height of tree in heatmap blocks, default = %d)
                        --output=FILE           (default = %s)\n""" % \
               (options['cluster-biom'],
                options['phylogeny-tree'],
                options['heatmap-scale'],
                options['heatmap-tree-height'],
                options['heatmap-pdf'])

    if command in ('wasabi','all') :
        print >> stderr, """    Wasabi options:
                        --xml=FILE              (default = %s)
                        --user=USER
                        --url=URL               (default = %s)\n""" % \
               (options['phylogeny-xml'],
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
    options = get_default_options(command)

    try :
        opts,args = getopt.getopt(
                        args,
                        "o:dp:f:r:ke:g:l:q:w:x:nm:a:b:c:t:vh",
                        [   "outdir=",
                            "prefix=",
                            "denoise", 
                            "forwardprimer=", 
                            "reverseprimer=", 
                            "clipprimers", 
                            "miderrors=", 
                            "midlength=", 
                            "length=",
                            "qualmethod=",
                            "quality=",
                            "windowlength=",
                            "maxhomopolymer=",
                            "keepambiguous",
                            "metadata=",
                            "totalduplicates=",
                            "samples=",
                            "duplicates=",
                            "similarity=",
                            "refalignment=",
                            "reftree=",
                            "output=",
                            "verbose",
                            "help",
                            "labels=",
                            "mergeclusters",
                            "chimeras",
                            "nohomopolymer",
                            "output=",
                            "biom=",
                            "tree=",
                            "notree",
                            "clusters=",
                            "xml=",
                            "url=",
                            "user=",
                            "subset=",
                            "outtree=",
                            "fliptree",
                            "treescale=",
                            "treeheight=",
                            "labelclip=",
                            "labeltokens=",
                            "denovo",
                            "delimiter=",
                            "missing",
                            "primererrors=",
                            "ladderise"
                        ]
                    )

    except getopt.GetoptError, err :
        print >> stderr, str(err) + "\n"
        usage(command)
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

        elif o in ('-f', '--forwardprimer') :
            options['forwardprimer'] = expect_iupac("forwardprimer", a)

        elif o in ('-r', '--reverseprimer') :
            options['reverseprimer'] = expect_iupac("reverseprimer", a)

        elif o in ('-k', '--clipprimers') :
            options['clipprimers'] = True

        elif o in ('-e', '--miderrors') :
            options['miderrors'] = expect_int("miderrors", a)

        elif o in ('-g', '--midlength') :
            options['midlength'] = expect_int("midlength", a)

        elif o in ('-l', '--length') :
            options['length'] = expect_int("length", a)
        
        elif o in ('--qualmethod',) :
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

        elif o in ('--refalignment',) :
            options['silva-fasta'] = a

        elif o in ('--reftree',) :
            options['silva-tree'] = a

        elif o in ('--labels',) :
            methods = ['none', 'blast', 'taxonomy']
            if a in methods :
                if a == 'none' :
                    options['labels'] = None
                else :
                    options['labels'] = a
            else :
                print >> stderr, "ERROR %s is not a valid labelling method (valid options: %s)" % \
                        (bold(a), list_sentence(bold_all(methods)))
                exit(1)

        elif o in ('--mergeclusters',) :
            options['merge-blast-hits'] = True

        elif o in ('--chimeras',) :
            options['chimeras'] = True

        elif o in ('--nohomopolymer',) :
            options['no-homopolymer-correction'] = True

        elif o in ('--clusters',) :
            options['cluster-fasta'] = a

        elif o in ('--biom',) :
            options['cluster-biom'] = a

        elif o in ('--tree',) :
            options['phylogeny-tree'] = a

        elif o in ('-p', '--prefix',) :
            options['prefix'] = a

        elif o in ('--output',) :
            options['heatmap-pdf'] = a

        elif o in ('--xml',) :
            options['phylogeny-xml'] = a

        elif o in ('--url',) :
            options['wasabi-url'] = a

        elif o in ('--user',) :
            options['wasabi-user'] = a

        elif o in ('--notree',) :
            options['heatmap-no-tree'] = True

        elif o in ('--outtree',) :
            options['heatmap-out-tree'] = a

        elif o in ('--subset',) :
            options['subset'] = a

        elif o in ('--fliptree',) :
            options['heatmap-flip-tree'] = True

        elif o in ('--treescale',) :
            options['heatmap-scale'] = expect_float("treescale", a)

        elif o in ('--treeheight',) :
            options['heatmap-tree-height'] = expect_int("treeheight", a)

        elif o in ('--labelclip',) :
            options['heatmap-label-clip'] = [ i.lower() for i in a.split(',') ]

        elif o in ('--labeltokens',) :
            options['heatmap-label-tokens'] = expect_int("labeltokens", a)
        
        elif o in ('--denovo',) :
            options['denovo'] = True
    
        elif o in ('--delimiter',) :
            if not a :
                options['delimiter'] = " "
            else :
                options['delimiter'] = a 

        elif o in ('--missing',) :
            options['label-missing'] = True

        elif o in ('--primererrors',) :
            options['primererrors'] = expect_int("primererrors", a)

        elif o in ('--ladderise',) :
            options['heatmap-ladderise'] = True

        else :
            assert False, "unhandled option %s" % o


    options['input-files'] = args

    return options

def check_options(command, options, log) :
    system = System()

    apply_prefix(options, command)

#    for i in options :
#        print i, options[i]

    if (command != 'preprocess') and (not system.check_directory(options['outdir'])) :
        exit(1)

    if command == 'preprocess' :
        if not system.check_directory(options['outdir'], create=True) :
            exit(1)

        if not system.check_files(options['input-files']) :
            exit(1)

        if options['denoise'] :
            if not options['forwardprimer'] :
                log.error("for denoising you must specify the forward primer!")
                exit(1)

    elif command == 'cluster' :
        #if options['metadata'] is None :
        #    print >> stderr, "Error: you must specify a metadata file"
        #    exit(1)

        if (options['metadata'] is not None) and (not system.check_file(options['metadata'])) :
            exit(1)

        for i in ('duplicate-threshold', 'total-duplicate-threshold', 'sample-threshold') :
            if options[i] <= 0 :
                log.error("%s must be > 0 (read %d)" % (i, options[i]))
                exit(1)

        # i think pagan's sensitivity limit is ~80%
        if options['otu-similarity'] < 0.8 or options['otu-similarity'] > 1.0 :
            log.error("similarity must be between 0.8 and 1.0 (read %.2f)" % options['otu-similarity'])
            exit(1)

    elif command == 'summary' :
        if not system.check_file(options['summary-file']) :
            log.error("could not find %s, did you run the 'preprocess' command yet?" % options['summary-file'])
            exit(1)

    elif command == 'label' :
        fasta_check = system.check_file(options['cluster-fasta'])
        biom_check  = system.check_file(options['cluster-biom'])

        def filename_or() :
            if fasta_check ^ biom_check :
                return "'%s'" % (options['cluster-fasta'] if not fasta_check else options['cluster-biom'])
            else :
                return "'%s' or '%s'" % (options['cluster-fasta'], options['cluster-biom'])

        if not fasta_check or not biom_check :
            log.error("could not find %s, did you run the 'cluster' command yet?" % filename_or())
            exit(1)

        if not options['label-centroids'] :
            log.error("you must specify a labelling method")
            exit(1)

    elif command == 'showcounts' or command == 'showlabels' :
        if not system.check_file(options['cluster-biom']) :
            log.error("could not find '%s', did you run the 'cluster' command yet?" % options['cluster-biom'])
            exit(1)

    elif command == 'phylogeny' :
        if not system.check_file(options['cluster-fasta']) :
            exit(1)

        if not options['denovo'] :
            if not options['silva-fasta'] or not options['silva-tree'] :
                log.error("you must either specify the location of the reference alignment and phylogeny or else use the --denovo option")
                exit(1)
            else :
                if not system.check_files([options['silva-fasta'], options['silva-tree']]) :
                    exit(1)

    elif command == 'heatmap' :
        if not system.check_files([options['cluster-biom']]) :
            exit(1)

    elif command == 'wasabi' :
        if not options['wasabi-user'] :
            log.error("you must specify your wasabi username!")
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

    if command == 'test' :
        test_system(output=True)
        return 0

    options = parse_args(command, argv[2:])
    log = setup_logging(options['verbose'])
    check_options(command, options, log)

    test_system(command, options, exit_on_failure=True)
    System.tempdir(options['outdir']) # some objects need this set

    wf = WorkFlow(options)

    if command == 'preprocess' :
        return wf.preprocess()
 
    elif command == 'summary' :
        return wf.summary()

    elif command == 'cluster' :
        return wf.cluster()

    elif command == 'label' :
        return wf.label()

    elif command == 'showcounts' :
        return wf.showcounts()

    elif command == 'showlabels' :
        return wf.showlabels()

    elif command == 'phylogeny' :
        return wf.phylogeny()
    
    elif command == 'heatmap' :
        return wf.heatmap()

    elif command == 'wasabi' :
        return wf.wasabi()

    else :
        print >> stderr, "'%s' appears to be partially implemented!" % command

    return 1

if __name__ == '__main__' :
    try :
        exit(main())
    except KeyboardInterrupt :
        print >> stderr, "Killed by user"

