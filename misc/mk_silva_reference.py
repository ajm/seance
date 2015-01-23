from sys import stderr, argv, exit
from os.path import splitext


def fix_seq(s) :
    return s.replace(' ', '').upper().replace('U', 'T')

def main() :
    if len(argv) != 2 :
        print >> stderr, "Usage: %s <SILVA.tgz>\n" % argv[0]
        exit(1)

    seq = ""
    for line in open(argv[1]) :
        line = line.rstrip()
        
        if not line :
            continue
        
        if line.startswith('>') :
            if seq :
                print fix_seq(seq)
            seq = ""
            print ">" + line.split()[1]
        else :
            seq += line

    print fix_seq(seq)



if __name__ == '__main__' :
    try :
        main()

    except KeyboardInterrupt :
        print >> stderr, "Killed by user..."
        exit(1)

