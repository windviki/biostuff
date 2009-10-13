from fasta import *
import sys

def main():
    help = """
    available actions:
        `extract`: extract sequences from a fasta file

    to view the help for a particular action, use:
        pyfasta [action] --help
    e.g.:
        pyfasta extract --help
    """        
    if len(sys.argv) == 1:
        print help
        sys.exit()

    action = sys.argv[1]

    sglobals = globals()
    if not action in sglobals:
        print "%s not a valid action" % action
        print help
        sys.exit()
    
    globals()[action](sys.argv[2:])

def extract(args):
    """
    >>> extract(['--fasta', 'tests/data/three_chrs.fasta', 'chr2'])
    TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT
    """

    import optparse
    parser = optparse.OptionParser("""... send in options + a list of sequences to extract
               e.g.:
               pyfasta extract --fasta some.fasta --header at2g26540 at3g45640""") 
    parser.add_option("--fasta", dest="fasta", help="path to the fasta file")
    parser.add_option("--header", dest="header", help="include headers", action="store_true", default=False)
    parser.add_option("--file", dest="file", help="if this flag is used, the sequences to extract" \
                                                  + "are read from the file specified in args"
                      , action="store_true", default=False)
    options, seqs = parser.parse_args(args)
    if not (options.fasta and len(seqs)):
        sys.exit(parser.print_help())

    f = Fasta(options.fasta)
    if options.file:
        seqs = (x.strip() for x in open(seqs[0]))
     
    for seqname in seqs:
        seq = f[seqname]
        if options.header:
            print ">%s" % seqname
        print seq


if __name__ == "__main__":
    main()
