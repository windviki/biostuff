from fasta import *
import sys

def main():
    help = """
    available actions:
        `extract`: extract sequences from a fasta file
        `info`: show info about the fasta file and exit.

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

def info(args):
    """
    >>> info(['--fasta', 'tests/data/three_chrs.fasta'])
    """
    import optparse
    parser = optparse.OptionParser("... a fasta file and print out the results in order of length."
                                   "by default only the first 30 sequences are printed.")

    parser.add_option("--all", dest="all", help="include headers",
                      action="store_true", default=False)
    parser.add_option("--fasta", dest="fasta", help="path to the fasta file")
    options, _ = parser.parse_args(args)
    if not (options.fasta):
        sys.exit(parser.print_help())
    import operator
    
    f = Fasta(options.fasta)
    info = sorted([(k, len(seq)) for k, seq in f.iteritems()], 
                  key=operator.itemgetter(1), reverse=True)
    
    total_len = sum(l for k, l in info)
    nseqs = len(info)
    if not options.all: info = info[:30]
    for k, l in info:
        print ">%s length: %i" % (k, l)

    if total_len > 1000000:
        total_len = "%.3fM" % (total_len / 1000000.)
    print
    print "%s basepairs in %i sequences" % (total_len, nseqs)


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
