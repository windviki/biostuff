from cnwalign import global_align, global_align_no_matrix
import cnwalign
__doc__ = cnwalign.__doc__


def main():
    import sys
    import optparse
    parser = optparse.OptionParser(usage="""
    %prog [options] seq1 seq2 
    """)
    parser.add_option("--gap", dest="gap", help="gap extend penalty (must be integer <= 0)", type="int", default=-1)
    parser.add_option("--gap_init", dest="gap_init", help="gap start penalty (must be integer <= 0)", type="int", default=-1)
    parser.add_option("--match", dest="match", help="match score (must be integer > 0)", type="int", default=1)
    parser.add_option("--mismatch", dest="mismatch", help="gap penalty (must be integer < 0)", type="int", default=-1)
    parser.add_option("--matrix", dest="matrix", help="scoring matrix in ncbi/data/ format,\
                                      if not specificied, match/mismatch are used", default=None)

    try:
        options, args = parser.parse_args()
    except:
        sys.exit(parser.print_help())
    if len(args) != 2:
        sys.exit(parser.print_help())
    print "\n".join(global_align(args[0], args[1], options.gap, options.match,
                                 options.mismatch, int(options.gap_init), options.matrix))

if __name__ == "__main__":
    main()
