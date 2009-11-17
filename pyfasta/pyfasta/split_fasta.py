from pyfasta import Fasta
import operator
import collections
import string
import sys
import optparse


def newnames(oldname, n, kmers=None, overlap=None):
    """
    >>> newnames('some.fasta', 1)
    ['some.split.fasta']

    >>> newnames('some.fasta', 2)
    ['some.a.fasta', 'some.b.fasta']

    >>> newnames('some', 2)
    ['some.a', 'some.b']

    >>> newnames('some.fasta', 2, kmers=1000)
    ['some.a.1Kmer.fasta', 'some.b.1Kmer.fasta']

    >>> newnames('some.fasta', 2, kmers=10000, overlap=2000)
    ['some.a.10Kmer.2Koverlap.fasta', 'some.b.10Kmer.2Koverlap.fasta']

    >>> newnames('some.fasta', 1, kmers=100000, overlap=2000)
    ['some.split.100Kmer.2Koverlap.fasta']

    """
    if kmers and kmers % 1000 == 0: kmers = "%iK" % (kmers/1000)
    if overlap and overlap % 1000 == 0: overlap = "%iK" % (overlap/1000)

    p = oldname.rfind("fa")
    kstr = ("%smer." % kmers) if kmers is not None else ""
    ostr = ("%soverlap." % overlap) if overlap is not None else ""
    if p != -1:
        pattern = oldname[:p] + "%s." + kstr + ostr + oldname[p:]
    else:
        pattern = oldname + kstr + ostr + ".%s"

    
    

    if n == 1:
        names = [pattern % "split"]
    else:
        names = [pattern % string.letters[i] for i in range(n)]
    print >>sys.stderr, "creating new files:"
    print >>sys.stderr, "\n".join(names)
    return names


def print_to_fh(fh, fasta, lens, seqinfo):
    key, seqlen = seqinfo
    lens[fh.name] += seqlen
    f = fasta
    assert len(str(f[key])) == seqlen, (key, seqlen, len(str(f[key])))
    print >>fh, ">%s" % key
    print >>fh, str(f[key])


def format_kmer(seqid, start):
    """
    prints out a header with 1-based indexing.

    >>> format_kmer('chr3', 1000)
    'chr3_1001'
    """
    return "%s_%i" % (seqid, start + 1)

def split(args):
    parser = optparse.OptionParser("""\
   split a fasta file into separated files.
        pyfasta split -n 6 [-k 5000 ] some.fasta
    the output will be some.1.fasta, some.2.fasta ... some.6.fasta
    the sizes will be as even as reasonable.
   """)

    parser.add_option("-n", "--n", type="int", dest="nsplits", 
                            help="number of new files to create")
    parser.add_option("-o", "--overlap", type="int", dest="overlap", 
                            help="overlap in basepairs", default=0)
    parser.add_option("-k", "--kmers", type="int", dest="kmers", default=-1,
                     help="""\
    split big files into pieces of this size in basepairs. default
    default of -1 means do not split a reasonable value would be 10Kbp""")
    options, fasta = parser.parse_args(args)
    if not (fasta and options.nsplits):
        sys.exit(parser.print_help())
    fasta, = fasta

    kmer = options.kmers if options.kmers != -1 else None
    overlap = options.overlap if options.overlap != 0 else None
    names = newnames(fasta, options.nsplits, kmers=kmer, overlap=overlap)

    fhs = [open(n, 'wb') for n in names]
    f = Fasta(fasta)
    if options.kmers == -1:
        return without_kmers(f, fhs)
    else: 
        return with_kmers(f, fhs, options.kmers, options.overlap)

def with_kmers(f, fhs, k, overlap):
    """
    splice the sequences in Fasta object `f` into pieces of length `k` 
    with the given `overlap` the results are written to the array of files
    `fhs`
    """
    i = 0
    for seqid in f.keys():
        seq = f[seqid]
        for (start0, subseq) in Fasta.as_kmers(seq, k, overlap=overlap):
            fh = fhs[i % len(fhs)]
            print >>fh, ">%s" % format_kmer(seqid, start0)
            print >>fh, subseq
            i += 1

def without_kmers(f, fhs):
    """
    long crappy function that does not solve the bin-packing problem.
    but attempts to distribute the sequences in Fasta object `f` evenly
    among the file handles in fhs.
    """
    name2fh = dict([(fh.name, fh) for fh in fhs])
    items = sorted([(key, len(f[key])) for key in f.keys()], 
                   key=operator.itemgetter(1))

    l1 = len(items) - 1
    l0 = 0
    lens = collections.defaultdict(int)

    n_added = 0
    while l0 < l1:
        fh = fhs[n_added % len(fhs)]
        added = False
        if n_added >= len(fhs):

            while l1 > l0:
                lmin = min(lens.itervalues())
                lmax = max(lens.itervalues())
                if float(lmin) / lmax < 0.80:
                    # it's way off, so add a large (l1)
                    name = find_name_from_len(lmin, lens)
                    fh = name2fh[name]
                    print_to_fh(fh, f, lens, items[l1])
                    l1 -= 1
                    added = True
                    n_added += 1
                                

                elif float(lmin) / lmax < 0.94:
                    # it's just a little off so add a small (l0)
                    name = find_name_from_len(lmin, lens)
                    fh = name2fh[name]
                    print_to_fh(fh, f, lens, items[l0])
                    l0 += 1
                    added = True
                    n_added += 1
                else:
                    break

                if not added:
                    break
                # TODO: test this on glycine
                #added = False
        if added:
            continue

        print_to_fh(fh, f, lens, items[l1])
        l1 -= 1
        n_added += 1

    if l0 == l1:
        fh = fhs[l0 % len(fhs)]
        print_to_fh(fh, f, lens, items[l0])


def find_name_from_len(lmin, lens):
    """
    reverse lookup, get name from dict
    """
    for fname, l in lens.iteritems():
        if l == lmin: 
            return fname
    raise Exception('name not found')
