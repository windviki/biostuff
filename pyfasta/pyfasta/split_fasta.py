from pyfasta import Fasta
import operator
import collections
import string
import sys
import optparse

def newnames(oldname, n):

    p = oldname.rfind("fa")
    if p != -1:
        pattern = oldname[:p] + "%s." + oldname[p:]
    else:
        pattern = oldname + ".%s"
    if n == 1:
        names = [pattern % "split"]
    else:
        names = [pattern % string.letters[i] for i in range(n)]
    print "creating new files:"
    print "\n".join(names)
    return names


def print_to_fh(fh, fasta, lens, seqinfo):
    key, seqlen = seqinfo
    lens[fh.name] += seqlen
    f = fasta
    assert len(str(f[key])) == seqlen, (key, seqlen, len(str(f[key])))
    print >>fh, ">%s" % key
    print >>fh, str(f[key])


def format_kmer(seqid, start, kmer_len):
    return "%s_%i (%i-mers)" % (seqid, start + 1, kmer_len)

def split(args):
    parser = optparse.OptionParser("""\
   split a fasta file into separated files.
        pyfasta split -n 6 [-k 5000 ] some.fasta
    the output will be some.1.fasta, some.2.fasta ... some.6.fasta
    the sizes will be as even as reasonable.
   """)

    parser.add_option("-n", "--n", type="int", dest="nsplits", 
                            help="number of new files to create")
    parser.add_option("-k", "--kmers", type="int", dest="kmers", default=-1,
                     help="""\
    split big files into pieces of this size default
    default of -1 means do not split a reasonable value would be 10K""")
    options, (fasta, ) = parser.parse_args(args)
    if not (fasta):
        sys.exit(parser.print_help())

    names = newnames(fasta, options.nsplits)

    fhs = [open(n, 'wb') for n in names]
    f = Fasta(fasta)
    if options.kmers == -1:
        return without_kmers(f, fhs)
    else: 
        return with_kmers(f, fhs, options.kmers)

def with_kmers(f, fhs, k):
    i = 0
    for seqid in f.keys():
        seq = f[seqid]
        for (start0, subseq) in seq.as_kmers(k):
            fh = fhs[i % len(fhs)]
            print >>fh, ">%s" % format_kmer(seqid, start0, k)
            print >>fh, subseq
            i += 1

def without_kmers(f, fhs):
    name2fh = dict([(fh.name, fh) for fh in fhs])
    # sort small to large.
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
        if added:
            continue

        print_to_fh(fh, f, lens, items[l1])
        l1 -= 1
        n_added += 1

    if l0 == l1:
        fh = fhs[l0 % len(fhs)]
        print_to_fh(fh, f, lens, items[l0])


def find_name_from_len(lmin, lens):
    for fname, l in lens.iteritems():
        if l == lmin: 
            return fname
    raise Exception('name not found')
