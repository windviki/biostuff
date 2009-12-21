from pyfasta import Fasta
from pyfasta.records import NpyFastaRecord, MemoryRecord, FastaRecord
record_classes = [NpyFastaRecord, MemoryRecord, FastaRecord]
try:
    from pyfasta.records import TCRecord
    record_classes.append(TCRecord)
except ImportError:
    pass

import os
import shutil
from nose.tools import assert_raises
import numpy as np
import glob

def _cleanup():
    for f in glob.glob("tests/data/three_chrs_t.fasta*"):
        os.unlink(f)
    for f in glob.glob("tests/data/three_chrs.fasta.*"):
        os.unlink(f)


def test_classes():

    for inplace in (True, False):
        for klass in record_classes:
            shutil.copyfile('tests/data/three_chrs.fasta', 'tests/data/three_chrs_t.fasta')
            f = Fasta('tests/data/three_chrs_t.fasta', record_class=klass, flatten_inplace=inplace)
            yield check_keys, f
            yield check_misc, f, klass
            yield check_contains, f
            yield check_shape, f
            yield check_bounds, f
            yield check_tostring, f
            yield check_kmers, f
            yield check_kmer_overlap, f
            yield check_slice_size, f
            yield check_slice, f
            yield check_full_slice, f
            yield check_array_copy, f
            yield check_array, f

            del f

            yield check_reload, klass

            _cleanup()


def check_keys(f):
    assert sorted(f.keys()) == ['chr1', 'chr2', 'chr3']
    assert sorted(f.iterkeys()) == ['chr1', 'chr2', 'chr3']


def check_reload(klass):
    f = Fasta('tests/data/three_chrs_t.fasta', record_class=klass)
    assert f

def check_full_slice(f):
    for k in f.keys():
        assert str(f[k]) == f[k][:]
        assert str(f[k]) == f[k][0:]


        assert str(f[k])[::2] == f[k][0:][::2]
        assert str(f[k])[::2] == f[k][:][::2]

def check_contains(f):

    assert not '_________' in f
    assert 'chr2' in f
    
def check_misc(f, klass):
    seq = f['chr2']
    assert seq.__class__ == klass

    assert (seq[0] == 'T')
    assert (seq[-1] == 'T')

    assert (f['chr3'][0:5][::-1] == 'ACGCA')


    for i in (1, len(seq) -2):
        assert (seq[i] == 'A')

    for i in (1, len(seq) -3):
        assert (seq[i: i + 2] == 'AA')
    
    assert (seq[1] == 'A')

    assert (seq[-1] == 'T')
    assert (seq[-2] == 'A')

    assert (seq[0] == 'T')
    assert (seq[6:9] == 'AAA')

    seq = 'TACGCACGCTAC'
    assert (seq == f['chr3'][-12:])
    assert (seq[:4] == f['chr3'][-12:-8])

def check_shape(f):
    assert (len(f['chr2']) == 80)
    assert (len(f['chr3']) == 3600)
    assert (f['chr2'][:2] == 'TA')

    assert (f['chr3'][:5] == 'ACGCA')

def check_bounds(f):
    c2 = f['chr2']
    assert (len(str(c2[0:900])) == 80)

    assert (c2[99:99] == "")
    assert (c2[99:99] == "")
    assert (c2[80:81] == "")
    assert (c2[79:81] == "T")

    assert (c2[800:810] == "")
    assert_raises(IndexError, lambda: c2[-800])


    assert (c2[-800:-810] == "")

def check_array(f):

    s = f.sequence({'chr': 'chr2', 'start': 1, 'stop':3}, asstring=False)
    assert np.all(s == np.array(['T', 'A', 'A'], dtype="|S1"))


def check_tostring(f):
    s = 'TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT'
    assert (str(f['chr2']) == s)

def check_kmers(f):
    seq = str(f['chr2'])

    kmers = list(Fasta.as_kmers(f['chr2'], 10))
    assert (len(kmers) == len(seq) / 10)
    assert (kmers[0] == (0, seq[:10]))

    seqs = [k[1] for k in kmers]
    assert ("".join(seqs) == seq)
    last_pair = kmers[-1]
    assert (seqs[-1][-1] == 'T')

    seq = str(f['chr3'])
    kmers = list(Fasta.as_kmers(f['chr3'], 1))
    assert (kmers[2][0] == 2)
    seqs = [k[1] for k in kmers]
    assert ("".join(seqs) == seq)


def check_kmer_overlap(f):
    chr2 = f['chr2']

    kmers = Fasta.as_kmers(chr2, 10, overlap=2)
    for i, k in enumerate(list(kmers)[:-1]):
        assert (len(k[1]) == 10)
        assert (k[0] == (i * (10 - 2)))

    kmers = Fasta.as_kmers(chr2, 10, overlap=4)
    seqs = [k[1] for k in kmers]
    paired_seqs = zip(seqs[0:-1], seqs[1:])
    for a, b in paired_seqs:
        if len(a) < 4 or len(b) < 4: continue
        assert (a[-4:] == b[:4])

def check_slice_size(f):
    assert (f['chr3'][:7] == 'ACGCATT')
    # take the first basepair of each codon
    assert (f['chr3'][0:7:3] == 'ACT')
    # take the 2nd basepair of each codon
    assert (f['chr3'][1:7:3] == 'CA')

def check_slice(f):
    assert (str(f['chr3'][:4]) == 'ACGC')


    seq = 'TACGCACGCTAC'
    assert (seq == f['chr3'][-12:])

    assert (seq[:4] == f['chr3'][-12:-8])

    assert (f['chr3'][0:5][::-1] == 'ACGCA')
    assert (f['chr3'][0:5] == 'ACGCA')


def check_array_copy(f):
    # the array is definitely a copy...
    a = np.array(f['chr3'])
    old = f['chr3'][1:5]
    assert (a.shape[0] == 3600)
    a[1:5] = np.array('N', dtype='c')
    c = f['chr3'][1:5]
    assert c == old

    assert a[1:5].tostring() == 'NNNN', a[1:5].tostring()


if __name__ == "__main__":
    import nose
    nose.main()
