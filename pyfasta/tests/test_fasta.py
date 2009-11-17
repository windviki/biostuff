from pyfasta import Fasta, NpyFastaRecord, MemoryRecord, FastaRecord
import unittest
import os

class FastaTest(unittest.TestCase):

    def setUp(self):
        self.f = Fasta('tests/data/three_chrs.fasta')

    def test_keys(self):

        self.assertEqual(sorted(self.f.keys()), ['chr1', 'chr2', 'chr3'])

    
    def test_mmap(self):
        seq = self.f['chr2']
        self.assertEqual(seq.__class__, NpyFastaRecord)

        self.assertEqual(seq[0], 'T')
        self.assertEqual(seq[-1], 'T')



        for i in (1, len(seq) -2):
            self.assertEqual(seq[i], 'A')

        for i in (1, len(seq) -3):
            self.assertEqual(seq[i: i + 2], 'AA')
        
        self.assertEqual(seq[1], 'A')

        self.assertEqual(seq[-1], 'T')
        self.assertEqual(seq[-2], 'A')

        self.assertEqual(seq[0], 'T')
        self.assertEqual(seq[6:9], 'AAA')

    def test_shape(self):
        self.assertEqual(len(self.f['chr2']), 80)
        self.assertEqual(len(self.f['chr3']), 3600)
        self.assertEqual(self.f['chr2'][:2], 'TA')

        self.assertEqual(self.f['chr3'][:5], 'ACGCA')

    def test_bounds(self):
        c2 = self.f['chr2']
        self.assertEquals(len(str(c2[0:900])), 80)

        self.assertEquals(c2[99:99], "")
        self.assertEquals(c2[99:99], "")
        self.assertEquals(c2[80:81], "")
        self.assertEquals(c2[79:81], "T")


    def test_tostring(self):
        s = 'TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT'
        self.assertEqual(str(self.f['chr2']), s)

    def test_kmers(self):
        seq = str(self.f['chr2'])

        kmers = list(Fasta.as_kmers(self.f['chr2'], 10))
        self.assertEqual(len(kmers), len(seq) / 10)
        self.assertEqual(kmers[0], (0, seq[:10]))

        seqs = [k[1] for k in kmers]
        self.assertEqual("".join(seqs), seq)
        last_pair = kmers[-1]
        self.assert_(seqs[-1][-1], 'T')

        seq = str(self.f['chr3'])
        kmers = list(Fasta.as_kmers(self.f['chr3'], 1))
        self.assertEquals(kmers[2][0], 2)
        seqs = [k[1] for k in kmers]
        self.assertEqual("".join(seqs), seq)

    def test_kmer_overlap(self):
        chr2 = self.f['chr2']

        kmers = Fasta.as_kmers(chr2, 10, overlap=2)
        for i, k in enumerate(list(kmers)[:-1]):
            self.assertEquals(len(k[1]), 10)
            self.assertEquals(k[0], (i * (10 - 2)))

        kmers = Fasta.as_kmers(chr2, 10, overlap=4)
        seqs = [k[1] for k in kmers]
        paired_seqs = zip(seqs[0:-1], seqs[1:])
        for a, b in paired_seqs:
            if len(a) < 4 or len(b) < 4: continue
            self.assertEqual(a[-4:], b[:4])

    def test_slice_size(self):
        self.assertEqual(self.f['chr3'][:7], 'ACGCATT')
        # take the first basepair of each codon
        self.assertEqual(self.f['chr3'][0:7:3], 'ACT')
        # take the 2nd basepair of each codon
        self.assertEqual(self.f['chr3'][1:7:3], 'CA')

    def test_slice(self):
        f = self.f
        self.assertEqual(str(f['chr3'][:4]), 'ACGC')

    
        seq = 'TACGCACGCTAC'
        self.assertEqual(seq, f['chr3'][-12:])

        self.assertEqual(seq[:4], f['chr3'][-12:-8])

        self.assertEqual(f['chr3'][0:5][::-1], 'ACGCA')
        self.assertEqual(f['chr3'][0:5], 'ACGCA')

    def tearDown(self):
        os.unlink('tests/data/three_chrs.fasta.flat')
        os.unlink('tests/data/three_chrs.fasta.gdx')

class FSeekClassTest(unittest.TestCase):
    def setUp(self):
        self.f = Fasta('tests/data/three_chrs.fasta', record_class=FastaRecord)

    def test_bounds_fseek(self):
        c2 = self.f['chr2']
        self.assertEquals(len(str(c2[0:900])), 80)

        self.assertEquals(c2[99:99], "")
        self.assertEquals(c2[99:99], "")
        self.assertEquals(c2[80:81], "")
        self.assertEquals(c2[79:81], "T")
    def test_get(self):
        f = self.f
        self.assertEqual(f['chr3'][0:5][::-1], 'ACGCA')
        seq = 'TACGCACGCTAC'
        self.assertEqual(seq, f['chr3'][-12:])
        self.assertEqual(seq[:4], f['chr3'][-12:-8])

        seq = self.f['chr2']
        self.assertEqual(seq.__class__, FastaRecord)

        self.assertEqual(seq[0], 'T')
        self.assertEqual(seq[-1], 'T')

        for i in (1, len(seq) -2):
            self.assertEqual(seq[i], 'A')

        for i in (1, len(seq) -3):
            self.assertEqual(seq[i: i + 2], 'AA')
        
        self.assertEqual(seq[1], 'A')
    def tearDown(self):
        os.unlink('tests/data/three_chrs.fasta.flat')
        os.unlink('tests/data/three_chrs.fasta.gdx')


class RecordClassTest(unittest.TestCase):
    def setUp(self):
        self.f = Fasta('tests/data/three_chrs.fasta', record_class=MemoryRecord)

    def test_get(self):
        f = self.f
        self.assertEqual(f['chr3'][0:5][::-1], 'ACGCA')
        seq = 'TACGCACGCTAC'
        self.assertEqual(seq, f['chr3'][-12:])
        self.assertEqual(seq[:4], f['chr3'][-12:-8])

        seq = self.f['chr2']
        self.assertEqual(seq.__class__, MemoryRecord)

        self.assertEqual(seq[0], 'T')
        self.assertEqual(seq[-1], 'T')

        for i in (1, len(seq) -2):
            self.assertEqual(seq[i], 'A')

        for i in (1, len(seq) -3):
            self.assertEqual(seq[i: i + 2], 'AA')
        
        self.assertEqual(seq[1], 'A')

    def test_bounds_fseek(self):
        c2 = self.f['chr2']
        self.assertEquals(len(str(c2[0:900])), 80)

        self.assertEquals(c2[99:99], "")
        self.assertEquals(c2[99:99], "")
        self.assertEquals(c2[80:81], "")
        self.assertEquals(c2[79:81], "T")

import numpy as np
class ArrayInterfaceTest(unittest.TestCase):

    def setUp(self):
        self.f = Fasta('tests/data/three_chrs.fasta')

    def test_len(self):
        a = np.array(self.f['chr3'])
        assert a.shape[0] == len(self.f['chr3'])

    def test_copy(self):
        # the array is definitely a copy...
        a = np.array(self.f['chr3'])
        old = self.f['chr3'][1:5]
        self.assertEqual(a.shape[0], 3600)
        a[1:5] = np.array('N', dtype='c')
        c = self.f['chr3'][1:5]
        assert c == old

        assert a[1:5].tostring() == 'NNNN', a[1:5].tostring()


if __name__ == "__main__":
    unittest.main()
