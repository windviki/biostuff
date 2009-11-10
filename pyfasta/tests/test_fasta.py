from pyfasta import Fasta, NpyFastaRecord
import unittest

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

    def test_tostring(self):
        s = 'TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT'
        self.assertEqual(str(self.f['chr2']), s)


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
        import os
        os.unlink('tests/data/three_chrs.fasta.npy')
        os.unlink('tests/data/three_chrs.fasta.gdx')


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
