import os
import string
import mmap
import numpy as np

from records import NpyFastaRecord

_complement = string.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')
complement  = lambda s: s.translate(_complement)

class Fasta(dict):
    def __init__(self, fasta_name, record_class=NpyFastaRecord):
        """
            >>> from pyfasta import Fasta, FastaRecord

            >>> f = Fasta('tests/data/three_chrs.fasta', 
            ...                          record_class=FastaRecord)
            >>> sorted(f.keys())
            ['chr1', 'chr2', 'chr3']

        slicing returns an object.
            >>> f['chr1']
            FastaRecord('tests/data/three_chrs.fasta.flat', 0..80)

        extract sequence with normal python syntax
            >>> f['chr1'][:10]
            'ACTGACTGAC'

        take the first basepair in each codon...
            >>> f['chr1'][0:10:3]
            'AGTC'

        """
        self.fasta_name = fasta_name
        self.record_class = record_class
        self.index, self.prepared = self.record_class.prepare(self, self.gen_seq_info())

        self.chr = {}

    @classmethod
    def as_kmers(klass, seq, k, overlap=0):
        kmax = len(seq)
        assert overlap < k, ('overlap must be < kmer length')
        i = 0
        while i < kmax:
            yield i, seq[i:i + k]
            i += k - overlap

    def gen_seq_info(self):
        """remove all newlines from the sequence in a fasta file
        and generate starts, stops to be used by the record class"""
        fh = open(self.fasta_name, 'r+')
        mm = mmap.mmap(fh.fileno(), os.path.getsize(self.fasta_name))

        # do the flattening (remove newlines)
        sheader = mm.find('>')
        snewline = mm.find('\n', sheader)
        idx = {}
        start = 0
        len_mm = len(mm)
        while sheader < len_mm:
            header = mm[sheader:snewline + 1]
            sheader = mm.find('>', snewline)
            if sheader == -1: sheader = len(mm)

            seq  = mm[snewline + 1: sheader].replace('\n','')

            stop = start + len(seq)
            yield header[1:].strip(), start, stop, seq

            start = stop
            snewline = mm.find('\n', sheader)
        fh.close()

    def __len__(self):
        # might not work for all backends?
        return len(self.index)

    def iterkeys(self):
        for k in self.index.iterkeys(): yield k

    def keys(self):
        return self.index.keys()
    
    def __contains__(self, key):
        return key in self.index

    def __getitem__(self, i):
        # this implements the lazy loading
        if i in self.chr:
            return self.chr[i]

        c = self.index[i]
        self.chr[i] = self.record_class(self.prepared, c[0], c[1])
        return self.chr[i]

    def sequence(self, f, asstring=True, auto_rc=True
            , exon_keys=None):
        """
        take a feature and use the start/stop or exon_keys to return
        the sequence from the assocatied fasta file:
        f: a feature
        asstring: if true, return the sequence as a string
                : if false, return as a numpy array
        auto_rc : if True and the strand of the feature == -1, return
                  the reverse complement of the sequence

            >>> from pyfasta import Fasta
            >>> f = Fasta('tests/data/three_chrs.fasta')
            >>> f.sequence({'start':1, 'stop':2, 'strand':1, 'chr': 'chr1'})
            'AC'

            >>> f.sequence({'start':1, 'stop':2, 'strand': -1, 'chr': 'chr1'})
            'GT'

            >>> f.index
            {'chr3': (160, 3760), 'chr2': (80, 160), 'chr1': (0, 80)}

        NOTE: these 2 are reverse-complement-ary because of strand
        #>>> f.sequence({'start':10, 'stop':12, 'strand': -1, 'chr': 'chr1'})
            'CAG'
            >>> f.sequence({'start':10, 'stop':12, 'strand': 1, 'chr': 'chr1'})
            'CTG'


            >>> f.sequence({'start':10, 'stop':12, 'strand': -1, 'chr': 'chr3'})
            'TGC'
            >>> f.sequence({'start':10, 'stop':12, 'strand': 1, 'chr': 'chr3'})
            'GCA'

            >>> f['chr3'][:][-10:]
            'CGCACGCTAC'

        
        a feature can have exons:
            >>> feat = dict(start=9, stop=19, strand=1, chr='chr1'
            ...    , exons=[(9,11), (13, 15), (17, 19)])

        by default, it just returns the full sequence between start
        and stop.
            >>> f.sequence(feat)
            'ACTGACTGACT'

        but if exon_keys is set to an iterable, it will search for
        those keys and will use the first to create a sequence and
        return the concatenated result.
            >>> f.sequence(feat, exon_keys=('rnas', 'exons'))
            'ACTACTACT'

        Note that sequence is 2 characters shorter than the entire
        feature, to account for the introns at base-pairs 12 and 16.

        Also note, it looks for an item with key of 'rnas', and didn't
        fine one, so it continued on to 'exons'. If it doesn't find
        any of the exon keys, it will fall back on the start, stop of
        the feature:
            >>> f.sequence(feat, exon_keys=('fake', 'also_fake'))
            'ACTGACTGACT'
        """
        assert 'chr' in f and f['chr'] in self, (f, f['chr'], self.keys())
        fasta    = self[f['chr']]
        sequence = None
        if not exon_keys is None:
            sequence = self._seq_from_keys(f, fasta, exon_keys)

        if sequence is None:
            sequence = fasta[(f['start'] - 1): f['stop']]

        if auto_rc and f.get('strand') in (-1, '-1', '-'):
            sequence = complement(sequence)[::-1]

        if asstring: return sequence
        return np.array(sequence, dtype='c')

    def _seq_from_keys(self, f, fasta, exon_keys, base='locations'):
        """Internal:
        f: a feature dict
        fasta: a Fasta object
        exon_keys: an iterable of keys, to look for start/stop
                   arrays to get sequence.
        base: if base ('locations') exists, look there fore the
        exon_keys, not in the base of the object:
            {'name': 'OS11G42690', 'stop': 25210251, 'locations':
            {'CDS': [(25210018, 25210251)]}, 'start': 25210018, 'chr':
            '11', 'strand': -1} set(['TRNA', 'CDS'])
        """
        fbase = f.get(base, f)
        for ek in exon_keys:
            if not ek in fbase: continue
            locs = fbase[ek]
            seq = ""
            for start, stop in locs:
                seq += fasta[start -1:stop]
            return seq
        return None

    def iteritems(self):
        for k in self.keys():
            yield k, self[k]
