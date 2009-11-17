import os
import cPickle
import string
import mmap
import numpy as np
import sys

_complement = string.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')
complement  = lambda s: s.translate(_complement)


class FastaRecord(object):
    __slots__ = ('fh', 'start', 'stop')
    ext = ".flat"
    idx = ".gdx"

    @classmethod
    def is_current(klass, fasta_name):
        utd = Fasta.is_up_to_date(fasta_name + klass.ext, fasta_name) 
        if not utd: return False
        return Fasta.is_up_to_date(fasta_name + klass.idx, fasta_name)

    def __init__(self, fh, start, stop):

        self.fh      = fh
        self.stop    = stop
        self.start   = start

    def __len__(self):
        return self.stop - self.start

    @classmethod
    def prepare(klass, fasta_obj, seqinfo_generator):
        """
        returns the __getitem__'able index. and the thing to get the seqs from.
        """
        f = fasta_obj.fasta_name
        if klass.is_current(f):
            return cPickle.load(open(f + klass.idx)), \
                               klass.modify_flat(f + klass.ext)

        idx = {}
        flatfh = open(f + klass.ext, 'wb')
        for seqid, start, stop, seq in seqinfo_generator:
            flatfh.write(seq)
            idx[seqid] = (start, stop)
            
        cPickle.dump(idx, open(f + klass.idx, 'wb'), -1)
        flatfh.close()

        return idx, klass.modify_flat(f + klass.ext)
    
    @classmethod
    def modify_flat(klass, flat_file):
        return open(flat_file, 'rb')
    
    def _adjust_slice(self, islice):
        if not islice.start is None and islice.start < 0:
            istart = self.stop + islice.start
        else:
            istart = self.start + (0 if islice.start is None else islice.start)


        if not islice.stop is None and islice.stop < 0:
            istop = self.stop + islice.stop
        else:
            istop = self.stop if islice.stop is None else (self.start + islice.stop)

        # this will give empty string
        if istart > self.stop: return self.stop, self.stop 
        if istop > self.stop: 
            istop = self.stop
        return istart, istop

    def __getitem__(self, islice):
        fh = self.fh
        fh.seek(self.start)
        if isinstance(islice, (int, long)):
            if islice < 0:
                fh.seek(self.stop + islice)
            else:
                fh.seek(self.start + islice)
            return fh.read(1)

        if islice.start == 0 and islice.stop == sys.maxint:
            if islice.step in (1, None):
                return fh.read(self.stop - self.start)
            return fh.read(self.stop - self.start)[::islice.step]
        
        istart, istop = self._adjust_slice(islice)
        if istart is None: return ""
        l = istop - istart
        if l == 0: return ""


        fh.seek(istart)

        if islice.step in (1, None):
            try:
                return fh.read(l)
            except:
                print "\n", l, fh
                raise

        return fh.read(l)[::islice.step]


    def __str__(self):
        return self[:]

    def __repr__(self):
        return "%s('%s', %i..%i)" % (self.__class__.__name__, self.fh.name,
                                   self.start, self.stop)

    @property
    def __array_interface__(self):
        return {
            'shape': (len(self), ),
            'typestr': '|S1',
            'version': 3,
            'data': buffer(self)
        }

class NpyFastaRecord(FastaRecord):
    __slots__ = ('start', 'stop', 'mm', 'tostring')

    def __init__(self, mm, start, stop, tostring=True):
        self.mm = mm
        self.start = start
        self.stop = stop
        self.tostring = tostring

    def __repr__(self):
        return "%s(%i..%i)" % (self.__class__.__name__,
                                   self.start, self.stop)

    @classmethod
    def modify_flat(klass, flat_file):
        mm = np.memmap(flat_file, dtype="S1", mode="r+")
        return mm

    def getdata(self, islice):
        if isinstance(islice, (int, long)):
            if islice >= 0:
                islice += self.start
            else:
                islice += self.stop
            return self.mm[islice]

        start, stop = self._adjust_slice(islice)
        return self.mm[start:stop:islice.step]

    def __getitem__(self, islice):
        d = self.getdata(islice)
        return d.tostring() if self.tostring else d

    @property
    def __array_interface__(self):
        return {
            'shape': (len(self), ),
            'typestr': '|S1',
            'version': 3,
            'data': self[:]
        }


class MemoryRecord(FastaRecord):
    @classmethod
    def prepare(klass, fasta_obj, seqinfo_generator):
        f = fasta_obj.fasta_name
        seqs = {}
        idx = {}
        for seqid, start, stop, seq in seqinfo_generator:
            seqs[seqid] = (seq, None)
            
        return seqs, seqs

    def __init__(self, _, seq, _none):
        self.seq = seq

    def __getitem__(self, slice):
        return self.seq.__getitem__(slice)

    def __len__(self):
        return len(self.seq)

class Fasta(dict):
    def __init__(self, fasta_name, record_class=NpyFastaRecord):
        """
            >>> from pyfasta import Fasta

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
        #self.gdx = fasta_name + record_class.idx
        #self.flat = fasta_name + record_class.ext
        self.index, self.prepared = self.record_class.prepare(self, self.gen_seq_info())

        self.chr = {}

    @classmethod
    def is_up_to_date(klass, a, b):
        return os.path.exists(a) and os.stat(a).st_mtime > os.stat(b).st_mtime

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

            # TODO: make this more efficient.
            stop = start + len(seq)
            yield header[1:].strip(), start, stop, seq

            start = stop
            snewline = mm.find('\n', sheader)


    def __len__(self):
        # might not work for all backends?
        return len(self.index)

    def iterkeys(self):
        for k in self.keys(): yield k
    def keys(self):
        return self.index.keys()
    
    def __contains__(self, key):
        return key in self.index.keys()

    def __getitem__(self, i):
        # this implements the lazy loading
        if i in self.chr:
            return self.chr[i]

        c = self.index[i]
        self.chr[i] = self.record_class(self.prepared, c[0], c[1])
        return self.chr[i]

    @classmethod
    def _load_index(self, path):
        """ """
        gdx = open(path, 'rb')
        try:
            return cPickle.load(gdx)
        finally:
            gdx.close()

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



if __name__ == "__main__":
    import doctest
    doctest.testmod()

