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

    def __init__(self, fh, start, stop):

        self.fh      = fh
        self.stop    = stop
        self.start   = start

    def __len__(self):
        return self.stop - self.start

    @classmethod
    def prepare(self, filename):
        return open(filename, "r")

    @classmethod
    def modify_flat(self, flat_file):
        return None
    
    def _adjust_slice(self, islice):
        if not islice.start is None and islice.start < 0:
            istart = self.stop + islice.start
        else:
            istart = self.start + (0 if islice.start is None else islice.start)


        if not islice.stop is None and islice.stop < 0:
            istop = self.stop + islice.stop
        else:
            istop = self.stop if islice.stop is None else (self.start + islice.stop)
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

        fh.seek(istart)

        l = istop - istart

        if islice.step in (1, None):
            return fh.read(l)

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
    ext = ".npy"
    def __init__(self, mm, start, stop, tostring=True):
        self.mm = mm
        self.start = start
        self.stop = stop
        self.tostring = tostring

    def __repr__(self):
        return "%s(%i..%i)" % (self.__class__.__name__,
                                   self.start, self.stop)

    @classmethod
    def prepare(self, filename):
        return np.load(filename, mmap_mode="r")

    @classmethod
    def modify_flat(self, flat_file):
        mm = np.memmap(flat_file, dtype="S1", mode="r+")
        a = np.array(mm, dtype="S1")
        del mm
        np.save(flat_file, a)

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
        self.gdx = fasta_name + ".gdx"
        self.flat = fasta_name + record_class.ext
        self.index = self.flatten()

        self.chr = {}

    @classmethod
    def is_up_to_date(klass, a, b):
        return os.path.exists(a) and os.stat(a).st_mtime > os.stat(b).st_mtime


    def flatten(self):
        """remove all newlines from the sequence in a fasta file"""

        if Fasta.is_up_to_date(self.gdx, self.fasta_name) and \
            Fasta.is_up_to_date(self.flat, self.fasta_name):
                self.prepared = self.record_class.prepare(self.flat)
                return Fasta._load_index(self.gdx)

        fh = open(self.fasta_name, 'r+')
        mm = mmap.mmap(fh.fileno(), os.path.getsize(self.fasta_name))

        out = open(self.flat, 'wb')

        # do the flattening (remove newlines)
        sheader = mm.find('>')
        snewline = mm.find('\n', sheader)
        idx = {}
        pos = 0
        len_mm = len(mm)
        while sheader < len_mm:
            header = mm[sheader:snewline + 1]
            #out.write(header)

            sheader = mm.find('>', snewline)
            if sheader == -1: sheader = len(mm)

            seq  = mm[snewline + 1: sheader].replace('\n','')
            out.write(seq)
            out.flush()
            p1 = out.tell()
            idx[header[1:].strip()] = (pos, p1)
            pos = p1 

            snewline = mm.find('\n', sheader)

        p = open(self.gdx, 'wb')
        cPickle.dump(idx, p)
        p.close(); fh.close(); out.close()

        self.record_class.modify_flat(self.flat)

        self.prepared = self.record_class.prepare(self.flat)
        #self.prepared = np.load(self.flat, mmap_mode="r")

        return idx




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
        """ old
        else:
            self.chr[i] = FastaRecord(self.fh, c[0], c[1])
        """
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

