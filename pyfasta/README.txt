==================================================
pyfasta: pythonic access to fasta sequence files.
==================================================


:Author: Brent Pedersen (brentp)
:Email: bpederse@gmail.com
:License: MIT


Implementation
==============

Requires Python >= 2.5. Stores a flattened version of the fasta file without 
spaces or headers and uses either a mmap of numpy binary format or fseek/fread so the
*sequence data is never read into memory*. Saves a pickle (.gdx) of the start, stop 
(for fseek/mmap) locations of each header in the fasta file for internal use.
Now supports the numpy array interface.
When the underlying sequence file contains fewer than 150 headers (e.g. fewer than 150 
chromosomes), the numpy binary format will be used and access will be significantly faster.
For greater than 150 sequences, fseek/fread are used.

Usage
=====

::

    >>> from pyfasta import Fasta

    >>> f = Fasta('tests/data/three_chrs.fasta')
    >>> sorted(f.keys())
    ['chr1', 'chr2', 'chr3']

    >>> f['chr1']
    NpyFastaRecord('tests/data/three_chrs.fasta.flat.npy', 0..80)

Slicing
-------
::

    >>> f['chr1'][:10]
    'ACTGACTGAC'

    # get the 1st basepair in every codon (it's python yo)
    >>> f['chr1'][::3]
    'AGTCAGTCAGTCAGTCAGTCAGTCAGT'


    # the index stores the start and stop of each header from the flattened 
    # fasta file. (you should never need this)
    >>> f.index
    {'chr3': (160, 3760), 'chr2': (80, 160), 'chr1': (0, 80)}


    # can query by a 'feature' dictionary
    >>> f.sequence({'chr': 'chr1', 'start': 2, 'stop': 9})
    'CTGACTGA'

    # same as:
    >>> f['chr1'][1:9]
    'CTGACTGA'

    # with reverse complement for - strand
    >>> f.sequence({'chr': 'chr1', 'start': 2, 'stop': 9, 'strand': '-'})
    'TCAGTCAG'

    # for files with < 150 sequences, it's possible to get back a numpy array directly
    >>> f['chr1'].tostring = False
    >>> f['chr1'][:10] # doctest: +NORMALIZE_WHITESPACE
    memmap(['A', 'C', 'T', 'G', 'A', 'C', 'T', 'G', 'A', 'C'], dtype='|S1')


---------------------
Numpy Array Interface
---------------------
::

    # FastaRecords support the numpy array interface.
    >>> import numpy as np
    >>> a = np.array(f['chr2'])
    >>> a.shape[0] == len(f['chr2'])
    True

    >>> a[10:14]
    array(['A', 'A', 'A', 'A'], 
          dtype='|S1')

    # mask a sub-sequence:
    >>> a[11:13] = np.array('N', dtype='c')
    >>> a[10:14].tostring()
    'ANNA'

   

    # cleanup (though for real use these will remain for faster access)
    >>> import os
    >>> os.unlink('tests/data/three_chrs.fasta.gdx')
    >>> os.unlink('tests/data/three_chrs.fasta.flat.npy')
