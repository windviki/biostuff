==================================================
pyfasta: pythonic access to fasta sequence files.
==================================================


:Author: Brent Pedersen (brentp)
:Email: bpederse@gmail.com
:License: MIT

.. contents ::

Implementation
==============

Requires Python >= 2.5. Stores a flattened version of the fasta file without 
spaces or headers and uses either a mmap of numpy binary format or fseek/fread so the
*sequence data is never read into memory*. Saves a pickle (.gdx) of the start, stop 
(for fseek/mmap) locations of each header in the fasta file for internal use.

Usage
=====
::

    >>> from pyfasta import Fasta

    >>> f = Fasta('tests/data/three_chrs.fasta')
    >>> sorted(f.keys())
    ['chr1', 'chr2', 'chr3']

    >>> f['chr1']
    NpyFastaRecord(0..80)


Slicing
-------
::

    >>> f['chr1'][:10]
    'ACTGACTGAC'

    # get the 1st basepair in every codon (it's python yo)
    >>> f['chr1'][::3]
    'AGTCAGTCAGTCAGTCAGTCAGTCAGT'

    # can query by a 'feature' dictionary
    >>> f.sequence({'chr': 'chr1', 'start': 2, 'stop': 9})
    'CTGACTGA'

    # same as:
    >>> f['chr1'][1:9]
    'CTGACTGA'

    # with reverse complement (automatic for - strand)
    >>> f.sequence({'chr': 'chr1', 'start': 2, 'stop': 9, 'strand': '-'})
    'TCAGTCAG'

Numpy
=====

The default is to use a memmaped numpy array as the backend. In which case it's possible to
get back an array directly...
::

    >>> f['chr1'].tostring = False
    >>> f['chr1'][:10] # doctest: +NORMALIZE_WHITESPACE
    memmap(['A', 'C', 'T', 'G', 'A', 'C', 'T', 'G', 'A', 'C'], dtype='|S1')

    >>> import numpy as np
    >>> a = np.array(f['chr2'])
    >>> a.shape[0] == len(f['chr2'])
    True

    >>> a[10:14] # doctest: +NORMALIZE_WHITESPACE
    array(['A', 'A', 'A', 'A'], dtype='|S1')

mask a sub-sequence
::

    >>> a[11:13] = np.array('N', dtype='c')
    >>> a[10:14].tostring()
    'ANNA'


Backends (Record class)
=======================
It's also possible to specify another record class as the underlying work-horse
for slicing and reading. Currently, there's just the default: 

  * NpyFastaRecord which uses numpy memmap
  * FastaRecord, which uses using fseek/fread
  * MemoryRecord which reads everything into memory and must reparse the original
    fasta every time.

it's possible to create your own using a sub-class of FastaRecord. see the source 
for details.
Next addition will be a pytables/hdf5 backend.
::

    >>> from pyfasta import FastaRecord # default is NpyFastaRecord
    >>> f = Fasta('tests/data/three_chrs.fasta', record_class=FastaRecord)
    >>> f['chr1']
    FastaRecord('tests/data/three_chrs.fasta.flat', 0..80)

other than the repr, it should behave exactly like the Npy record class backend


Command Line Interface
======================
there's also a command line interface to manipulate / view fasta files.
the `pyfasta` executable is installed via setuptools, running it will show
help text.

split a fasta file into 6 new files of relatively even size:

  $ pyfasta **split** -n 6 original.fasta

create 1 new fasta file with the sequence split into 10K-mers:

  $ pyfasta **split** -n 1 -k 10000 original.fasta

2 new fasta files with the sequence split into 10K-mers with 2K overlap:

  $ pyfasta **split** -n 2 -k 10000 -o 2000 original.fasta


show some info about the file (and show gc content):

  $ pyfasta **info** --gc test/data/three_chrs.fasta


extract sequence from the file. use the header flag to make
a new fasta file. the args are a list of sequences to extract.

  $ pyfasta **extract** --header --fasta test/data/three_chrs.fasta seqa seqb seqc


cleanup 
=======
(though for real use these will remain for faster access)
::

    >>> import os
    >>> os.unlink('tests/data/three_chrs.fasta.gdx')
    >>> os.unlink('tests/data/three_chrs.fasta.flat')
