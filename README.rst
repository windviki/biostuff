===============================================================================
biostuff 
===============================================================================

miscellaneous modules for bioinformatics with tests and documentation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Skidmarks_
----------
find runs (non-randomness) in sequences 

::

    >>> from skidmarks import gap_test, wald_wolfowitz, auto_correlation, serial_test
    >>> serial_test('110000000000000111111111111')
    {'chi': 18.615384615384617, 'p': 0.00032831021826061683}



nwalign_ 
--------
fast Needleman-Wunsch_ global alignment in cython. command-line and python usage
::

    >>> TODO

pyfasta_
--------
pythonic access to fasta sequence files
::

    >>> from pyfasta import Fasta
    >>> f = Fasta('some.fasta')
    >>> f.keys()
    ['chr1', 'chr2', 'chr3']

    >>> f['chr1'][10:20]
    'actgatcgga'



simpletable_
------------
pytables_ wrapper for easy access to s structured data.
::

    >> TODO



.. _Skidmarks: http://pypi.python.org/pypi/skidmarks/
.. _SimpleTable: http://pypi.python.org/pypi/simpletable/
.. _nwalign: http://pypi.python.org/pypi/nwalign/
.. _pyfasta: http://pypi.python.org/pypi/pyfasta/
.. _pytables: http://pytables.org/
.. _Needleman-Wunsch: http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm 
