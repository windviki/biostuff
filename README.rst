===============================================================================
biostuff 
===============================================================================

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

    >>> import nwalign as nw
    >>> nw.global_align("CEELECANTH", "PELICAN", matrix='PAM250')
    ('CEELECANTH', '-PELICA--N')


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

    >>> class ATable(SimpleTable):
    ...     x = tables.Float32Col()
    ...     y = tables.Float32Col()
    ...     name = tables.StringCol(16)


    >>> tbl = ATable('test_docs.h5', 'atable1')

    # insert as with pytables. 
    >>> row = tbl.row
    >>> for i in range(50):
    ...    row['x'], row['y'] = i, i * 10
    ...    row['name'] = "name_%i" % i
    ...    row.append()
    >>> tbl.flush()

    #access the entire array via the numpy array interface
    >>> import numpy as np
    >>> np.asarray(tbl)



    #query the data (query() alias of tables' readWhere()
    >>> tbl.query('(x > 4) & (y < 70)') #doctest: +NORMALIZE_WHITESPACE
    array([('name_5', 5.0, 50.0), ('name_6', 6.0, 60.0)],
        dtype=[('name', '|S16'), ('x', '<f4'), ('y', '<f4')])




.. _Skidmarks: http://pypi.python.org/pypi/skidmarks/
.. _SimpleTable: http://pypi.python.org/pypi/simpletable/
.. _nwalign: http://pypi.python.org/pypi/nwalign/
.. _pyfasta: http://pypi.python.org/pypi/pyfasta/
.. _pytables: http://pytables.org/
.. _Needleman-Wunsch: http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm 
