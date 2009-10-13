"""
------------------------------------------------------------------------------
nwalign: fast `cython`_  - `Needleman-Wunsch`_ alignment
------------------------------------------------------------------------------

.. _`Needleman-Wunsch`: http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm 
.. _`scoring matrix`: http://en.wikipedia.org/wiki/Substitution_matrix
.. _`cython`: http://cython.org

This module provides a python module and a command-line interface to do global-
sequence alignment using the `Needleman-Wunsch` algorithm. It uses `cython`_ 
and numpy for speed.

Example Command-Line Usage 
==========================


the nwalign executable is installed to the PATH by setuptools::

    $ nwalign alphabet alpet
    alphabet
    alp---et

specify an alignment `scoring matrix`_ ::

    $ nwalign --matrix /usr/share/ncbi/data/BLOSUM62 EEAEE EEEEG
    EEAEE-
    EE-EEG


Usage as a python module
========================
::

    >>> import nwalign as nw
    >>> nw.global_align("CEELECANTH", "PELICAN", matrix='PAM250')
    ('CEELECANTH', '-PELICA--N')


the matrix is specified as the full path to an `scoring matrix`_ as
is distributed with the NCBI toolset.
"""
import numpy as np
cimport numpy as np


cimport cython
import sys
import os.path


cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char *s)

cdef extern from "Python.h":
    ctypedef void PyObject
    PyObject *PyString_FromStringAndSize(char *, size_t)
    int _PyString_Resize(PyObject **, size_t)
    char * PyString_AS_STRING(PyObject *)

ctypedef np.int_t DTYPE_INT
ctypedef np.uint_t DTYPE_UINT
ctypedef np.int8_t DTYPE_BOOL

cdef size_t UP = 0, LEFT = 1, DIAG = 2, NONE = 3

MATRIX = ['BLOSUM62']

cdef inline int imax2(int a, int b):
    if a >= b: return a
    return b


@cython.boundscheck(False)
def read_matrix(path):
    cdef np.ndarray[DTYPE_INT, ndim=2] a
    cdef size_t ai = 0, i
    cdef int v

    #if os.path.basename(path) in MATRIX:
    #    import cPickle
    #    op = os.path
    #    fpath = op.join(op.dirname(op.abspath(__file__)), 'data', path + '.pkl')
    #    fh = open(fpath, 'rb')
    #    a = cPickle.load(fh)
    #    h = cPickle.load(fh)
    #    fh.close()
    #    return h, a

    fh = open(path)
    headers = None
    while headers is None:
        line = fh.readline().strip()
        if line[0] == '#': continue
        headers = [x for x in line.split(' ') if x]

    line = fh.readline()
    a = np.ndarray((len(headers), len(headers)), dtype=np.int)

    while line:
        line = [int(x) for x in line[:-1].split(' ')[1:] if x]
        for i in range(len(line)):
            v = line[i]
            a[ai, i] = v
        ai += 1
        line = fh.readline()
    assert ai == len(headers), (ai, len(headers))
    return "".join(headers), a

cdef inline size_t strpos(char *tstr, char check):
    cdef size_t i = 0
    cdef size_t slen = strlen(tstr)
    while i < slen:
        if tstr[i] == check: return i
        i += 1
    return -1


@cython.boundscheck(False)
cpdef global_align(object _seqj, object _seqi, int gap=-1, int match=1, int mismatch=-1, int gap_init=-1, object matrix=None):
    """
    perform a global sequence alignment (needleman-wunsch) on seq and and 2. using
    the matrix for nucleotide transition from matrix if available.
    where matrix is of the format provided in the ncbi/data directory.

    >>> from nwalign import global_align
    >>> global_align('COELANCANTH', 'PELICAN')
    ('COELANCANTH', '-PEL-ICAN--')

    """
    cdef char* seqj = _seqj
    cdef char* seqi = _seqi

    cdef size_t max_j = strlen(seqj)
    cdef size_t max_i = strlen(seqi)
    cdef size_t i, j, seqlen, align_counter = 0, p
    cdef int diag_score, up_score, left_score, tscore

    cdef char *align_j, *align_i, *sheader
    cdef char ci, cj
    cdef size_t ii, jj
    cdef PyObject *ai, *aj
    cdef int zero=0, one=1


    cdef np.ndarray[DTYPE_BOOL, ndim=2] agap = np.zeros((max_i + 1, max_j + 1), dtype=np.int8)
    cdef np.ndarray[DTYPE_INT, ndim=2] score = np.zeros((max_i + 1, max_j + 1), dtype=np.int)
    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)
    cdef np.ndarray[DTYPE_INT, ndim=2] amatrix



    if matrix is not None:
        pyheader, amatrix = read_matrix(matrix)
        sheader = pyheader
  
    pointer[<size_t>0, <size_t>0] = NONE
    score[<size_t>0, <size_t>0] = 0
    
    pointer[<size_t>0, <size_t>1:] = LEFT
    pointer[<size_t>1:, <size_t>0] = UP

    score[<size_t>0, <size_t>1:] = gap * np.arange(max_j, dtype=np.int)
    score[<size_t>1:, <size_t>0] = gap * np.arange(max_i, dtype=np.int)

    agap[<size_t>0, <size_t>1:] = one
    agap[<size_t>1:, <size_t>0] = one
    agap[<size_t>0, <size_t>0] = zero

    for i in range(1, max_i + 1):
        ci = seqi[i - 1]
        for j in range(1, max_j + 1):
            cj = seqj[j - 1]

            if matrix is None:
                diag_score = score[i - 1, j - 1] + (cj == ci and match or mismatch)
            else:
                ii = strpos(sheader, ci)
                jj = strpos(sheader, cj)
                if ii == -1 or jj == -1:
                    if ii == -1:
                        py_ci = <object>PyString_FromStringAndSize(&ci, <size_t>1)
                        #raise Exception(py_ci + "from: " + seqi  + " not in scoring matrix")
                        sys.stderr.write("'" + py_ci + "' from: " + seqi  + " not in scoring matrix\n")
                    if jj == -1:
                        py_cj = <object>PyString_FromStringAndSize(&cj, <size_t>1)
                        #raise Exception(py_cj + "from: " + seqj + " not in scoring matrix")
                        sys.stderr.write("'" + py_cj + "' from: " + seqj + " not in scoring matrix\n")
                    sys.stderr.write("using score of score of 0\n")
                    tscore = 0
                else:
                    tscore = amatrix[ii, jj]

                diag_score = score[i - 1, j - 1] + tscore

            up_score = score[i - 1, j] + (gap if agap[i - 1, j] == one else gap_init)
            left_score   = score[i, j - 1] + (gap if agap[i, j - 1] == one else gap_init)
            
            if diag_score >= up_score:
                if diag_score >= left_score:
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap[i, j] = zero
                else:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                    agap[i, j] = one
            else:
                agap[i, j] = one
                if up_score > left_score:
                    score[i, j]  = up_score
                    pointer[i, j] = UP
                else:
                    score[i, j]   = left_score
                    pointer[i, j] = LEFT

    seqlen = max_i + max_j
    ai = PyString_FromStringAndSize(NULL, seqlen)
    aj = PyString_FromStringAndSize(NULL, seqlen)

    # had to use this and PyObject instead of assigning directly...
    align_j = PyString_AS_STRING(aj)
    align_i = PyString_AS_STRING(ai)
        
    while True:
        p = pointer[i, j]
        if p == NONE: break
        if p == DIAG:
            i -= 1
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = seqi[i]
        elif p == LEFT:
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = c"-"
        elif p == UP:
            i -= 1
            align_j[align_counter] = c"-"
            align_i[align_counter] = seqi[i]
        else:
            raise Exception('wtf!')
        align_counter += 1

    _PyString_Resize(&aj, align_counter)
    _PyString_Resize(&ai, align_counter)
    return (<object>aj)[::-1], (<object>ai)[::-1]
            
def main():
    import optparse
    parser = optparse.OptionParser(usage="""
    %prog [options] seq1 seq2 
    """)
    parser.add_option("--gap", dest="gap", help="gap extend penalty (must be integer <= 0)", type="int", default=-1)
    parser.add_option("--gap_init", dest="gap_init", help="gap start penalty (must be integer <= 0)", type="int", default=-1)
    parser.add_option("--match", dest="match", help="match score (must be integer > 0)", type="int", default=1)
    parser.add_option("--mismatch", dest="mismatch", help="gap penalty (must be integer < 0)", type="int", default=-1)
    parser.add_option("--matrix", dest="matrix", help="scoring matrix in ncbi/data/ format,\
                                      if not specificied, match/mismatch are used", default=None)

    try:
        options, args = parser.parse_args()
    except:
        sys.exit(parser.print_help())

    if len(args) != 2:
        sys.exit(parser.print_help())
    print "\n".join(global_align(args[0], args[1], options.gap, options.match, options.mismatch, options.gap_init, options.matrix))
        
