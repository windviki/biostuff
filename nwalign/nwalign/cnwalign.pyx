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
cdef read_matrix(path, dict cache={}):
    if path in cache: return cache[path]
    cdef np.ndarray[DTYPE_INT, ndim=2] a
    #cdef np.ndarray[DTYPE_UINT, ndim=1] header
    cdef size_t ai = 0, i
    cdef int v

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
    result = "".join(headers), a
    cache[path] = result
    return result


cdef inline size_t strpos(char *tstr, char check):
    cdef size_t i = 0
    cdef size_t slen = strlen(tstr)
    while i < slen:
        if tstr[i] == check: return i
        i += 1
    return -1


@cython.boundscheck(False)
@cython.nonecheck(False)
def global_align(object _seqj, object _seqi, int gap=-1, int match=1, int
                 mismatch=-1, int gap_init=-1, object matrix=None):
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

    cdef bint last_match=1
    assert gap <= 0, "gap penalty must be <= 0"
    assert mismatch <= 0, "mismatch must be <= 0"
    assert gap_init <= 0, "gap_init must be <= 0"


    cdef np.ndarray[DTYPE_BOOL, ndim=2] agap = np.zeros((max_i + 1, max_j + 1), dtype=np.int8)
    cdef np.ndarray[DTYPE_INT, ndim=2] score = np.zeros((max_i + 1, max_j + 1), dtype=np.int)
    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)
    cdef np.ndarray[DTYPE_INT, ndim=2] amatrix
    cdef bint has_matrix = 0



    if matrix is not None:
        pyheader, amatrix = read_matrix(matrix)
        sheader = pyheader
        has_matrix = 1
  
    pointer[0, 0] = NONE
    score[0, 0] = 0
    
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

    score[0, 1:] = gap * np.arange(max_j, dtype=np.int)
    score[1:, 0] = gap * np.arange(max_i, dtype=np.int)

    agap[0, 1:] = one
    agap[1:, 0] = one
    agap[0, 0] = zero


    for i in range(1, max_i + 1):
        ci = seqi[<size_t>(i - 1)]
        if has_matrix:
            ii = strpos(sheader, ci)
            if ii == -1:
                py_ci = <object>PyString_FromStringAndSize(&ci, <size_t>1)
                raise Exception(py_ci + "from: " + seqi + " not in scoring matrix")

        for j in range(1, max_j + 1):
            cj = seqj[<size_t>(j - 1)]
            # TODO: move this to separate function.
            if not has_matrix:
                if cj == ci:
                    diag_score = score[i - 1, j - 1] + match
                    last_match = 1
                else:
                    diag_score = score[<size_t>(i - 1), <size_t>(j - 1)] + (gap_init if last_match
                                                        else mismatch)
                    last_match = 0
            else:
                jj = strpos(sheader, cj)
                if jj == -1:
                    py_cj = <object>PyString_FromStringAndSize(&cj, <size_t>1)
                    sys.stderr.write("'" + py_cj + "' from: " + seqj + " not in scoring matrix\n")
                    sys.stderr.write("using score of score of 0\n")
                    tscore = 0
                else:
                    tscore = amatrix[<size_t>ii, <size_t>jj]

                diag_score = score[<size_t>(i - 1), <size_t>(j - 1)] + tscore

            up_score = score[<size_t>(i - 1), <size_t>j] + (gap if agap[<size_t>(i - 1), j] == one else gap_init)
            left_score   = score[<size_t>i, <size_t>(j - 1)] + (gap if agap[i, <size_t>(j - 1)] == one else gap_init)
            
            if diag_score >= up_score:
                if diag_score >= left_score:
                    score[<size_t>i, <size_t>j] = diag_score
                    pointer[<size_t>i, <size_t>j] = DIAG
                    agap[<size_t>i, <size_t>j] = zero
                else:
                    score[<size_t>i, <size_t>j] = left_score
                    pointer[<size_t>i, <size_t>j] = LEFT
                    agap[<size_t>i, <size_t>j] = one
            else:
                agap[i, j] = one
                if up_score > left_score:
                    score[<size_t>i, <size_t>j]  = up_score
                    pointer[<size_t>i, <size_t>j] = UP
                else:
                    score[<size_t>i, <size_t>j]   = left_score
                    pointer[<size_t>i, <size_t>j] = LEFT

    seqlen = max_i + max_j
    ai = PyString_FromStringAndSize(NULL, seqlen)
    aj = PyString_FromStringAndSize(NULL, seqlen)

    # had to use this and PyObject instead of assigning directly...
    align_j = PyString_AS_STRING(aj)
    align_i = PyString_AS_STRING(ai)
        
    p = pointer[i, j]
    while p != NONE:
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
        p = pointer[i, j]

    _PyString_Resize(&aj, align_counter)
    _PyString_Resize(&ai, align_counter)
    return (<object>aj)[::-1], (<object>ai)[::-1]
