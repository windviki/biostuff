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

cdef size_t UP = 1, LEFT = 2, DIAG = 3, NONE = 4

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
def global_align(object _seqj, object _seqi, int match=1, 
                 int gap_open=-1, int gap_extend=-1, object matrix=None):
    """
    perform a global sequence alignment (needleman-wunsch) on seq and and 2. using
    the matrix for nucleotide transition from matrix if available.
    where matrix is of the format provided in the ncbi/data directory.

    >>> from nwalign import global_align
    >>> global_align('COELANCANTH', 'PELICAN')
    ('COELANCANTH', '-PEL-ICAN--')

    """
    if matrix is None:
        return global_align_no_matrix(_seqj, _seqi, match, gap_open, gap_extend)

    cdef char* seqj = _seqj
    cdef char* seqi = _seqi

    cdef size_t max_j = strlen(seqj)
    cdef size_t max_i = strlen(seqi)
    cdef size_t i, j, seqlen, align_counter = 0, p, ib, jb
    cdef int diag_score, up_score, left_score, tscore

    cdef char *align_j, *align_i, *sheader
    cdef char ci, cj
    cdef size_t ii, jj
    cdef PyObject *ai, *aj
    cdef int zero=0, one=1

    assert gap_extend <= 0, "gap_extend penalty must be <= 0"
    assert gap_open <= 0, "gap_open must be <= 0"

    cdef np.ndarray[DTYPE_BOOL, ndim=2] agap = np.ones((max_i + 1, max_j + 1), dtype=np.int8)
    cdef np.ndarray[DTYPE_INT, ndim=2] score = np.empty((max_i + 1, max_j + 1), dtype=np.int)
    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.empty((max_i + 1, max_j + 1), dtype=np.uint)
    cdef np.ndarray[DTYPE_INT, ndim=2] amatrix

    pyheader, amatrix = read_matrix(matrix)
    sheader = pyheader
  
    pointer[0, 0] = NONE
    score[0, 0] = 0
    
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

    score[0, 1:] = gap_open * np.arange(1, max_j + 1, dtype=np.int)
    score[1:, 0] = gap_open * np.arange(1, max_i + 1, dtype=np.int)

    agap[0, 0] = zero

    for i in range(1, max_i + 1):
        ci = seqi[<size_t>(i - 1)]
        ii = strpos(sheader, ci)

        for j in range(1, max_j + 1):
            cj = seqj[<size_t>(j - 1)]
            jj = strpos(sheader, cj)

            if jj == -1 or ii ==  -1:
                sys.stderr.write("'" + chr(ci) + " or " + chr(cj) + " not in scoring matrix\n")
                sys.stderr.write("using score of score of 0\n")
                tscore = 0
            else:
                tscore = amatrix[ii, jj]

            diag_score = score[<size_t>(i - 1), <size_t>(j - 1)] + tscore
            up_score   = score[<size_t>(i - 1), j] + (gap_open if agap[<size_t>(i - 1), j] == zero else gap_extend)
            left_score = score[i, <size_t>(j - 1)] + (gap_open if agap[i, <size_t>(j - 1)] == zero else gap_extend)
            
            if up_score > diag_score:
                if up_score > left_score:
                    score[i, j]  = up_score
                    pointer[i, j] = UP
                else:
                    score[i, j]   = left_score
                    pointer[i, j] = LEFT
            else:
                if left_score > diag_score: # or (diag_score < 0 and i == max_i):
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else:
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap[i, j] = zero# if tscore > 0 else one

    seqlen = max_i + max_j
    ai = PyString_FromStringAndSize(NULL, seqlen)
    aj = PyString_FromStringAndSize(NULL, seqlen)

    cdef int score_max #, = score[:, -1].max()

    # had to use this and PyObject instead of assigning directly...
    align_j = PyString_AS_STRING(aj)
    align_i = PyString_AS_STRING(ai)

    # here, the final pt could be a DIAG, even though it's not 
    # at the highest score... back track to get to the highest
    # score.
    #######################################################
    # so here, given 2 seqs:
    # "AGEBANAN"
    # "ACEBAN"
    # the alignments:
    # AGEBANN
    # ACEBAN--
    # and:
    # AGEBANAN
    # ACEB--AN
    # score equally, but we assume that the former is preferred.
    #######################################################
    if max_j > max_i:
        score_max = score[-1, :].max()
        while score[i, j] < score_max:
            j -= 1
            align_i[align_counter] = c"-"
            align_j[align_counter] = seqj[j]
            align_counter += 1
    elif max_i > max_j:
        score_max = score[:, -1].max()
        while score[i, j] < score_max:
            i -= 1
            align_i[align_counter] = seqi[i]
            align_j[align_counter] = c"-"
            align_counter += 1

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
    return (<object>aj)[::-1], (<object>ai)[::-1] #, score.max()


@cython.boundscheck(False)
@cython.nonecheck(False)
cpdef global_align_no_matrix(object _seqj, object _seqi, int match, int gap_open, int gap_extend):
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
    cdef int diag_score, up_score, left_score

    cdef char *align_j, *align_i
    cdef char ci, cj
    cdef size_t ii, jj
    cdef PyObject *ai, *aj
    cdef int zero=0, one=1

    assert gap_extend <= 0, "gap penalty must be <= 0"
    assert gap_open <= 0, "gap_open must be <= 0"


    cdef np.ndarray[DTYPE_BOOL, ndim=2] agap = np.zeros((max_i + 1, max_j + 1), dtype=np.int8)
    cdef np.ndarray[DTYPE_INT, ndim=2] score = np.zeros((max_i + 1, max_j + 1), dtype=np.int)
    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)

    pointer[0, 0] = NONE
    score[0, 0] = 0
    
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

    score[0, 1:] = gap_open + gap_extend * np.arange(max_j, dtype=np.int)
    score[1:, 0] = gap_open + gap_extend * np.arange(max_i, dtype=np.int)

    agap[0, 1:] = one
    agap[1:, 0] = one
    agap[0, 0] = zero

    for i in range(1, max_i + 1):
        ci = seqi[<size_t>(i - 1)]
        for j in range(1, max_j + 1):
            cj = seqj[<size_t>(j - 1)]
            if cj == ci:
                diag_score = score[i - 1, j - 1] + match
                agap[i, j] = zero
            else:
                diag_score = score[i - 1, j - 1] + (gap_extend if agap[i - 1, j - 1] == one else gap_open)
                agap[i, j] = one

            up_score = score[<size_t>(i - 1), <size_t>j] + (gap_extend if agap[<size_t>(i - 1), j] == one else gap_open)
            left_score   = score[<size_t>i, <size_t>(j - 1)] + (gap_extend if agap[i, <size_t>(j - 1)] == one else gap_open)

            if diag_score >= up_score:
                if diag_score >= left_score:
                    score[<size_t>i, <size_t>j] = diag_score
                    pointer[<size_t>i, <size_t>j] = DIAG
                else:
                    score[<size_t>i, <size_t>j] = left_score
                    pointer[<size_t>i, <size_t>j] = LEFT
            else:
                if up_score >= left_score:
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
    return (<object>aj)[::-1], (<object>ai)[::-1] #, score.max()
