# -*- coding: utf-8 -*-
"""
The MIT License (MIT)

Copyright (c) 2014 Leonardo Kewitz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
__version__ = "0.1.00"
__author__ = "kewitz"
__license__ = "MIT"
DEBUG = True

from ctypes import cdll, byref
from numpy import zeros, matrix, ctypeslib


def debug(s, pref="%"):
    if DEBUG:
        print "[%s] %s" % (pref, s)

try:
    lib = cdll.LoadLibrary('./src/lib.so')
except:
    raise
finally:
    CUDAdevices = lib.GetCUDADevices()
    CUDAcapable = CUDAdevices >= 1
    if CUDAcapable:
        debug("Hardware compatível com CUDA detectado.")
        kernels = cdll.LoadLibrary('./src/kernels.so')


def LUCroutDecompose(A):
    """
    Implementação do método de Crout para decomposição LU.
    """
    assert A.shape[0] == A.shape[1] and type(A) is matrix, "'A' deve ser NxN."
    L = zeros(A.shape)
    n = A.shape[0]
    U = L.copy()
    lib.LUDec(n, byref(ctypeslib.as_ctypes(A)),
              byref(ctypeslib.as_ctypes(L)),
              byref(ctypeslib.as_ctypes(U)))
    return L, U


def LUCroutInplaceDecompose(A):
    """
    Implementação do método de Crout para decomposição LU sobrescrevendo a
    matriz original.
    """
    assert A.shape[0] == A.shape[1] and type(A) is matrix, "'A' deve ser NxN."
    LU = A.copy()
    n = A.shape[0]
    lib.LUInDec(n, byref(ctypeslib.as_ctypes(LU)))
    return LU


def LUSolve(LU, b):
    """
    Resolve LU = b.
    """
    n = LU.shape[0]
    x = b.copy()
    lib.LUSolve(n, byref(ctypeslib.as_ctypes(LU)),
                byref(ctypeslib.as_ctypes(x)),
                byref(ctypeslib.as_ctypes(b)))
    return x


def GaussJordan(A, b):
    """
    Resolve Ax = b através do método de eliminação de Gauss-Jordan com
    pivotação.
    """
    assert A.shape[0] == A.shape[1] and type(A) is matrix, "'A' deve ser NxN."
    assert b.shape[0] == A.shape[0], "'b' deve ser compatível com A."
    n = A.shape[0]
    x = zeros(n)
    try:
        lib.SolveGJ(n, byref(ctypeslib.as_ctypes(A)),
                    byref(ctypeslib.as_ctypes(x)),
                    byref(ctypeslib.as_ctypes(b)))
    except:
        raise
    return x


def Jacobi(A, b, ks=1000, cuda=CUDAcapable):
    """
    Resolve Ax = b através do método iterativo de Jacobi.
    """
    assert A.shape[0] == A.shape[1] and type(A) is matrix, "'A' deve ser NxN."
    assert b.shape[0] == A.shape[0], "'b' deve ser compatível com A."
    n = A.shape[0]
    x = zeros(n)
    func = kernels.CUJacobi if cuda else lib.SolveJacobi
    try:
        func(n, ks, byref(ctypeslib.as_ctypes(A)),
             byref(ctypeslib.as_ctypes(x)),
             byref(ctypeslib.as_ctypes(b)))
    except:
        raise
    return x


def GaussSeidel(A, b, ks=50):
    """
    Resolve Ax = b através do método iterativo de Jacobi.
    """
    assert A.shape[0] == A.shape[1] and type(A) is matrix, "'A' deve ser NxN."
    assert b.shape[0] == A.shape[0], "'b' deve ser compatível com A."
    n = A.shape[0]
    x = zeros(n)
    try:
        lib.SolveGaussSeidel(n, ks, byref(ctypeslib.as_ctypes(A)),
                             byref(ctypeslib.as_ctypes(x)),
                             byref(ctypeslib.as_ctypes(b)))
    except:
        raise
    return x
