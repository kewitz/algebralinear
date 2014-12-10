# -*- coding: utf-8 -*-
import numpy as np
from ctypes import *
from numpy import *

try:
    lib = cdll.LoadLibrary('./src/lib.so')
except e:
    raise e


def LUCroutDecompose(A):
    """
    Implementação do método de Crout para decomposição LU.
    """
    assert A.shape[0] == A.shape[1] and type(A) is matrix, "'A' deve ser NxN."
    n = A.shape[0]
    L = zeros(A.shape)
    U = L.copy()
    lib.LUDec(n, byref(np.ctypeslib.as_ctypes(A)),
              byref(np.ctypeslib.as_ctypes(L)),
              byref(np.ctypeslib.as_ctypes(U)))
    return L, U


def LUCroutInplaceDecompose(A):
    """
    Implementação do método de Crout para decomposição LU sobrescrevendo a
    matriz original.
    """
    assert A.shape[0] == A.shape[1] and type(A) is matrix, "'A' deve ser NxN."
    n = A.shape[0]
    lib.LUInDec(n, byref(np.ctypeslib.as_ctypes(A)))
    return A


def GaussJordan(A, b):
    """
    Resolve Ax = b através do método de eliminação de Gauss-Jordan com pivotação.
    """
    assert A.shape[0] == A.shape[1] and type(A) is matrix, "'A' deve ser NxN."
    assert b.shape[0] == A.shape[0], "'b' deve ser compatível com A."
    n = A.shape[0]
    x = zeros(n)
    try:
        lib.SolveGJ(n, byref(ctypeslib.as_ctypes(A)),
                    byref(ctypeslib.as_ctypes(x)),
                    byref(ctypeslib.as_ctypes(b)))
    except e:
        raise e
    return x
