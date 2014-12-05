# -*- coding: utf-8 -*-
import numpy as np
from ctypes import *
from numpy import *

try:
    lib = cdll.LoadLibrary('./matrix.so')
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
    lib.LUCroutDecompose(n, byref(np.ctypeslib.as_ctypes(A)),
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
    lib.LUCroutInplaceDecompose(n, byref(np.ctypeslib.as_ctypes(A)))
    return A
