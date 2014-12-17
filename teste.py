#!/usr/bin/python
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

from ctypes import *
from numpy import *
import matrixlib as ml

A = matrix([[10., -1., 2., 0.],
              [-1., 11., -1., 3.],
              [2., -1., 10., -1.],
              [0.0, 3., -1., 8.]])
b = array([6., 25., -11., 15.])


def testaLU(A, b):
    L, U = ml.LUCroutDecompose(A)
    A2 = L.dot(U)
    e = abs(A - A2)
    print "A\n", A
    print "L\n", L
    print "U\n", U
    print "LU\n", A2
    print "A - LU\n", e
    LU = ml.LUCroutInplaceDecompose(A)
    print LU
    x = ml.LUSolve(LU, b)
    print x


def testaGaussJordan(A, b):
    x = ml.GaussJordan(A, b)
    print x


testaLU(A, b)
testaGaussJordan(A, b)
