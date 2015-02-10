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

from numpy import *
import unittest
import matrixlib as ml
import numpy.testing


def getAB():
    A = matrix([[10., -1., 2., 0.],
                [-1., 11., -1., 3.],
                [2., -1., 10., -1.],
                [0.0, 3., -1., 8.]])
    b = array([6., 25., -11., 15.])
    return A, b


class LUTest(unittest.TestCase):
    def setUp(self):
        self.A, self.b = getAB()

    def testLUeA(self):
        L, U = ml.LUCroutDecompose(self.A)
        A2 = L.dot(U)
        numpy.testing.assert_almost_equal(self.A, A2)

    def testSolver(self):
        LU = ml.LUCroutInplaceDecompose(self.A)
        x = ml.LUSolve(LU, self.b)
        numpy.testing.assert_almost_equal(x, [1., 2., -1., 1.])


class GaussJordanTest(unittest.TestCase):
    def setUp(self):
        self.A, self.b = getAB()

    def test(self):
        x = ml.GaussJordan(self.A, self.b)
        numpy.testing.assert_almost_equal(x, [1., 2., -1., 1.])


class JacobiTest(unittest.TestCase):
    def setUp(self):
        self.A, self.b = getAB()

    def test(self):
        x = ml.Jacobi(self.A, self.b, 20)
        numpy.testing.assert_almost_equal(x, [1., 2., -1., 1.])


class GaussSeidelTest(unittest.TestCase):
    def setUp(self):
        self.A, self.b = getAB()

    def test(self):
        x = ml.GaussSeidel(self.A, self.b, 10)
        numpy.testing.assert_almost_equal(x, [1., 2., -1., 1.])

if __name__ == "__main__":
    unittest.main()