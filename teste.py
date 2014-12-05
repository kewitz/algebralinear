#!/usr/bin/python
# -*- coding: utf-8 -*-
from ctypes import *
from numpy import *
import matrixlib as ml

m = matrix(random.rand(3, 3))

L, U = ml.LUCroutDecompose(m)
m2 = L.dot(U)
e = abs(m - m2)

print "A\n", m
print "L\n", L
print "U\n", U
print "LU\n", m2
print "A - LU\n", e

m = ml.LUCroutInplaceDecompose(m)
print m
