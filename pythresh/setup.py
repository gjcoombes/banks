#! /usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Module: setup.py
Created on Sun Oct 13 16:00:50 2013
@author: gav
Description:

"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("cyshape", ["cyshape.pyx"])]
)