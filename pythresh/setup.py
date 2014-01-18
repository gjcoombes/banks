#! /usr/bin/env python2
# -*- coding: utf-8 -*-
# distutils: include_dirs = C:/WinPython-64bit-2.7/python-2.7.5.amd64/Lib/site-packages/numpy/core/include
"""
Module: setup.py
Created on Sun Oct 13 16:00:50 2013
@author: gav
Description:

"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Compiler import Options
from Cython.Build import cythonize

import numpy

Options.annotate = True

copt = ['/openmp', '/Ox', '/fp:fast','/favor:INTEL64','/Og']

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
        Extension("cythresh_timeseries",
                  sources = ["cythresh_v0_8.pyx"],
                  include_dirs = [
                      "C:/WinPython-64bit-2.7/python-2.7.5.amd64/Lib/site-packages/numpy/core/include",
                      "C:/Program Files (x86)/Microsoft Visual Studio 9.0/VC/include",
                  ],
                  extra_compile_args = copt,
        ),
    ],
)
