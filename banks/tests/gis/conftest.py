#1 /usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Module: conftest.py
Package: banks.tests
Created on Sun Feb 09 10:14:00 2014
@author: gav
Description:
This conftest file contains fixtures for any test_modules within subdirectories
"""
### Imports
from __future__ import print_function, division

import os.path as osp
import pytest

from context import banks

### Logging
import logging
logging.basicConfig(level=logging.DEBUG)
debug, info, error = logging.debug, logging.info, logging.error
### Constants

### Classes

### Functions
@pytest.fixture(scope='module')
def shark_bay_record():
    txt = """\
 21
Shark Bay - Kalbarri

POLYGON
1 1 1
RECTANGLE
 1







0.0
 15
 113.827202986249           -26.917153608139
 113.738670759538           -26.9654503372097
 113.815310597586           -27.0855132986786
 113.865522905273           -27.1619569843129
 113.919699342514           -27.2465720170934
 113.973875779755           -27.3534240082475
 114.018802581369           -27.4449294701931
 114.034659099586           -27.4754144470855
 114.063729382984           -27.5551045039212
 114.090156913345           -27.6265420015409
 114.086192783791           -27.6710205043739
 114.198509787827           -27.668679981573
 114.039944605658           -27.2512709681449
 113.827202986249           -26.9183318232656
 113.827202986249           -26.917153608139"""
    return txt

@pytest.fixture(scope='module')
def shark_bay_tie_filepath():
    base_dir = osp.abspath(osp.dirname(banks.__file__))
    fn = "Shark Bay-Kalbarri.tie"
    fp = osp.join(base_dir, 'tests', 'resources', fn)
    return fp

@pytest.fixture(scope='module')
def shark_bay_record_generator():
    base_dir = osp.abspath(osp.dirname(banks.__file__))
    fn = "Shark Bay-Kalbarri.tie"
    fp = osp.join(base_dir, 'tests', 'resources', fn)
    def gen():
        for line in open(fp):
            yield line
    return gen

@pytest.fixture(scope='module')
def j0232_receptor_tie_filepath():
    base_dir = osp.abspath(osp.dirname(banks.__file__))
    fn = "J0232_receptors.tie"
    fp = osp.join(base_dir, 'tests', 'resources', fn)
    return fp

### Tests

if __name__ == "__main__":
    shark_bay_tie_filepath()





    print("Done __main__")
