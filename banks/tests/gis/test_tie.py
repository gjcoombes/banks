#1 /usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Module: test_tie.py
Package: banks.git
Created on Sun Feb 09 09:59:21 2014
@author: gav
Description:

"""
### Imports
from __future__ import print_function, division

import pytest
from pprint import pprint, pformat
from itertools import islice
from context import banks
from banks.gis.tie import *
#from banks.gis.tie_file_utils import *
### Logging
import logging
logging.basicConfig(level=logging.DEBUG)
debug, info, error = logging.debug, logging.info, logging.error
### Constants

### Classes

### Fixtures
#  Available from upward contest.py
# for example in banks.test.conftest.py
# * shark_bay - the text from a single tie record
### Tests

def test_shark_bay_fix(shark_bay_record):
    print(shark_bay_record)
    assert 1 == 1

def test_shark_bay_tie(shark_bay_tie_filepath):
    with open(shark_bay_tie_filepath) as fh:
        data = fh.read()
        print(data)
    assert 1 == 1

def test_j0232_receptor_tie(j0232_receptor_tie_filepath):
    with open(j0232_receptor_tie_filepath) as fh:
        data = fh.read()
        print(data)
    assert 1 == 1

def test_shark_bay_record_generator(shark_bay_record_generator):
    gen = shark_bay_record_generator()
    map(print, islice(gen, 16))
    map(print, islice(gen, 15))
    assert 1 == 1

def test_tie_properties(shark_bay_tie_filepath):
    with open(shark_bay_tie_filepath) as source:
        props = tie_properties(islice(source, 16))
    assert 'n_verts' in props

def test_tie_feature(shark_bay_record):
    feat = tie_feature(shark_bay_record)
    assert feat['properties']['names'] == 'Shark Bay - Kalbarri'
    assert len(feat['geometry']['coordinates'][0]) == 15

def test_tie_feature_collection(shark_bay_tie_filepath):
    feats = tie_feature_collection(shark_bay_tie_filepath)
    pprint(feats)
    assert len(feats['features']) == 1

def test_tie_feature_collection2(j0232_receptor_tie_filepath):
    feats = tie_feature_collection(j0232_receptor_tie_filepath)
    pprint(feats)
    assert len(feats['features']) == 20

if __name__ == "__main__":
    pytest.main(['-x'])
    print("Done __main__")
