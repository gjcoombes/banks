#! /usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Module: grid_search.py
Created on Mon Dec 02 14:49:54 2013
@author: gcoombes
Description:
Catalogue all the depth grids on the shared file system

1. Find each depth grid.
2. Read and parse the header
    err - write filepath to grid_error_fps.log
3. Compute the required fields
4. Sink to database

"""
### Imports
from __future__ import print_function

import os
import os.path as osp
import fnmatch

import psycopg2 as pg
import shapely
import pyproj

from bonaparte.utils.general import checksum

### Logging ###
import logging
logging.basicConfig(level=logging.DEBUG)
debug, info, error = logging.debug, logging.info, logging.error

### Constants

### Classes

### Functions
def insert_record(record_d, conn=None,
                  conn_str="host=tornado dbname=joseph user=postgres password=Please123"):
    """Insert a grid record to the db"""
    _conn = conn or pg.connect(conn_str)
    cur = _conn.cursor()
    ins_str = """
    INSERT INTO simap_grids
    (name, description, abbreviation, filename, source_dirs, sha1,
     project, origin, delta, extent, habitat_id, scenario_id)
    VALUES
    (%(name)s, %(description)s, %(abbreviation)s, %(filename)s, %(source_dirs)s,
    %(sha1)s, %(project)s, %(origin)s, %(delta)s, %(extent)s,
    %(habitat_id)s, %(scenario_id)s);
    """
    cur.execute(ins_str, record_d)
    _conn.commit()
    cur.close()

def grid_spider(base_dirs=None):
    """Yield paths to grid files"""
    if base_dirs is None:
        _base_dirs = [r"\\manta\e$\Loc_Data"]
    else:
        _base_dirs = list(base_dirs)
    patt = "*.dep"

    for bd in _base_dirs:
        for root, dirs, files in os.walk(bd):
            for f in fnmatch.filter(files, patt):
                yield osp.join(root, f)

### Fixtures
def connection():
    conn_str = "host=tornado dbname=db user=postgres password=Please123"
    conn = pg.connect(conn_str)
    return conn

### Tests
def test_imports():
    assert pg
    assert gluon

def test_connection():
    conn_str = "host=tornado dbname=joseph user=postgres password=Please123"
    conn = pg.connect(conn_str)
    assert conn
    conn.close()



def test_insert_table():
    conn_str = "host=tornado dbname=joseph user=postgres password=Please123"
    conn = pg.connect(conn_str)
    cur = conn.cursor()
#    tup = (3,murphy,davo,new record,,,,"105, -24, 125, -12",,)
    cur.execute("""insert into grid (name, abbreviation, description, path, sha1, project, extent, origin, delta)
    values (%s, %s, %s, %s, %s, %s, %s, %s, %s); """,
     ('woodside_trial2', 'wel', 'new record', ' ', ' ', ' ', '105, -24, 125, -12', ' ', ' '))
    conn.commit()
    cur.close()
    conn.close()

def test_insert_record():
    """Sample record for the simap_grids table"""
    d = {
        "name": "Balnaves leak",
        "description": """ Slow sub surface leak """,
        "abbreviation": "balnaves_800m",
        "filename": "MidWA_NWS_1km.DEP",
        "source_dirs": [
            r"\\manta\e$\Loc_Data\j0275_apache_eastspar\GRIDS",
            r"\\apasa-server\p$\SIMAPGRIDS\GeneralisedNW_WA_Setups"
        ],
        "project": "apache balnaves",
        "sha1": "abbcceef",
        "origin": "POINT(114.0 -32)",
        "delta": "POINT(0.007 0.006)",
        "extent": None,
        "habitat_id": '13',
        "scenario_id": '65',
    }
    insert_record(d)
    return d

def test_grid_table():
    conn_str = "host=tornado dbname=joseph user=postgres password=Please123"
    conn = pg.connect(conn_str)
    cur = conn.cursor()
    q = """ SELECT * FROM grid; """
    cur.execute(q)
    result = cur.fetchall()
    for r in result:
        print(r)
    cur.close()
    conn.close()

def test_spider():
    base_dir = r"\\READYNAS_DUO\backup\Project_Archive"
    for grid_fp in grid_spider(base_dir ):
        print(grid_fp)


if __name__ == "__main__":

    test_imports()
#    test_connection()

#    test_DAL_connection()
#    test_DAL_grid_table()
#    test_DAL_person_table()
#    test_insert_table()
#    test_grid_table()
    test_spider()





























    print("Done __main__")
