#! /usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Module: murohy_tie_file.py
Created on Fri Sep 27 13:05:43 2013
@author: gcoombes
Description:

"""
### Imports
from __future__ import print_function

from pprint import pprint

import shapely

### Logging ###
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
debug, info, error = logger.debug, logger.info, logger.error

### Constants
HEADER_LINES = 16
### Classes

### Functions
def read_tie_record(s):
    """Return a tie dict from a complete string record"""
    d = {}
    lines = s.split("\n")
    d['id']          = int(lines[0].strip())
    d['name']        = lines[1].strip()
    d['description'] = lines[2].strip()
    d['type']        = lines[3].strip()
    d['color_width_height'] = [ int(i) for i in lines[4].strip().split() ]
    d['icon_name']   = lines[5].strip()
    d['layer_num']   = int(lines[6].strip())
    d['attr1']       = lines[7].strip()
    d['attr2']       = lines[8].strip()
    d['attr3']       = lines[9].strip()
    d['attr4']       = lines[10].strip()
    d['attr5']       = lines[11].strip()
    d['attr6']       = lines[12].strip()
    d['link_fn']     = lines[13].strip()
    d['real']        = float(lines[14].strip())
    d['n_verts']     = int(lines[15].strip())
    d['verts']       = verts(lines[16: 16 + d['n_verts']])
    return d

def tie_header(lines):
    """Return a tie record (dict)"""
    d = {}
#    for i, line in enumerat e(lines):
#        print(i, line)
    d['id']          = int(lines[0].strip())
    d['name']        = lines[1].strip()
    d['description'] = lines[2].strip()
    d['type']        = lines[3].strip()
    d['color_width_height'] = [ int(i) for i in lines[4].strip().split() ]
    d['icon_name']   = lines[5].strip()
    d['layer_num']   = int(lines[6].strip())
    d['attr1']       = lines[7].strip()
    d['attr2']       = lines[8].strip()
    d['attr3']       = lines[9].strip()
    d['attr4']       = lines[10].strip()
    d['attr5']       = lines[11].strip()
    d['attr6']       = lines[12].strip()
    d['link_fn']     = lines[13].strip()
    d['real']        = float(lines[14].strip())
    d['n_verts']     = int(lines[15].strip())
    d['verts']       = None
    return d

def tie_verts(vert_lines):
    """Return a list of tuples"""
    verts_ls = []
    for line in vert_lines:
        verts_ls.append(tuple(float(i) for i in line.split()))
    return verts_ls

def read_tie_file(fp):
    """Return a list of dicts"""
    records = []
    with open(fp) as f:
        lines = f.readlines()
    while lines:
        header = lines[:HEADER_LINES]
        del lines[:HEADER_LINES]
        record =  (header)
        vs     = lines[:record['n_verts']]
        del lines[:record['n_verts']]
        record['verts'] = tie_verts(vs)
        records.append(record)
    return records

def verts(lines):
    verts_ls = []
    for line in lines:
        verts_ls.append([float(i) for i in line.split()])
    return verts_ls



def tie_record(tie_dict):
    """Return a complete tie record"""
    d = tie_dict
    color_width_height = "{} {} {}".format(*[int(i) for i in d['color_width_height']])
    verts = "\n".join(
        [ " {:16.12f} {:26.13f}".format(*tup) for tup in d['verts']]
    )

    record = """\
 {d[id]}
{d[name]}
{d[description]}
{d[type]}
{color_width_height}
{d[icon_name]}
 {d[layer_num]}
{d[attr1]}
{d[attr2]}
{d[attr3]}
{d[attr4]}
{d[attr5]}
{d[attr6]}
{d[link_fn]}
{d[real]}
 {d[n_verts]}
{verts}""".format(d=tie_dict,
    color_width_height=color_width_height,
    verts=verts)
    return record

def shark_bay_example():
    return """\
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

### Tests

if __name__ == "__main__":
    tie_fp = r"E:\Loc_Data\j0232_murphy\GIS\GISDATA\J0232_receptors.tie"
    new_tie_fp = r"E:\Loc_Data\j0232_murphy\GIS\GISDATA\J0232_split_receptors.tie"
    inner_shark_fp = r"E:\Loc_Data\j0232_murphy\GIS\GISDATA\Inner Shark Bay.tie"
    outer_shark_fp = r"E:\Loc_Data\j0232_murphy\GIS\GISDATA\Outer Shark Bay.tie"
    north_metro_fp = r"E:\Loc_Data\j0232_murphy\GIS\GISDATA\Northern Perth.tie"
    south_metro_fp = r"E:\Loc_Data\j0232_murphy\GIS\GISDATA\Southern Metro.tie"
    recs = read_tie_file(tie_fp)
    inner_recs = read_tie_file(inner_shark_fp)
    outer_recs = read_tie_file(outer_shark_fp)
    north_recs = read_tie_file(north_metro_fp)
    south_recs = read_tie_file(south_metro_fp)

    print(len(recs))

    discard = ['Shark Bay',
               'Yanchep - Mandurah']
    wanted = filter(lambda r: r['name'] not in discard, recs)

    wanted.append(inner_recs[0])
    wanted.append(outer_recs[0])
    wanted.append(north_recs[0])
    wanted.append(south_recs[0])
    for i, r in enumerate(wanted):
        print("({}, '{}'),".format(i+1, r['name']))

    tuples = [
        (14, 'Pelsaert Group'),
        (15, 'Wallabi Group'),
        (16, 'Easter Group'),
        (17, 'Abrolhos Shoals'),
        (1, 'Northern Coast'),
        (2, 'Ningaloo Coast'),
        (3, 'Ningaloo Coast South'),
        (19, 'Inner Shark Bay'),
        (20, 'Outer Shark Bay'),
        (18, 'Shark Bay - Kalbarri'),
        (4, 'Kalbarri - Geraldton'),
        (5, 'Geraldton - Jurien Bay'),
        (6, 'Jurien Bay - Yanchep'),
        (21, 'Northern Metro'),
        (7, 'Rottnest Island'),
        (22, 'Southern Metro'),
        (8, 'Mandurah - Geographe Bay'),
        (9, 'Geographe bay - Augusta'),
        (10, 'Augusta - Walpole'),
        (11, 'Walpole - Albany'),
        (12, 'Recherche Archipelago'),
        (13, 'Albany - Esperance'),

    ]
    new_recs = []
    for i, (idx, name) in enumerate(tuples, 1):
        r = wanted[idx-1]
        assert r['name'] == name, "wrong name: %s -> %s" % (r['name'], name)
        r['id'] = i
        new_recs.append(r)

    for r in new_recs:
        print(r['id'], r['name'])

    with open(new_tie_fp, "w") as sink:
        for r in new_recs:
            sink.write(tie_record(r))

#    s = shark_bay_example()
#    d = read_tie_record(s)
#    pprint.pprint(d)
#    s2 = tie_record(d)
#    print(s2)
#    print(s)
    print("Done __main__")
