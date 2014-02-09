# -*- coding: utf-8 -*-
"""
Module: tie.py
Created on Mon Jan 27 17:12:19 2014
@author: gav
Description:
Utilities for converting to and from ASA's tie polygon format

The central common format will be geojson records suitable for the
shapely/fiona gis stack.
"""
### Imports
from __future__ import print_function

from pprint import pprint, pformat
import shapely
import geojson
from itertools import islice, izip, chain


### Logging ###
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
debug, info, error = logger.debug, logger.info, logger.error
### Constants
HEADER_LINES = 16
FIELDS = ['id', 'names', 'description', 'type',
          'color_width_height', 'icon_name', 'layer_num',
          'attr1','attr2','attr3','attr4','attr5', 'attr6',
          'link_fn', 'real', 'n_verts']
# Validators
def _str(s):
    return s.strip()

def _int(s):
    return int(_str(s))

def _float(s):
    return float(_str(s))

def _point(s):
    return tuple( float(i) for i in s.strip().split() )

def list_of_ints(s):
    return [ int(i) for i in s.strip().split() ]

FIELD_CONVERTERS = {
    'id'         : _int,
    'names'      : _str,
    'description': _str,
    'type'       : _str,
    'color_width_height': list_of_ints,
    'icon_name'  : _str,
    'layer_num'  : _int,
    'attr1'      : _str,
    'attr2'      : _str,
    'attr3'      : _str,
    'attr4'      : _str,
    'attr5'      : _str,
    'attr6'      : _str,
    'link_fn'    : _str,
    'real'       : _float,
    'n_verts'    : _int,
    'point'      : _point,
}
### Classes


### Functions
def tie_feature(txt):
    """
    Return a complete geojson feature from a single tie record

    Args:

    Need CRS but should probably go in the feature collection level
    """
    lines = chain(txt.split('\n'))
    props = tie_properties(islice(lines, HEADER_LINES))

    coords = tie_coordinates(islice(lines, props['n_verts']))
    polygon = geojson.Polygon(coordinates=coords)
    feat = geojson.Feature(id=props['id'],
                           geometry=polygon,
                           properties=props,
                           crs=tie_crs())
    return feat

def tie_feature_gen(lines, crs=None):
    """
    Return a complete geojson feature from a single tie record

    Need CRS but shoul probably go in the feature collection level
    """
    crs = crs or tie_crs()
    try:
        while True:
            props = tie_properties(islice(lines, HEADER_LINES))
            debug(props['n_verts'])
            coords = tie_coordinates(islice(lines, props['n_verts']))
            polygon = geojson.Polygon(coordinates=coords)
            feat = geojson.Feature(id=props['id'],
                                   geometry=polygon,
                                   properties=props,
                                   crs=crs)
            yield feat
    except StopIteration:
        pass


def tie_properties(prop_lines, conv=FIELD_CONVERTERS):
    """
    Return the geojson properties dictionary from a tie record header
    """
    prop_lines = list(prop_lines)
    if len(prop_lines) < 16:
        raise StopIteration
    props = { k: conv[k](s) for k, s in izip(FIELDS, prop_lines) }
    return props

def tie_coordinates(vert_lines,conv=FIELD_CONVERTERS):
    """
    Return the list of tuple points
    Note the geojson coordinates call for a nested list
    [[(x1, y1), (x2, y2)], # exterior
     [(..), (..)],         # first hole
     [(..), (..)],         # second hole
    ]
    """
    exterior = [ conv['point'](s) for s in vert_lines ]
    verts = [ exterior ]
    return verts

def tie_crs():
    crs = {"crs": {
            "type": "name",
            "properties": {
            "name": "EPSG:4326"}}}
    return crs

def tie_feature_collection(fp):
    """
    Return a geojson feature collection of receptor polygons
    """
    with open(fp) as source:
        features = [ f for f in tie_feature_gen(source) ]
        feat_coll = geojson.FeatureCollection(features)
    return feat_coll

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
### Tests

if __name__ == "__main__":

    print("Done __main__")
