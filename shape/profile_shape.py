#! /usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Module: netcdf_to_shape.py
Created on Sun Jul 28 15:51:04 2013
@author: gcoombes
Dexcription: Create a shapefile of contours or probability regions from teh netcdf data

"""
### IMPORTS ###
from __future__ import print_function

import sys, os
import time
from itertools import product

import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, mapping
import fiona

import logging
logging.basicConfig(level=logging.INFO)
debug, info, warning = logging.debug, logging.info, logging.warning

### Constants ###

NC_FN_TPL = "SC8_GEE_DIESEL_136M3_{season}.nc"
SHP_FN_TPL = "j0232_sc8_{season}_{var}.shp"
NC_DIR = r"J:\data\j0232_murphy\netcdf\sc8"
SHP_DIR = r"J:\data\j0232_murphy\gis\gisdata\sc8"
SURF_PROB_STR = "concprob"
ENTR_PROB_STR = "wcmaxprob"
AROM_PROB_STR = "arommaxprob"

CONFIG = {
    "threshold_idx": {
        "surf_prob_1gm":    0,
        "surf_prob_10gm":   1,
        "surf_prob_25gm":   2,
        "entr_prob_10ppb":  0,
        "entr_prob_100ppb": 1,
        "entr_prob_500ppb": 2,
        "arom_prob_6ppb":   0,
        "arom_prob_50ppb":  1,
        "arom_prob_400ppb": 2,
    },
    "nc_var_str": {
        "surf_prob_1gm":    SURF_PROB_STR,
        "surf_prob_10gm":   SURF_PROB_STR,
        "surf_prob_25gm":   SURF_PROB_STR,
        "entr_prob_10ppb":  ENTR_PROB_STR,
        "entr_prob_100ppb": ENTR_PROB_STR,
        "entr_prob_500ppb": ENTR_PROB_STR,
        "arom_prob_6ppb":   AROM_PROB_STR,
        "arom_prob_50ppb":  AROM_PROB_STR,
        "arom_prob_400ppb": AROM_PROB_STR,
    },
}
PROPORTION_TO_PERCENT = 100
### CLASSES ###

### FUNCTIONS ###
def make_shape_file(season, nc_var, cfg):
    """Write a shape file to disk"""
    nc_var_str = cfg['nc_var_str'][nc_var]
    nc_fn = NC_FN_TPL.format(season=season.upper())
    shp_fn = SHP_FN_TPL.format(season=season, var=nc_var)
    nc_fp = os.path.join(NC_DIR, nc_fn)
    shp_fp = os.path.join(SHP_DIR, shp_fn)
    prob_thresholds = prob_threshold_array()
    threshold_idx = cfg['threshold_idx'][nc_var]
    nc_var_str = cfg['nc_var_str'][nc_var]
    print("nc_var_str is {}".format(nc_var_str))

    root = Dataset(nc_fp, 'r')
    dlon, dlat = grid_edge_length(root)
    lon = root.variables['lon']
    lat = root.variables['lat']
    var_array = root.variables[nc_var_str]

    schema = {'geometry':'Polygon',
              'properties': {nc_var_str: 'float'}}
    print("Schema is {}".format(schema))
    driver = "ESRI Shapefile"

    with fiona.collection(shp_fp, 'w', driver, schema) as sink:
        for lower_prop, upper_prop in prob_thresholds[1:]:  # discard first bin
            prob_percent = lower_prop * PROPORTION_TO_PERCENT
            threshold_arr = var_array[threshold_idx, :, :]
            log_arr = np.logical_and(threshold_arr >= lower_prop,
                                     threshold_arr <  upper_prop)
            idx_arr = np.argwhere(log_arr)
            for coord_idx in idx_arr:
                coord = (lon[coord_idx[1]], lat[coord_idx[0]])
                polygon = grid_polygon(coord, dlon, dlat)
                record = {'geometry': mapping(polygon),
                          'properties': {nc_var_str: prob_percent }}
                sink.write(record)
    del root

def mk_shp_v2(season, nc_var, cfg):
    """Write a shape file to disk"""
    nc_var_str = cfg['nc_var_str'][nc_var]
    nc_fn = NC_FN_TPL.format(season=season.upper())
    shp_fn = SHP_FN_TPL.format(season=season, var=nc_var)
    nc_fp = os.path.join(NC_DIR, nc_fn)
    shp_fp = os.path.join(SHP_DIR, shp_fn)
    prob_thresholds = prob_threshold_array()
    threshold_idx = cfg['threshold_idx'][nc_var]
    nc_var_str = cfg['nc_var_str'][nc_var]
    print("nc_var_str is {}".format(nc_var_str))

    root = Dataset(nc_fp, 'r')
    lon = root.variables['lon']
    lat = root.variables['lat']
    dlon = lon[1] - lon[0]
    dlat = lat[1] - lat[0]
    var_array = root.variables[nc_var_str]
    threshold_arr = var_array[threshold_idx, :, :]

    schema = {'geometry':'Polygon',
              'properties': {nc_var_str: 'float'}}
    print("Schema is {}".format(schema))
    driver = "ESRI Shapefile"

    with fiona.collection(shp_fp, 'w', driver, schema) as sink:
        for lower_prop, upper_prop in prob_thresholds[1:]:  # discard first bin
            prob_percent = lower_prop * PROPORTION_TO_PERCENT

            log_arr = np.logical_and(threshold_arr >= lower_prop,
                                     threshold_arr <  upper_prop)
            idx_arr = np.argwhere(log_arr)
            for lat_idx, lon_idx in idx_arr:
                coord = (lon[lon_idx], lat[lat_idx])
                polygon = grid_polygon(coord, dlon, dlat)
                record = {'geometry': mapping(polygon),
                          'properties': {nc_var_str: prob_percent }}
                sink.write(record)
    del root

def mk_shp_v3(season, nc_var, cfg):
    """Write a shape file to disk"""
    nc_var_str = cfg['nc_var_str'][nc_var]
    nc_fn = NC_FN_TPL.format(season=season.upper())
    shp_fn = SHP_FN_TPL.format(season=season, var=nc_var)
    nc_fp = os.path.join(NC_DIR, nc_fn)
    shp_fp = os.path.join(SHP_DIR, shp_fn)
    prob_thresholds = prob_threshold_array()
    threshold_idx = cfg['threshold_idx'][nc_var]
    nc_var_str = cfg['nc_var_str'][nc_var]
    print("nc_var_str is {}".format(nc_var_str))

    root = Dataset(nc_fp, 'r')
    lon = np.array(root.variables['lon'])
    lat = np.array(root.variables['lat'])
    dlon = lon[1] - lon[0]
    dlat = lat[1] - lat[0]
    var_array = root.variables[nc_var_str]
    threshold_arr = var_array[threshold_idx, :, :]

    schema = {'geometry':'Polygon',
              'properties': {nc_var_str: 'float'}}
    print("Schema is {}".format(schema))
    driver = "ESRI Shapefile"

    with fiona.collection(shp_fp, 'w', driver, schema) as sink:
        for lower_prop, upper_prop in prob_thresholds[1:]:  # discard first bin
            prob_percent = lower_prop * PROPORTION_TO_PERCENT

            log_arr = np.logical_and(threshold_arr >= lower_prop,
                                     threshold_arr <  upper_prop)
            lat_idx, lon_idx = np.nonzero(log_arr)
            lons = lon[lon_idx].reshape(-1,1)
            lats = lat[lat_idx].reshape(-1,1)

            for lon_datum, lat_datum in np.hstack((lons, lats)):
                polygon = grid_polygon((lon_datum, lat_datum), dlon, dlat)
                record = {'geometry': mapping(polygon),
                          'properties': {nc_var_str: prob_percent }}
                sink.write(record)
    del root

def grid_edge_length(root):
    """Return a tuple of (delta_lon, delta_lat)
    """
    delta_lon = np.mean(np.diff(root.variables["lon"]))
    delta_lat = np.mean(np.diff(root.variables["lat"]))
    return (delta_lon, delta_lat)

def grid_polygon(coord, dlon, dlat):
    """Given a lon, lat point, construct a grid square polygon"""
    lonc, latc = coord
    lon1, lon2 = lonc - dlon / 2, lonc + dlon / 2
    lat1, lat2 = latc - dlat / 2, latc + dlat / 2
    coords = [(lon1, lat1),
              (lon2, lat1),
              (lon2, lat2),
              (lon1, lat2),
              (lon1, lat1)]
    return Polygon(coords)

def dot(c="."):
    sys.stdout.write(c)
    sys.stdout.flush()

def prob_threshold_array(step=0.02):
    """Return a array of tuples

    eg tuple (0.02, 0.04)
    """
    prob_vector = np.arange(0., 1 + 1.5 * step, step) # Including end points
    prob_thresholds = np.zeros((len(prob_vector) -1, 2))
    prob_thresholds[:,0] = prob_vector[:-1]
    prob_thresholds[:,1] = prob_vector[1:]
    return prob_thresholds


### TESTS ###

if __name__ == "__main__":
    nc_vars = ["surf_prob_1gm", "surf_prob_10gm", "surf_prob_25gm",
               "entr_prob_10ppb", "entr_prob_100ppb", "entr_prob_500ppb",
               "arom_prob_6ppb",  "arom_prob_50ppb", "arom_prob_400ppb"]
    seasons = ["sum", "win"]
    tuples = product(seasons, nc_vars)
    start_time = lap_time = time.time()
#    for season, nc_var in tuples:
    for season, nc_var in [("sum", "entr_prob_100ppb")]:
        print("*** {} {} ***".format(season, nc_var))
        make_shape_file(season, nc_var, CONFIG)

    elapsed = time.time() - start_time
    print("Full run in {} seconds".format(elapsed))



    print("Done __main__")

