#! /usr/bin/env python2
# cython: profile=False
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
cimport numpy as np
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
#PROPORTION_TO_PERCENT = 100
### CLASSES ###

### FUNCTIONS ###
def mk_shp_v4(season, nc_var, cfg):
    """Write a shape file to disk"""
    cdef:
        float [:] lon_corners
        float [:] lat_corners
        float [:, :] prob_bins, prob_arr
        int threshold_idx

    nc_fp = os.path.join(NC_DIR, NC_FN_TPL.format(season=season.upper()))
    shp_fp = os.path.join(SHP_DIR, SHP_FN_TPL.format(season=season, var=nc_var))
    prob_bins = prob_threshold_array()
    threshold_idx = cfg['threshold_idx'][nc_var]
    nc_var_str = cfg['nc_var_str'][nc_var]

    root = Dataset(nc_fp, 'r')
    lon_corners = compute_corners(root.variables['lon'])
    lat_corners = compute_corners(root.variables['lat'])
    prob_arr = root.variables[nc_var_str][threshold_idx, :, :]

    schema = {'geometry':'Polygon',
              'properties': {nc_var_str: 'float'}}
    driver = "ESRI Shapefile"

    with fiona.collection(shp_fp, 'w', driver, schema) as sink:
        compute_and_write_records(sink, prob_bins, prob_arr,
                                  lon_corners, lat_corners, nc_var_str)
    del root

cdef float [:] compute_corners(float [:] vec):
    """Return the corner values of each grid cell

    A vector shifted by half a step plus one element
    """
    cdef int n = vec.shape[0]
    cdef int i 
    cdef float delta = (vec[1] - vec[0])
    cdef float [:] corners = np.zeros((vec.shape[0] + 1))

    for i in range(n):
        corners[i] = vec[i] - delta
    corners[n] = vec[n-1] + delta
    return corners

cdef compute_and_write_records(
        sink, 
        float [:, :] prob_bins, 
        prob_arr,
        float [:] lon,
        float [:] lat,
        nc_var_str):
    """Write a polygon for each grid cell by bins
    """
    cdef:
        int [:, :] log_arr
        int i, j, lon_idx, lat_idx, PROPORTION_TO_PERCENT
        int [:] lon_ixs, lat_ixs 
        float lower_prob, upper_prob, lon1, lon2, lat1, lat2     
        float [:, :] verts

    PROPORTION_TO_PERCENT = 100
    for i in range(1, prob_bins.shape[0]): # discard first bin
        lower_prob = prob_bins[i, 0]
        upper_prob = prob_bins[i, 1]
        log_arr = np.logical_and(prob_arr >= lower_prob,
                                 prob_arr <  upper_prob)
        lat_ixs, lon_ixs = np.nonzero(log_arr)  # Find the indexes of non zero elements
        for j in range(lon_ixs.shape[0]):       # Len of either lat_ixs or lon_ixs
            lon_idx = lon_ixs[j]                # find the ith index
            lat_idx = lat_ixs[j]
            lon1 = lon[lon_idx]
            lon2 = lon[lon_idx+1]
            lat1 = lat[lat_idx]
            lat2 = lat[lat_idx+1]
            # lon1, lon2 = lon[lon_idx: lon_idx+2]  # Find the actual corner values
            # lat1, lat2 = lat[lat_idx: lat_idx+2]
            verts = ((lon1, lat1),
                    (lon2, lat1),
                    (lon2, lat2),
                    (lon1, lat2),
                    (lon1, lat1))
            record = {'geometry': {"type": "Polygon",
                                   "coordinates": ((verts,))},
                      'properties': {
                          nc_var_str: lower_prob * PROPORTION_TO_PERCENT
                       }
                     }
            sink.write(record)

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

cdef float [:, :] prob_threshold_array(float step=0.02):
    """Return a array of tuples

    eg tuple (0.02, 0.04)
    """
    cdef:
        float [:] prob_vector
        float [:, :] prob_thresholds

    prob_vector = np.arange(0., 1 + 1.5 * step, step) # Including end points
    prob_thresholds = np.zeros((len(prob_vector)-1, 2))
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
        mk_shp_v4(season, nc_var, CONFIG)

    elapsed = time.time() - start_time
    print("Full run in {} seconds".format(elapsed))



    print("Done __main__")

