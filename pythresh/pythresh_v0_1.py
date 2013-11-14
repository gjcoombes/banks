# -*- coding: utf-8 -*-
"""
Module: pythresh_v0_1.py
Created on Thu Nov 14 21:34:28 2013
@author: gav
Description: Incremental imporovement of pythresh

test_particles runs in 7.22 seconds per loop

"""
### Imports
from __future__ import print_function

import os.path as osp
import numpy as np
from time import time

### Logging
import logging
logging.basicConfig(level=logging.DEBUG)
debug, info, error = logging.debug, logging.info, logging.error
### Constants
SURFACE_TYPE = 0
AROMATIC_TYPE = -1
ENTRAINED_TYPE = -2
SURFACTANT_TYPE = -3

TR3_RECORD_SIZE = 40
EMPTY_HEADER_LINES = 7
### Classes

### Functions

### Stoch_lib module functions here ###
def sl_grid_extent(h):
    """Return the upper left and lower right corners as lon, lat pairs

    Maybe return a dictionary instead?
    """
    if __debug__:
        def log(log_msg):
            debug("sl_grid_extent: {}".format(log_msg))
        log("header is {}".format(h))

    # The upper left corner is closest to the 0, 0 the geographical origin
    upper_left_lon = h['lon_lower_left']
    upper_left_lat = h['lat_lower_left'] + h['n_rows'] * h['lat_delta']
    lower_right_lon = h['lon_lower_left'] + h['n_cols'] * h['lon_delta']
    lower_right_lat = h['lat_lower_left']
    return {'upper_left' : [float(x) for x in [upper_left_lon,  upper_left_lat ]],
            'lower_right': [float(x) for x in [lower_right_lon, lower_right_lat]]}

def sl_grid_header(fp):
    """Grab the header from the grid file, return as numpy record array"""
    if __debug__:
        def log(log_msg):
            debug("sl_grid_header: {}".format(log_msg))
        log("grid_fp is {}".format(fp))
        assert osp.isfile(fp), "No such file %s" % fp

    dt_names = ['h1', 'h2', 'lon_lower_left', 'lat_lower_left',
                'lon_delta', 'lat_delta', 'n_cols', 'n_rows']
    dt_formats = ['<i2'] * 2 + ['<f4'] * 4 + ['<i2'] * 2
    dt = np.dtype(zip(dt_names, dt_formats))
    with open(fp, 'rb') as fh:
        header = np.fromfile(fh, dtype=dt, count=1)

    if __debug__:
        log("header is {}".format(header))
    return header

def sl_lu3_data(lu3_fp):
    """Return binary data from file"""
    if __debug__:
        def log(log_msg):
            debug("sl_lu3_data: {}".format(log_msg))
        log("lu3_fp is {}".format(lu3_fp))
        assert osp.isfile(lu3_fp), "No such file %s" % lu3_fp

    dtype_ls = [('iver', '<i4'), ('SIMTIME', '<i4'),
                ('time', '<f4'), ('rec1st', '<i4'), ('rec2st', '<i4'),
                ('rec2end', '<i4'), ('rec3st', '<i4'), ('rec3end', '<i4'),
                ('rec5st', '<i4'), ('rec5end', '<i4'), ('sed2st', '<i4'),
                ('sed2end', '<i4'), ('rar1st', '<i4'), ('rar1end', '<i4'),
                ('rth1st', '<i4'), ('rth1end', '<i4'), ('rsf1st', '<i4'),
                ('rsf1end', '<i4'), ('rsp1st', '<i4'), ('rsp1end', '<i4'),
                ('rss1st', '<i4'), ('rss1end', '<i4'), ('rat1st', '<i4'),
                ('rat1end', '<i4')]
    data = np.fromfile(lu3_fp, np.dtype(dtype_ls))

    if __debug__:
        log("lu3 data starts \n {}".format(data[:6]))
        log("lu3 data ends \n {}".format(data[-5:]))
    return data


def sl_particles(tr3_fp, lu3_fp, grid_fp):
    """
    Yield surf, entr, arom, shor particles for each time step.

    Fields for surf, entr, arom records:
        'lon', 'lat', 'radius', 'prev_lon', 'prev_lat',
        'type', 'mass', 'density', 'viscosity', 'age'

    Fields for shor records:
        'igrid', 'jgrid', '_1', '_2', 'habitat_type', 'area',
        'shore_length', 'lon', 'lat', 'mass'
    """

    if __debug__:
        def log(log_msg): debug("sl_particles: {}".format(log_msg))
        log("TR3 file is {}".format(tr3_fp))
        log("LU3 file is {}".format(lu3_fp))
        log("Grid file is {}".format(grid_fp))

    grid_header = sl_grid_header(grid_fp)
    grid_record = sl_grid_extent(grid_header)
    lower_lon, upper_lat = grid_record['upper_left']
    upper_lon, lower_lat = grid_record['lower_right']

    lu3_arr = sl_lu3_data(lu3_fp)
    particle_names_ls = ['lon', 'lat', 'radius', 'prev_lon', 'prev_lat',
                         'type', 'mass', 'density', 'viscosity', 'age']
    particle_formats_ls = ['<f4'] * 5 + ['<i4'] + ['<f4'] * 4
    particle_dtype = np.dtype(zip(particle_names_ls, particle_formats_ls))
    shore_names_ls = ['igrid', 'jgrid', '_1', '_2', 'habitat_type', 'area',
                      'shore_length', 'lon', 'lat', 'mass' ]
    shore_formats_ls = ['<i4'] * 5 + ['<f4'] * 5
    shore_dtype = np.dtype(zip(shore_names_ls, shore_formats_ls))

    with open(tr3_fp, 'rb') as fh:
        for row in lu3_arr:
            empty = np.fromfile(fh, dtype=particle_dtype, count=EMPTY_HEADER_LINES)
            particles = np.fromfile(fh, dtype=particle_dtype,
                                    count=row['rec2end'] - row['rec2st'])
            shore_cells = np.fromfile(fh, dtype=shore_dtype,
                                    count=row['rec3end'] - row['rec3st'] or 1)
            np_and = np.logical_and
            surf_mask = np.array(particles['type'] == SURFACE_TYPE)
            entr_mask = np.array(particles['type'] == ENTRAINED_TYPE)
            arom_mask = np.array(particles['type'] == AROMATIC_TYPE)
            lon_mask = np_and(np.array(particles['lon'] > lower_lon),
                              np.array(particles['lon'] < upper_lon))
            lat_mask = np_and(np.array(particles['lat'] > lower_lat),
                              np.array(particles['lat'] < upper_lat))
            bounds_mask = np_and(lon_mask, lat_mask)
            surf_p = particles[np_and(bounds_mask, surf_mask)]
            entr_p = particles[np_and(bounds_mask, entr_mask)]
            arom_p = particles[np_and(bounds_mask, arom_mask)]
            yield (row['time'], surf_p, entr_p, arom_p, shore_cells)

### Tests
def test_particles():
    def log(log_msg): debug("test_particles: {}".format(log_msg))

    start_time = time()
    project_dir = r"J:\data\j0267_nv_remodel"
    stem = "J0267_SC3_SBED_LEAK_TRA_001"
    grid_fn = "VanGogh_800m.DEP"

    tr3_fp = osp.join(project_dir, "modelout", stem + ".tr3" )
    lu3_fp = osp.join(project_dir, "modelout", stem + ".lu3" )
    grid_fp = osp.join(project_dir, "grids", grid_fn)

    for i, tup in enumerate(sl_particles(tr3_fp, lu3_fp, grid_fp)):
        sim_time, surf, entr, arom, shore = tup
#        print(sim_time)
    elapsed_time = time() - start_time
    log("test in {} seconds".format(elapsed_time))

if __name__ == "__main__":

    test_particles()





    print("Done __main__")
