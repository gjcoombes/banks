# -*- coding: utf-8 -*-
"""
Module: grids.py
Created on Mon Jan 27 17:21:40 2014
@author: gav
Description:
Main module for manipulation and parsing of ASA dep and hab grids
"""
### Imports
from __future__ import print_function

import os.past as osp
import numpy as np

### Logging
import logging
logging.basicConfig(level=logging.DEBUG)
debug, info, error = logging.debug, logging.info, logging.error
### Constants

### Classes

### Functions
def gridder_arr_factory(grid_fp=None,
                        grid_ext=None,
                        grid_spc=None):
    """Return a function to convert lon, lat to row, col

    Args of factory:
        grid_fp - String of the full path to a APASA .dep grid file
        grid_extent - a dictionary of the upper right and lower left corners
            see stoch_lib.grid_extent
        grid_spacing_arr - array of lon_delta and lat_delta
            see stoch.lib.grid_spacing
    """
    err_msg = "Incomplete args. see \n {}".format(gridder_arr_factory.__doc__)
    assert grid_fp or (grid_ext and grid_spc), err_msg

    if grid_fp is not None:
        assert osp.isfile(grid_fp), "{} is not a file".format(grid_fp)
        if grid_ext is None:
            grid_ext = grid_extent(grid_fp)
        if grid_spc is None:
            grid_spc = grid_spacing_arr(grid_fp)

    origin = np.array(grid_ext['upper_left'])
    print("origin is {}".format(origin))
    delta = np.array(grid_spc).view('<f4')
    print("delta is {}".format(delta))

    def _inner(particle_pos_arr):
        # A little inconvenient to have to cast back to ndarray for ops
        new_shape = (len(particle_pos_arr), 2)
        particle_pos_nd = particle_pos_arr.view('<f4').reshape(new_shape)
        res_arr = np.floor_divide((particle_pos_nd - origin), delta)
        # Cast to int ready for geocoding
        return res_arr.astype(np.uint16)
    return _inner

def grid_extent(fp):
    """Return the upper left and lower right corners as lon, lat pairs

    Maybe return a dictionary instead?
    """
    h = grid_header(fp)
    debug("grid header is")
    debug(h)
    # The upper left corner is closest to the 0, 0 the geographical origin
    upper_left_lon = h['lon_lower_left']
    upper_left_lat = h['lat_lower_left'] + h['n_rows'] * h['lat_delta']
    lower_right_lon = h['lon_lower_left'] + h['n_cols'] * h['lon_delta']
    lower_right_lat = h['lat_lower_left']
    return {'upper_left' : [float(x) for x in [upper_left_lon,  upper_left_lat ]],
            'lower_right': [float(x) for x in [lower_right_lon, lower_right_lat]]}

def grid_header(fp):
    """Grab the header from the grid file, return as numpy record array"""
    dt_names = ['h1', 'h2', 'lon_lower_left', 'lat_lower_left',
                'lon_delta', 'lat_delta', 'n_cols', 'n_rows']
    dt_formats = ['<i2'] * 2 + ['<f4'] * 4 + ['<i2'] * 2
    dt = np.dtype(zip(dt_names, dt_formats))
    with open(fp, 'rb') as fh:
        header = np.fromfile(fh, dtype=dt, count=1)
    return header

def grid_spacing_arr(fp):
    """Return the grid spacing as dictionary """
    h = grid_header(fp)
    res_arr = h[['lon_delta', 'lat_delta']].copy()
    # Change sign of lat_delta for decreasing (-ve) lat from origin
    res_arr.dtype.names = ['lon', 'lat']
    # TODO - Make it work for the northern hemisphere
    if SOUTHERN_HEMISPHERE:
        res_arr['lat'] = -np.abs(res_arr['lat'])
#    res_arr['lat'] = res_arr['lat'] * -1
    return  res_arr


### Tests

if __name__ == "__main__":

    print("Done __main__")
