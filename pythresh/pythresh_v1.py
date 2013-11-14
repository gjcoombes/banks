# -*- coding: utf-8 -*-
"""
Module: build_particle_array.py
Created on Sat Oct 19 15:34:52 2013
@author: gav
Description:

"""
### Imports
from __future__ import print_function

import os, sys
import time

import os.path as osp
import numpy as np
import numpy.ma as ma
import pickle

from numpy import array, newaxis
import scipy.sparse as sp
import matplotlib.pyplot as plt


from pyproj import Proj
from shapely.geometry import Polygon
#
#import bonaparte
#import bonaparte.stoch.stoch_lib as sl
#import bonaparte.utils.grid_cell_areas as ga

### Logging
import logging
logging.basicConfig(level=logging.DEBUG)
debug, info, warn, error = logging.debug, logging.info, logging.warn, logging.error
### Constants

SOUTHERN_HEMISPHERE = True

SURFACE_TYPE = 0
AROMATIC_TYPE = -1
ENTRAINED_TYPE = -2
SURFACTANT_TYPE = -3

TR3_RECORD_SIZE = 40
EMPTY_HEADER_LINES = 7

### Classes

### Functions

### Grid cell area module
def ga_grid_cell_areas(grid_header=None, grid_fp=None):
    """
    Return an array of the areas of a column of grid cells
    """
    assert bool(grid_header) != osp.isfile(grid_fp or "")
    _gh = grid_header or sl_grid_header(grid_fp)
    points_array = ga_grid_column_verts(_gh)
    utm_points_array = ga_reproject(points_array)
    area_array = ga_polygon_areas(utm_points_array)
    return area_array

def ga_grid_column_verts(grid_header=None, grid_fp=None):
    """
    Given a grid, return an array of the corners of a column of the grid
    Which column is unimportant as all will have the same set of areas.

    The array will have dimensions
        n_rows,  that is the length of the column
        n_points, that is 5, the number to specify a rectangle
        n_geo_dims, that is 2 lat,lon

    eg
    [[(ax1, ay1), (ax2, ay2), (ax3, ay3), (ax4, ay4), (ax5, ay5)],
     [(bx1, by1),..                                   (bx5, by5)],]

     Fastest way to fill numpy array
     http://stackoverflow.com/questions/5891410/numpy-array-initialization-fill-with-identical-values?lq=1
    """
    # Ensure  grid_header xor grid_fp exist
    assert bool(grid_header) != osp.isfile(grid_fp or "")
    _gh = grid_header or sl_grid_header(grid_fp)
    n_rows = _gh['n_rows']
    dy     = float(_gh['lon_delta'])
    lon_0  = float(_gh['lon_lower_left'])
    lat_0  = float(_gh['lat_lower_left'])

    # Make me a function to generate vertices
    vertices = ga_verts_factory(_gh)
    # Need a sequence of lower left points
    verts_array = np.empty(shape=(n_rows, 5, 2))
    ll_corners = np.empty(shape=(n_rows, 2)) # lower left corners (lon, lat)
    ll_corners[:,0] = lon_0
    ll_corners[:,1] = np.linspace(lat_0, lat_0 + n_rows * dy, n_rows)
    verts_array[:] = np.array(map(vertices, ll_corners))
    return verts_array

def ga_polygon_areas(arr):
    """
    Given an column of points, return a column of polygon areas
    """
    ps = map(Polygon, arr)
    areas = [p.area for p in ps]
    return np.array(areas)

def ga_reproject(arr, zone=None):
    """Given an aray of points, return the utm coordinates"""
    new_arr = np.empty_like(arr)
    _zone = zone or ga_utm_zone(None)
    proj = Proj(proj="utm", zone=_zone, ellps="WGS84")
    for i, grid_cell in enumerate(arr):
        for j, point in enumerate(grid_cell):
            new_arr[i, j, :] = np.array(proj(*point))
    return new_arr

def ga_utm_zone(point):
    """
    *** Warning stub only - fixed output
    Given a geographical point, return the appropriate utm zone

    Args:
        point - array of shape (1, 2) ie (lon, lat)

    Returns:
        zone - string of the form "50L"

    """
    warn("***Warning stub function - fixed return value***")
    return "50L"

def ga_verts_factory(grid_header):
    """
    Return a function that will calculate the five verts given the lower left corner
    """
    dx = np.array([float(grid_header['lon_delta']), 0])
    dy = np.array([0, float(grid_header['lat_delta'])])

    def verts(point):
        _verts = np.empty((5,2))
        _verts[0] = point
        _verts[1] = point + dx
        _verts[2] = point + dx + dy
        _verts[3] = point + dy
        _verts[4] = point
        return _verts

    return verts

###
### Stoch_lib module ###


def sl_gridder_arr_factory(grid_fp=None,
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
    err_msg = "Incomplete args. see \n {}".format(sl_gridder_arr_factory.__doc__)
    assert grid_fp or (grid_ext and grid_spc), err_msg

    if grid_fp is not None:
        assert osp.isfile(grid_fp), "{} is not a file".format(grid_fp)
        if grid_ext is None:
            grid_ext = sl_grid_extent(grid_fp)
        if grid_spc is None:
            grid_spc = sl_grid_spacing_arr(grid_fp)

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

def sl_grid_extent(fp):
    """Return the upper left and lower right corners as lon, lat pairs

    Maybe return a dictionary instead?
    """
    h = sl_grid_header(fp)
    debug("grid header is")
    debug(h)
    # The upper left corner is closest to the 0, 0 the geographical origin
    upper_left_lon = h['lon_lower_left']
    upper_left_lat = h['lat_lower_left'] + h['n_rows'] * h['lat_delta']
    lower_right_lon = h['lon_lower_left'] + h['n_cols'] * h['lon_delta']
    lower_right_lat = h['lat_lower_left']
    return {'upper_left' : [float(x) for x in [upper_left_lon,  upper_left_lat ]],
            'lower_right': [float(x) for x in [lower_right_lon, lower_right_lat]]}

def sl_grid_header(fp):
    """Grab the header from the grid file, return as numpy record array"""
    dt_names = ['h1', 'h2', 'lon_lower_left', 'lat_lower_left',
                'lon_delta', 'lat_delta', 'n_cols', 'n_rows']
    dt_formats = ['<i2'] * 2 + ['<f4'] * 4 + ['<i2'] * 2
    dt = np.dtype(zip(dt_names, dt_formats))
    with open(fp, 'rb') as fh:
        header = np.fromfile(fh, dtype=dt, count=1)
    return header

def sl_grid_mass_dense(particles, gridder, grid_shape):
    """Grid the mass of the particles and return a dense array

    Args:
        particles - array from the chunker
        gridder - function to return the grid indices
        grid_shape - tuple of (n_rows, n_cols)
    """
    idxs = gridder(particles[['lon', 'lat']])
    _data = particles['mass']
    _col = idxs[:,0]     # lon
    _row = idxs[:,1]     # lat
    mass_coo = sp.coo_matrix((_data, (_row, _col)), shape=grid_shape)
    return mass_coo.todense()

def sl_grid_mass_sparse(particles, grid_extent, grid_spacing):
    """
    Take a full sets of particles, grid and return the particles
    origin is [ 105.  -15.]
    delta is [ 0.00764526 -0.00722456]
    """
    def log(log_msg): debug("sl_grid_mass_sparse: {}".format(log_msg))

    origin = np.array(grid_extent['upper_left'])
    print("origin is {}".format(origin))
    delta = np.array(grid_spacing).view('<f4')
    ijv = particles[['lon', 'lat', 'mass']].copy()
    log("ijv flags are %s" % ijv.flags)
    log("ijv dtype is %s" % ijv.dtype)
    log("ijv shape is %s" % ijv.shape)

    n_rows = ijv.shape[0]
    result = np.empty((n_rows, 3))
    for n in range(n_rows):
        result[n, 0] = (ijv[n, 0] - origin[0]) // delta[0]
        result[n, 1] = (ijv[n, 1] - origin[1]) // delta[1]
        result[n, 2] =  ijv[n, 2]








def sl_grid_mass_csr(particles, gridder, grid_shape):
    """Grid the mass of the particles and return a dense array

    Args:
        particles - array from the chunker
        gridder - function to return the grid indices
        grid_shape - tuple of (n_rows, n_cols)
    """
    idxs = gridder(particles[['lon', 'lat']])
    _data = particles['mass']
    _col = idxs[:,0]     # lon
    _row = idxs[:,1]     # lat
    mass_coo = sp.coo_matrix((_data, (_row, _col)), shape=grid_shape)
    return mass_coo.tocsr()

def sl_grid_spacing_arr(fp):
    """Return the grid spacing as dictionary """
    h = sl_grid_header(fp)
    res_arr = h[['lon_delta', 'lat_delta']].copy()
    # Change sign of lat_delta for decreasing (-ve) lat from origin
    res_arr.dtype.names = ['lon', 'lat']
    # TODO - Make it work for the northern hemisphere
    # Should be easiest to adopt a lowerleft orgin for both hemispheres
    if SOUTHERN_HEMISPHERE:
        res_arr['lat'] = -np.abs(res_arr['lat'])
#    res_arr['lat'] = res_arr['lat'] * -1
    return  res_arr


def sl_lu3_data(lu3_fp):
    """Return binary data from file"""
    dtype_ls = [('iver', '<i4'), ('SIMTIME', '<i4'),
                ('time', '<f4'), ('rec1st', '<i4'), ('rec2st', '<i4'),
                ('rec2end', '<i4'), ('rec3st', '<i4'), ('rec3end', '<i4'),
                ('rec5st', '<i4'), ('rec5end', '<i4'), ('sed2st', '<i4'),
                ('sed2end', '<i4'), ('rar1st', '<i4'), ('rar1end', '<i4'),
                ('rth1st', '<i4'), ('rth1end', '<i4'), ('rsf1st', '<i4'),
                ('rsf1end', '<i4'), ('rsp1st', '<i4'), ('rsp1end', '<i4'),
                ('rss1st', '<i4'), ('rss1end', '<i4'), ('rat1st', '<i4'),
                ('rat1end', '<i4')]
    return np.fromfile(lu3_fp, np.dtype(dtype_ls))

def sl_particles_and_shore_generator(tr3_fp, lu3_fp, grid_fp):
    """Return a generator of particles and shore cells
    Note mismatch in dtype lengths between particles and shore
    Usage:
    > gen = particles_and_shore_generator(tr3_fp, lu3_fp, grid_fp)
    > for time, surf, entr, arom , shore in gen:
    >     ...

    Yields (time, surf_p, entr_p, arom_p, shore_c)
    """
    def log(log_msg): debug("particle_and_shore_generator: {}".format(log_msg))
    if __debug__:
        log("TR3 file is {}".format(tr3_fp))
        log("LU3 file is {}".format(lu3_fp))
        log("Grid file is {}".format(grid_fp))


    grid_record = sl_grid_extent(grid_fp)
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

    def inner():
        with open(tr3_fp, 'rb') as fh:
            for row in lu3_arr:
                empty = np.fromfile(fh, dtype=particle_dtype, count=EMPTY_HEADER_LINES)
                particles = np.fromfile(fh, dtype=particle_dtype,
                                        count=row['rec2end']-row['rec2st'])
                shore_cells = np.fromfile(fh, dtype=shore_dtype,
                                        count=row['rec3end']-row['rec3st'] or 1)
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

    return inner()

def main():
    """
    Main func
    """
    project_dir = r"J:\data\j0267_nv_remodel"
    stem = "J0267_SC3_SBED_LEAK_TRA_001"
    grid_fn = "VanGogh_800m.DEP"
    h5_fn = "j0267_data.h5"

    tr3_fp = osp.join(project_dir, "modelout", stem + ".tr3" )
    lu3_fp = osp.join(project_dir, "modelout", stem + ".lu3" )
    grid_fp = osp.join(project_dir, "grids", grid_fn)

    surf_threshold = 1e-6 # T/m2 or 1 g/m2

    header = sl_grid_header(grid_fp)
    grid_shape = (header['n_rows'], header['n_cols'])

    particles = sl_particles_and_shore_generator(tr3_fp, lu3_fp, grid_fp)
    gridder = sl_gridder_arr_factory(grid_fp)

    start_time = time.time()
    max_surf_mass = np.zeros(grid_shape, dtype=np.float32)

    for i, tup in enumerate(particles):
        sim_time, surf, entr, arom, shore = tup
        surf_dense = sl_grid_mass_dense(surf, gridder, grid_shape)
        max_surf_mass = np.maximum(max_surf_mass, surf_dense)

    # Now we need a threshold_matrix to find which cells have exceeded the threshold
    cell_areas = ga_grid_cell_areas(grid_fp=grid_fp)
    # the mass threshold is the threshold * area eg 0.001 kg/m2 * 640000 m2 = mass in Ts
    mass_threshold_T = cell_areas * surf_threshold
    exceedance = max_surf_mass >= mass_threshold_T[:, np.newaxis]
    max_mass = ma.array(max_surf_mass, mask=(max_surf_mass == 0.0))
    elapsed_time = time.time() - start_time
    print("Finished {} timesteps in {} seconds".format(i, elapsed_time))
#    plt.imshow(exceedance, origin="upper", interpolation="nearest")
#    plt.show()

### Tests

if __name__ == "__main__":

    main()


    print("Done __main__")