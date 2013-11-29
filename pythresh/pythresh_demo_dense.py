# -*- coding: utf-8 -*-
"""
Module: pythresh_demo_dense.py
Created on Sat Oct 19 15:34:52 2013
@author: gav
Description:
This is a an example of how stochastic threshold processing could be
implemented in python/numpy.

The chief difference between this implementation and the current matlab one
is (I believe) the individual processing of timesteps

Prelim test runs:
=================
For the Hyde Scenario 1 Summer 001 run:
Pythresh: Finished 2520 timesteps in 174.411999941 seconds
Rough comparison to StochIndiTest
pythresh:  Finished in 175 seconds
           mem: 0.17 GB cpu: 13% on 1 core
StochIndi: Finished in 265 seconds
           mem: 5.2 GB cpu: ~35% across 4 cores

For the Northern Endeavor Scenario 2a run:
pythresh:  Finished 2688 timesteps in 156 seconds
           0.15 GB RAM 13% CPU (one processor)
StochIndi: Finished all timesteps in ~240 seconds
           3.5 GB; ~30% CPU (peak 50% it think)
"""
### Imports
from __future__ import print_function

import os, sys
import time

import os.path as osp
import numpy as np
import numpy.ma as ma
import pickle
import tables as tb

from numpy import array, newaxis
from numpy import logical_and as np_and

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
logger = logging.getLogger("pythresh")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)
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
def parse_grid_header(grid_header):
    """
    Separate single parameters and cast to correct type for gridding

    Returns:
        olon, olat, dlon, dlat
    """
    gh = grid_header
    olon = float(gh['lon_lower_left'])
    olat = float(gh['lat_lower_left'])
    dlon = float(gh['lon_delta'])
    dlat = float(gh['lat_delta'])
    ni = int(gh['n_cols'])
    nj = int(gh['n_rows'])
    return (olon, olat, dlon, dlat, ni, nj)

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
#    print("origin is {}".format(origin))
    delta = np.array(grid_spc).view('<f4')
#    print("delta is {}".format(delta))

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
#    if len(mass_coo.data) > 0:
#        print("sl_grid_mass_dense max mass is {}".format(np.max(mass_coo.data)))
    result = mass_coo.todense()
#    print("dense arr max mass is {}".format(np.max(result)))
    return result

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

def sl_update_agg_arrays(sim_time, surf_dense,
                         max_mass_arr, min_time_arr, mass_threshold_col):
    """
    Update the aggregate arrays in-place

    Find where the new cells are higher than the aggregate mass
    Update those cells in the mass array
    Update those cells with the time in the time array
    """
    # Max mass for all time steps
    max_mass_arr = np.maximum(max_mass_arr, surf_dense)
    # Time of earliest exceedance
    exceedance_mask = surf_dense >= mass_threshold_col[:, np.newaxis]
    earliest_mask = sim_time < min_time_arr
    update_time_mask = np_and(exceedance_mask, earliest_mask)
    min_time_arr[update_time_mask] = sim_time
    return max_mass_arr, min_time_arr


def main(tr3_fp, lu3_fp, grid_fp, h5_fp):
    """
    Main func
    """

    tr3_fp = osp.join(project_dir, "modelout", stem + ".tr3" )
    lu3_fp = osp.join(project_dir, "modelout", stem + ".lu3" )
    grid_fp = osp.join(project_dir, "grids", grid_fn)

    header = sl_grid_header(grid_fp)
    grid_shape = (header['n_rows'], header['n_cols'])

    particles = sl_particles_and_shore_generator(tr3_fp, lu3_fp, grid_fp)
    gridder = sl_gridder_arr_factory(grid_fp)

    start_time = time.time()

    # Now we need a threshold_matrix to find which cells have exceeded the threshold
    surf_threshold = 1e-6 # T/m2 or 1 g/m2
    cell_areas = ga_grid_cell_areas(grid_fp=grid_fp)
    mass_threshold_col = cell_areas * surf_threshold
    # Build the result arrays
    grid_header = sl_grid_header(grid_fp)
    olon, olat, dlon, dlat, ni, nj = parse_grid_header(grid_header)
    max_mass_arr = np.zeros((nj, ni))
    min_time_arr = np.empty((nj, ni))
    min_time_arr.fill(1e9)

    for i, tup in enumerate(particles):
#        if i > 200:
#            warn("*** Debug: Breaking processing at step {}".format(i))
#            break
        sim_time, surf, entr, arom, shore = tup
        surf_dense = sl_grid_mass_dense(surf, gridder, grid_shape)
#        print("main surf_dense max is {}".format(np.max(surf_dense)))
        max_mass_arr, min_time_arr = sl_update_agg_arrays(sim_time, surf_dense,
                             max_mass_arr, min_time_arr, mass_threshold_col)

    exceedance_arr = max_mass_arr >= mass_threshold_col[:, np.newaxis]

    elapsed_time = time.time() - start_time
    print("Finished {} timesteps in {} seconds".format(i, elapsed_time))
    result = {
        "max_mass_arr"  : max_mass_arr,
        "min_time_arr"  : min_time_arr,
        "exceedance_arr": exceedance_arr,
    }
    return result

def store(result, fp, method="pickle"):
    """Store the results in a file and return the filepath"""
    if method == "pickle":
        with open(fp, "wb") as sink:
            pickle.dump(result, sink)
    elif method == "hdf5":
        with tb.open_file(fp, "w") as root:
            result_grp = root.create_group("/", "result", "pythresh results")
            root.create_carray(result_grp, "max_mass", obj=result['max_mass_arr'])
            root.create_carray(result_grp, "min_time", obj=result['min_time_arr'])
            root.create_carray(result_grp, "exceedance", obj=result['exceedance_arr'])


    return fp


### Tests

if __name__ == "__main__":

    project_dir = r"E:\Loc_Data\pythresh_trials"
    stem = "J0272_SC1_SURF_MUTI_255M3HR_SUM_001"
    grid_fn = "MidWA_NWS_1km.DEP"
#    stem = "J0266_SC2A_SUBS_18962M3_LAM_Q1_001"
#    grid_fn = "Penguin_1000m.DEP"

    h5_fn = "j0272_data.h5"
    pkl_fn = stem + ".pkl"

    tr3_fp = osp.join(project_dir, "modelout", stem + ".tr3" )
    lu3_fp = osp.join(project_dir, "modelout", stem + ".lu3" )
    grid_fp = osp.join(project_dir, "grids", grid_fn)
    h5_fp = osp.join(project_dir, "hdf5", h5_fn)
    pkl_fp = osp.join(project_dir, "hdf5", pkl_fn)

#    result = main(tr3_fp, lu3_fp, grid_fp, h5_fp)
#    store(result, pkl_fp, "pickle")

#    with open(pkl_fp, "rb") as source:
#        result = pickle.load(source)
#    store(result, h5_fp, "hdf5")

    with tb.open_file(h5_fp, "r") as h5_src:
        print(h5_src)
        max_mass = h5_src.get_node("/result", "max_mass").read()
        min_time = h5_src.get_node("/result", "min_time").read()
        exceedance = h5_src.get_node("/result", "exceedance").read()

#    exceedance = result['exceedance_arr']
#    max_mass = result['max_mass_arr']
#    min_time = result['min_time_arr']
    min_time[min_time == 1e9] = np.nan

    plt.imshow(exceedance, origin="upper", interpolation="nearest")
    plt.colorbar()
    plt.show()









    print("Done __main__")
