# -*- coding: utf-8 -*-
# cython: profile=True
"""
Module: cythresh_v0_2.py
Created on Thu Nov 14 21:34:28 2013
@author: gav
Description: Incremental improvement of pythresh

test_particles runs in 7.22 seconds per loop
adding contiguous (but -OO) raises to 10.2 s
runnng non-optimised but with silent logs adds 5s -> 15s

running main first 400 chunnks in < 15 s
Chief bottlenecks at grid_mass_sparse, update arrays.
Particularly divmod

%lprun on grid_mass_sparse confirm >70% of time in function spent on divmod
%lprun on update_agg_arrays shows even time on array lookups,
    need to check if separate vecs or single array best for memory locality

Typing variables in sl_grid_mass_sparse took func time from 7.6s to 0.7s
Now timeit reports best loop of 7.002s for first 400 chunks

Typing variables in update_agg took func time from 5.7 to 0s?
With profiling off timeit has 1.37s best loop for 400 chunks

Final optimisations make full run in 13s
or with ten thresholds in 16 s (9s on hyde)

Now with timeseries
"""
### Imports
from __future__ import print_function

import os.path as osp
import matplotlib.pyplot as plt

import numpy as np
cimport numpy as np

import cython
cimport cython

from libc.math cimport modf
from cpython.array cimport array, clone

import scipy.sparse as sp
from time import time
import tables as tb

from pyproj import Proj
from shapely.geometry import Polygon

### Logging
import logging
logger = logging.getLogger("pythresh")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)

debug, info, warn, error = logger.debug, logger.info, logger.warn, logger.error
### Constants
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

SURFACE_TYPE = 0
AROMATIC_TYPE = -1
ENTRAINED_TYPE = -2
SURFACTANT_TYPE = -3

TR3_RECORD_SIZE = 40
EMPTY_HEADER_LINES = 7
### Classes

### Functions


### Stoch_lib module functions here ###
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
    warn("***Warning ga_utm_zone:stub function - fixed return value***")
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

def sl_contiguous_doubles(p_rec):
    """
    Return a copied C contiguous array of lon, lats, mass (doubles) <f8 or np.float64
    """
    arr = np.empty((p_rec.shape[0], 3), dtype=np.float64)
    arr[:,0] = p_rec['lon']
    arr[:,1] = p_rec['lat']
    arr[:,2] = p_rec['mass']
    return arr

def sl_contiguous_singles(p_rec):
    """
    Return a copied C contiguous array of lon, lats, mass (single floats) <f4 or np.float32
    """
    arr = p_rec[['lon', 'lat', 'mass']]
    return arr

def sl_grid_extent(h):
    """Return the upper left and lower right corners as lon, lat pairs

    Maybe return a dictionary instead?
    """
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

#    if __debug__:
#        log("header is {}".format(header))
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

@cython.boundscheck(False)
@cython.wraparound(False)
cdef sl_grid_mass_sparse(double [:,::1] particles,
                        double olon,
                        double olat,
                        double dlon,
                        double dlat,
                        Py_ssize_t ni,
                        Py_ssize_t nj):
    """
    Grid the particles into i, j, mass triples

    Args:
        particles: a c-contiguous array of doubles lon, lat, mass
        olon: a double representing the longitude of the origin
        olat: a double representing the latitude of the origin
        dlon: a double representing the column spacing
        dlat: a double representing the row spacing
        ni: Do i need these?
        nj:

    Returns:
        arr_ij : int array of i, j grid indices
        arr_mass: double array of masses
            Note: the output array must have the same length
    """
    cdef Py_ssize_t N = particles.shape[0]
    cdef Py_ssize_t m
    cdef double i_dbl
    cdef double j_dbl

    np_arr   = np.empty((N,), dtype=np.int32)
    np_arr2 = np.empty_like(np_arr)

    cdef int [::1] particle_i_arr = np_arr
    cdef int [::1] particle_j_arr = np_arr2
    cdef double [:] particle_mass_arr = particles[:,2]

    for m in xrange(N):
        modf( (particles[m,0] - olon) / dlon, &i_dbl)
        modf( (particles[m,1] - olat) / dlat, &j_dbl)
        particle_i_arr[m] = <int>i_dbl
        particle_j_arr[m] = <int>j_dbl

    arr_i, arr_j, arr_mass = sl_sum_particles(particle_i_arr,
                                              particle_j_arr,
                                              particle_mass_arr,
                                              ni, nj)
    return arr_i, arr_j, arr_mass


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
    data = np.fromfile(lu3_fp, np.dtype(dtype_ls))
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

    grid_header = sl_grid_header(grid_fp)
    grid_record = sl_grid_extent(grid_header)
    lower_lon, upper_lat = grid_record['upper_left']
    upper_lon, lower_lat = grid_record['lower_right']

    lu3_arr = sl_lu3_data(lu3_fp)
    particle_names_ls = ['lon', 'lat', 'radius', 'prev_lon', 'prev_lat',
                         'type',
                         'mass', 'density', 'viscosity', 'age']
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

def sl_sum_particles(particle_i_arr,
                     particle_j_arr,
                     particle_mass_arr,
                     ni, nj):
    """
    Sum the masses of particles with the same grid indices

    The implementation is currently using the scipy sparse coo matrix,
    but could be done with a dictionary with a (i, j) tuple as key and mass as value

    Args:
        particle_ij_arr: array of ints with the i, j indices of each particle
        particle_mass_arr: array of doubles with the mass of each particle

    Returns:
        arr_ij: array of ints with the i'j indices of each summed grid cell
        arr_mass: array of doubles with mass of each grid cell
    """
    _data = particle_mass_arr
    _col = particle_i_arr # lon
    _row = particle_j_arr # lat
    mass_coo = sp.coo_matrix((_data, (_row, _col)), shape=(nj, ni))
    return mass_coo.col, mass_coo.row, mass_coo.data

@cython.boundscheck(False)
@cython.wraparound(False)
cdef sl_update_agg_arrays(int sim_time,
                         int [::1] i_vec,
                         int [::1] j_vec,
                         double [:] mass_vec,
                         double [:,:] max_mass_arr,
                         int    [:,:,:] min_time_arr,
                         double [:,:] mass_threshold_cols
                         ):
    """
    Update the aggregate arrays in-place

    Find where the new cells are higher than the aggregate mass
    Update those cells in the mass array
    Update those cells with the time in the time array
    """
    cdef Py_ssize_t n
    cdef Py_ssize_t i
    cdef Py_ssize_t j
    cdef Py_ssize_t k
    cdef double mass
    cdef double threshold
    cdef Py_ssize_t n_cells = i_vec.shape[0]
    cdef Py_ssize_t n_thresholds = mass_threshold_cols.shape[1]

    for n in xrange(n_cells):
        i = i_vec[n]
        j = j_vec[n]
        mass = mass_vec[n]
        if mass > max_mass_arr[j, i]:
            max_mass_arr[j, i] = mass

        for k in range(n_thresholds):
            try:
                if mass > mass_threshold_cols[j,k] and sim_time < min_time_arr[k,j,i]:
                    min_time_arr[k,j,i] = sim_time
            except:
                print("simtime: {}, k: {}, j: {}, i: {}".format(sim_time, k, j, i))

                raise
    return max_mass_arr, min_time_arr

def main(tr3_fp, lu3_fp, grid_fp, h5_fp, surf_thresholds_gm2):
    """
    """
    def log(log_msg): debug("main: {}".format(log_msg))

    start_time = time()

    cell_areas = ga_grid_cell_areas(grid_fp=grid_fp)
    mass_threshold_cols = cell_areas[:,None] * surf_thresholds_gm2[None,:]

    grid_header = sl_grid_header(grid_fp)
    olon, olat, dlon, dlat, ni, nj = parse_grid_header(grid_header)

    max_mass_arr = np.zeros((nj, ni))
    min_time_arr = np.empty((surf_thresholds_gm2.shape[0],nj, ni), dtype=np.int32)
    min_time_arr.fill(1e9)

    with tb.open_file(h5_fp, "w") as root:
        tb_int = tb.Int32Atom
        result_grp = root.create_group("/", "result", "pythresh results")
        i_vlarr = root.create_vlarray("/result", "i_vec", tb.Int32Atom())
        j_vlarr = root.create_vlarray("/result", "j_vec", tb.Int32Atom())
        mass_vlarr = root.create_vlarray("/result", "mass_vec", tb.Float64Atom())

        for i, tup in enumerate(sl_particles(tr3_fp, lu3_fp, grid_fp)):
            sim_time, surf, entr, arom, shore = tup
            surf_arr = sl_contiguous_doubles(surf)
            if len(surf_arr) == 0:
                continue
            i_vec, j_vec, mass_vec = sl_grid_mass_sparse(surf_arr, olon, olat, dlon, dlat, ni, nj)
            i_vlarr.append(i_vec)
            j_vlarr.append(j_vec)
            mass_vlarr.append(mass_vec)
            # Write the timestep array to hdf5 here
            max_mass_arr, min_time_arr = sl_update_agg_arrays(sim_time,
                                                              i_vec, j_vec, mass_vec,
                                                              max_mass_arr, min_time_arr, mass_threshold_cols)

        exceedance_arr = np.empty((surf_thresholds_gm2.shape[0], nj, ni), dtype=np.bool)
        for k, threshold in enumerate(surf_thresholds_gm2):
            threshold_cols = mass_threshold_cols[:,k]

            ex = max_mass_arr >= threshold_cols[:, np.newaxis]
            exceedance_arr[k, :, :] = ex
        # Write agg_arrays to hdf5 here
        root.create_carray(result_grp, "max_mass", obj=np.array(max_mass_arr))
        root.create_carray(result_grp, "min_time", obj=np.array(min_time_arr))
        root.create_carray(result_grp, "exceedance", obj=np.array(exceedance_arr))
    elapsed_time = time() - start_time
    info("test in {} seconds".format(elapsed_time))
    result = {
        "max_mass_arr"  : max_mass_arr,
        "min_time_arr"  : min_time_arr,
        "exceedance_arr": exceedance_arr,
    }
    return result

### Tests
def particles_fixture():
    project_dir = r"J:\data\j0267_nv_remodel"
    stem = "J0267_SC3_SBED_LEAK_TRA_001"
    grid_fn = "VanGogh_800m.DEP"

    tr3_fp = osp.join(project_dir, "modelout", stem + ".tr3" )
    lu3_fp = osp.join(project_dir, "modelout", stem + ".lu3" )
    grid_fp = osp.join(project_dir, "grids", grid_fn)
    return sl_particles(tr3_fp, lu3_fp, grid_fp)

def test_contiguous_doubles():
    particles = particles_fixture()
    for i, tup in enumerate(particles):
        sim_time, surf, entr, arom, shore = tup
        c_arr = sl_contiguous_doubles(surf)
        if i > 5:
            break

def test_contiguous_singles():
    particles = particles_fixture()
    for i, tup in enumerate(particles):
        sim_time, surf, entr, arom, shore = tup
        c_arr = sl_contiguous_singles(surf)
        if i > 3:
            break

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

def test_math():
    cdef double n = 1.25
    cdef double d = 0.1
    cdef double f
    cdef int i
    modf(n / d, &f)
    i = <int>f
    print(i)

def test_array():
    cdef int i
    cdef array template = array('d', [0])
    tpl_arr = clone(template, 12, False)
    print(tpl_arr)

if __name__ == "__main__":

#    test_contiguous_doubles()
#    test_contiguous_singles()

#    main()
    test_array()




    print("Done __main__")
