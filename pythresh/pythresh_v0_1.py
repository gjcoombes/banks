# -*- coding: utf-8 -*-
"""
Module: pythresh_v0_1.py
Created on Thu Nov 14 21:34:28 2013
@author: gav
Description: Incremental imporovement of pythresh

test_particles runs in 7.22 seconds per loop
adding contiguous (but -OO) raises to 10.2 s
runnng non-optimised but with silent logs adds 5s -> 15s

"""
### Imports
from __future__ import print_function

import os.path as osp
import numpy as np
importr scipy.sparse as sp
from time import time

### Logging
import logging
logger = logging.getLogger("pythresh")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)

debug, info, error = logger.debug, logger.info, logger.error
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
    if __debug__:
        def log(log_msg):
            debug("sl_contiguous_doubles: {}".format(log_msg))
        log("start of recarray is {}".format(p_rec[:6]))
        log("end of recarray is {}".format(p_rec[-5:]))
        log("shape of recarray is {}".format(p_rec.shape))
        log("dtype of recarray is {}".format(p_rec.dtype))
        log("flags of recarray is {}".format(p_rec.flags))

    arr = p_rec[['lon', 'lat', 'mass']].astype(np.float64)

    if __debug__:
        log("start of arr is {}".format(arr[:6]))
        log("end of arr is {}".format(arr[-5:]))
        log("shape of arr is {}".format(arr.shape))
        log("dtype of arr is {}".format(arr.dtype))
        log("flags of arr is {}".format(arr.flags))
    return arr

def sl_contiguous_singles(p_rec):
    """
    Return a copied C contiguous array of lon, lats, mass (single floats) <f4 or np.float32
    """
    if __debug__:
        def log(log_msg):
            debug("sl_contiguous_doubles: {}".format(log_msg))
        log("start of recarray is {}".format(p_rec[:6]))
        log("end of recarray is {}".format(p_rec[-5:]))
        log("shape of recarray is {}".format(p_rec.shape))
        log("dtype of recarray is {}".format(p_rec.dtype))
        log("flags of recarray is {}".format(p_rec.flags))

    arr = p_rec[['lon', 'lat', 'mass']]

    if __debug__:
        log("start of arr is {}".format(arr[:6]))
        log("end of arr is {}".format(arr[-5:]))
        log("shape of arr is {}".format(arr.shape))
        log("dtype of arr is {}".format(arr.dtype))
        log("flags of arr is {}".format(arr.flags))
    return arr

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

def sl_grid_mass_sparse(particles, olon, olat, dlon, dlat, ni, nj):
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
    if __debug__:
        def log(log_msg):
            info("sl_grid_mass_sparse: {}".format(log_msg))
        log("shape of particles is %s" % particles.shape)
        log("start of particles is %s" % particles[:6])
        log("end of particles is %s" % particles[-5:])

    N = particles.shape[1]
    for i in range(N):


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

def sl_update_agg_arrays(sim_time, ij_arr, mass_arr,
                         max_mass_arr, min_time_arr, mass_threshold_col):
    """
    Update the aggregate arrays in-place

    Find where the new cells are higher than the aggregate mass
    Update those cells in the mass array
    Update those cells with the time in the time array
    """
    pass

def main():
    """
    """
    def log(log_msg): debug("main: {}".format(log_msg))

    start_time = time()
    project_dir = r"J:\data\j0267_nv_remodel"
    stem = "J0267_SC3_SBED_LEAK_TRA_001"
    grid_fn = "VanGogh_800m.DEP"

    tr3_fp = osp.join(project_dir, "modelout", stem + ".tr3" )
    lu3_fp = osp.join(project_dir, "modelout", stem + ".lu3" )
    grid_fp = osp.join(project_dir, "grids", grid_fn)

    surf_threshold = 1e-6 # T/m2 or 1 g/m2
    cell_areas = ga_grid_cell_areas(grid_fp=grid_fp)
    mass_threshold_col = None

    grid_header = sl_grid_header(grid_fp)
    olon, olat, dlon, dlat, ni, nj = parse_grid_header(grid_header)
    max_mass_arr = np.zeros((nj, ni))
    min_time_arr = np.empty_like(max_mass_arr).filled(1e9)

    for i, tup in enumerate(sl_particles(tr3_fp, lu3_fp, grid_fp)):
        sim_time, surf, entr, arom, shore = tup
        surf_arr = sl_contiguous_doubles(surf)
#        # Original code
#        header = sl_grid_header(grid_fp)
#        grid_shape = (header['n_rows'], header['n_cols'])
#        gridder = sl_gridder_arr_factory(grid_fp)
#        surf_dense = sl_grid_mass_dense(surf, gridder, grid_shape)
#        max_surf_mass = np.maximum(max_surf_mass, surf_dense)

        # what we want
        ij_arr, mass_arr = sl_grid_mass_sparse(surf_arr, olon, olat, dlon, dlat)
        # Write the timestep array to hdf5 here
        sl_update_agg_arrays(sim_time, ij_arr, mass_arr,
                             max_mass_arr, min_time_arr, mass_threshold_col)

#        print(sim_time)
    elapsed_time = time() - start_time
    info("test in {} seconds".format(elapsed_time))

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
        if i >3:
            break

def test_contiguous_singles():
    particles = particles_fixture()
    for i, tup in enumerate(particles):
        sim_time, surf, entr, arom, shore = tup
        c_arr = sl_contiguous_singles(surf)
        if i >3:
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

if __name__ == "__main__":

#    test_contiguous_doubles()
#    test_contiguous_singles()

    main()




    print("Done __main__")
