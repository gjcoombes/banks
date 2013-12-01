# -*- coding: utf-8 -*-
"""
Module: compare_results.py
Created on Sat Nov 30 13:38:04 2013
@author: gav
Description:
Compare pythresh and stochastic_threshold results
"""
### Imports
from __future__ import print_function

import numpy as np
import numpy.ma as ma
import tables as tb
from netCDF4 import Dataset
import matplotlib.pyplot as plt

from pythresh_demo_dense import ga_grid_cell_areas
# Usage cell_areas = ga_grid_cell_areas(grid_fp=grid_fp)

### Logging
import logging
logging.basicConfig(level=logging.DEBUG)
debug, info, error = logging.debug, logging.info, logging.error
### Constants

### Classes

### Functions
def h5_data(h5_fp, compare_var="max_mass", grid_fp=None):
    """Return an array of the compared variable from the pythresh file"""
    v_dict = {
        "max_mass": "max_mass",
        "min_time": "min_time",
        "exceedance": "exceedance",
        }
    v = v_dict[compare_var]

    grp = "/result"
    with tb.open_file(h5_fp, "r") as h5_src:
        print(h5_src)
        arr = h5_src.get_node(grp, v).read()

    if compare_var == "max_mass":
        cell_areas = ga_grid_cell_areas(grid_fp=grid_fp)
        ma_arr = ma.masked_equal(arr, 0)
        conc_arr = ma_arr / cell_areas[:, np.newaxis]
        out_arr = conc_arr
    elif compare_var == "exceedance":
        out_arr = ma.masked_equal(arr, 0)
    elif compare_var == "min_time":
        out_arr = ma.masked_equal(arr, 1e9)
    return out_arr

def nc_data(nc_fp, compare_var="max_mass"):
    """Return an array of the compared variable from the stoch_thresh file"""
    v_dict = {
        "max_mass": "conc_max",
        "min_time": "surfthreshtime",
        "exceedance": "concprob",
        }
    v = v_dict[compare_var]
    with Dataset(nc_fp) as src:
        print(src)
        print(src.variables[v])
        if compare_var == "max_mass":
            arr = src.variables[v][:,:]
            out_arr = arr * 1e-3 # scale to T/m2
        elif compare_var in {"exceedance",}:
            out_arr = src.variables[v][0,:,:]
        elif compare_var in {"min_time",}:
            out_arr = src.variables[v][0,:,:]
    return out_arr

def comparison_plot(h5_arr, nc_arr):
    """Return a matplotlib figure"""
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    im1 = ax1.imshow(h5_arr, origin="upper", interpolation="nearest")
    im2 = ax2.imshow(nc_arr, origin="lower", interpolation="nearest")
#    cb1 = fig.colorbar(im1)
#    cb2 = fig.colorbar(im2)
    return fig

def comparison_stats(h5_arr, nc_arr):
    """Compare the vailes of the arrays"""

    print("Mean of h5 data is {}, and max is {}, num of non-zero is {}".format(
        ma.mean(h5_arr), ma.max(h5_arr), np.nonzero(h5_arr)[0].shape))
    print("Mean of nc data is {}, and max is {}, num of non-zero is {}".format(
        ma.mean(nc_arr), ma.max(nc_arr), np.nonzero(nc_arr)[0].shape))
### Tests

if __name__ == "__main__":
    project = "j0272_sc1_020"
    compare_var = "min_time"

    fp_dict = {
        "j0272_sc1_001": {
            "nc_fp":   r"J:\data\pythresh_trials\netcdf\J0272_SC1_SURF_MUTI_255M3HR_SUM_file_001.nc",
            "h5_fp":   r"J:\data\pythresh_trials\hdf5\j0272_data.h5",
            "grid_fp": r"J:\data\pythresh_trials\grids\MidWA_NWS_1km.DEP",
        },
        "j0272_sc1_010": {
            "nc_fp":   r"J:\data\pythresh_trials\netcdf\J0272_SC1_SURF_MUTI_255M3HR_SUM_spurious_file_010.nc",
            "h5_fp":   r"J:\data\pythresh_trials\hdf5\j0272_sc1_sum_010_data.h5",
            "grid_fp": r"J:\data\pythresh_trials\grids\MidWA_NWS_1km.DEP",
        },
        "j0272_sc1_016": {
            "nc_fp":   r"J:\data\pythresh_trials\netcdf\J0272_SC1_SURF_MUTI_255M3HR_SUM_spurious_file_016.nc",
            "h5_fp":   r"J:\data\pythresh_trials\hdf5\j0272_sc1_sum_016_data.h5",
            "grid_fp": r"J:\data\pythresh_trials\grids\MidWA_NWS_1km.DEP",
        },
        "j0272_sc1_020": {
            "nc_fp":   r"J:\data\pythresh_trials\netcdf\J0272_SC1_SURF_MUTI_255M3HR_SUM_spurious_file_020.nc",
            "h5_fp":   r"J:\data\pythresh_trials\hdf5\j0272_sc1_sum_020_data.h5",
            "grid_fp": r"J:\data\pythresh_trials\grids\MidWA_NWS_1km.DEP",
        },
        "j0266_sc2_001": {
            "nc_fp":   r"J:\data\pythresh_trials\netcdf\J0266_SC2A_SUBS_18962M3_LAM_Q1_file_001.nc",
            "h5_fp":   r"J:\data\pythresh_trials\hdf5\j0266_sc2_sum_001_data.h5",
            "grid_fp": r"J:\data\pythresh_trials\grids\Penguin_1000m.DEP",
        },
    }

    fps = fp_dict[project]

    h5_arr = h5_data(fps['h5_fp'], compare_var, fps['grid_fp'])
    nc_arr = nc_data(fps['nc_fp'], compare_var)
    fig = comparison_plot(h5_arr, nc_arr)
    stats = comparison_stats(h5_arr, nc_arr)

    plt.show()
    print("Done __main__")
