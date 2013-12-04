# -*- coding: utf-8 -*-
"""
Module: cythresh_harness.py
Created on Wed Dec 04 20:56:27 2013
@author: gav
Description:
Wrapper for cythresh

"""
### Imports
from __future__ import print_function
import os.path as osp
import numpy as np
import matplotlib.pyplot as plt

import cythresh

### Logging
import logging
logging.basicConfig(level=logging.DEBUG)
debug, info, error = logging.debug, logging.info, logging.error
### Constants

### Classes

### Functions

### Tests

if __name__ == "__main__":

    project_dir = r"J:\data\j0267_nv_remodel"
    stem = "J0267_SC3_SBED_LEAK_TRA_001"
    grid_fn = "VanGogh_800m.DEP"

    tr3_fp = osp.join(project_dir, "modelout", stem + ".tr3" )
    lu3_fp = osp.join(project_dir, "modelout", stem + ".lu3" )
    grid_fp = osp.join(project_dir, "grids", grid_fn)
    h5_fp = osp.join(project_dir, 'hdf5', stem + ".h5")
#    surf_thresholds_gm2 = np.arange(1e-6, 11e-6, 1e-6)
    surf_thresholds_gm2 = np.array([1e-7])

    result = cythresh.main(tr3_fp, lu3_fp, grid_fp, h5_fp, surf_thresholds_gm2)
    exceedance_arr = result['exceedance_arr']
    im = exceedance_arr[0,:,:]
    plt.imshow(im, origin="upper", interpolation="nearest")
    plt.show()
    print("Done __main__")
