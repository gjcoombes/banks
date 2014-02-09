# -*- coding: utf-8 -*-
"""
Module: cl_trial2.py
Created on Sat Oct 26 17:23:50 2013
@author: gav
Description:

"""
### Imports
from __future__ import print_function

import pyopencl as cl
import numpy as np

### Logging
import logging
logging.basicConfig(level=logging.DEBUG)
debug, info, error = logging.debug, logging.info, logging.error
### Constants

### Classes

### Functions
def print_info(device):


    s = """
    Name              : {0.name}
    Compute Units     : {0.max_compute_units}
    Local memory size : {0.local_mem_size} B
    Clock Frequency   : {0.max_clock_frequency} MHz
    {0.execution_capabilities}

    """
    print(s.format(device))

def info_summary():
    platforms = cl.get_platforms()

    print("\nNumber of OpenCL platforms:", len(platforms))

    print ("\n-------------------------")

    # Investigate each platform
    for p in platforms:
        # Print out some information about the platforms
        print ("Platform:", p.name)
        print ("Vendor:", p.vendor)
        print ("Version:", p.version)

        # Discover all devices
        devices = p.get_devices()
        print ("Number of devices:", len(devices))

        # Investigate each device
        for d in devices:
            print ("\t-------------------------")
            # Print out some information about the devices
            print ("\t\tName:", d.name)
            print ("\t\tVersion:", d.opencl_c_version)
            print ("\t\tMax. Compute Units:", d.max_compute_units)
            print ("\t\tLocal Memory Size:", d.local_mem_size/1024, "KB")
            print ("\t\tGlobal Memory Size:", d.global_mem_size/(1024*1024), "MB")
            print ("\t\tMax Alloc Size:", d.max_mem_alloc_size/(1024*1024), "MB")
            print ("\t\tMax Work-group Size:", d.max_work_group_size)

            # Find the maximum dimensions of the work-groups
            dim = d.max_work_item_sizes
            print ("\t\tMax Work-item Dims:(", dim[0], " ".join(map(str, dim[1:])), ")")

            print ("\t-------------------------")

        print ("\n-------------------------")
### Tests

if __name__ == "__main__":
    info_summary()
#    platform = cl.get_platforms()[0]
#
#    for d in platform.get_devices():
#        print_info(d)
#    ctx = cl.create_some_context()
#    print(ctx.properties)
#    queue = cl.CommandQueue(ctx)


    print("Done __main__")
