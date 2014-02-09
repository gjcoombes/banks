# -*- coding: utf-8 -*-
"""
Module: cl_trial1.py
Created on Sat Oct 26 17:09:11 2013
@author: gav
Description:
Smoke test for opencl routines

"""
### Imports
from __future__ import print_function

import time
### Logging
import logging
logging.basicConfig(level=logging.DEBUG)
debug, info, error = logging.debug, logging.info, logging.error
### Constants

### Classes

### Functions
import pyopencl as cl
import numpy
import numpy.linalg as la

a = numpy.random.rand(5000000).astype(numpy.float32)
b = numpy.random.rand(5000000).astype(numpy.float32)

ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)

mf = cl.mem_flags
a_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a)
b_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=b)
dest_buf = cl.Buffer(ctx, mf.WRITE_ONLY, b.nbytes)

prg = cl.Program(ctx, """
    __kernel void sum(
        __global const float *a,
        __global const float *b,
        __global float *c)
    {
      int gid = get_global_id(0);
      c[gid] = a[gid] + b[gid];
    }
    """).build()

prg.sum(queue, a.shape, None, a_buf, b_buf, dest_buf)

a_plus_b = numpy.empty_like(a)
start_time = time.time()
cl.enqueue_copy(queue, a_plus_b, dest_buf)
elapsed = time.time() - start_time
print(la.norm(a_plus_b - (a+b)), la.norm(a_plus_b))
print("CL Time elapsed {} s".format(elapsed))

ref_start = time.time()
ref = a + b
ref_elapsed = time.time() - ref_start
print("Ref Time elapsed {} s".format(ref_elapsed))
### Tests

if __name__ == "__main__":

    print("Done __main__")
