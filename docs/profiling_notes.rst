Profiling of make_shpe_file in netcdf_to_shape.py
=================================================

Timing runs for sc8 sum entr_prob_100ppb
Wall time 11.7s
Now removing print statements still 11.7s
Now moving to another file for iterative profiling - profile_shape.py

Plan of Action
--------------

1.  Line profile for silly inclusions to loop
2.  Inline polygon code - remove shapely
3.  Replace mask and argwhere with non_zero numpu funcs
3.  Replace grid_polygon computation with array lookup for corners
4.  Prepare for Cython
    a.  Numpy array version
    b.  Unrolled loop version

Line profiling
--------------
With

import profile_shape
from profile_shape import make_shape_file

nc_var = "entr_prob_100ppb"
season = "sum"
cfg = profile_shape.CONFIG

%lprun -f make_shape_file make_shape_file(season, nc_var, cfg)

Highlights:
    Overall run time from 11.7s -> 14.2s - ie profiling overhead
    20% on the mean diff line
    29% on the threshold_arr slice (inside loop doh!)
    19% on the coord lookup for tuple!

mk_shp_v2
---------

Changes:
    - Inline grid_edge_length and diff a single pair (confirmed all equal)
    - Move the slice of threshold threshold arr outside of loop
    - Unpack lat, lon idx with the loop

Results:
    %timeit mk_shp_v2(season, nc_var, cfg)

    5.66s per loop
    so 11.7 -> 5.66 or 2.06x speedup

Line_profile

    - 16% logical mask 16%
    - 15% argwhere  15%
    - 38% lookup of lat lon from idx array - 38% (tightloop)
    - 7% grid_polygon
    - 9% update of fiona record
    - 11% write of record

Next:
    Ensure local copy of lon/lat arrays
    Inline record update - mk_shp_v3

mk_shp_v3
---------
original timeit 5.67s

Now change lon and lat to local numpy arrays
Now timeit 5.73s
Reshape array and stack in columns
5.72s

Time for cython array lookups

mk_shp_v4 from cyshape.py
-------------------------
stly typed

Refactor into functions in prep for cython
Improved time of 5.09s
Changes
    - made corner arrays for lookup
    - using index-range loops in prep for c-speed
    - added new function compute_corners
    - inlined the polygon vert

Time for Hello World cython
Got to hello world after some trouble with compilers
Follow directions here
http://blog.victorjabur.com/2011/06/05/compiling-python-2-7-modules-on-windows-32-and-64-using-msvc-2008-express/

Compiled cyshape in to cyshape.pyd using cython and msvc compiler
Working code as compiled
% timeit in 5.12s

Typed variables inside compute_corners - 5.14s
Typed variables inside prob_threshold_array - 5.14s
cimport numpy                               - 5.16s
Mostly type compute_and_write_records but still - 5.08s
