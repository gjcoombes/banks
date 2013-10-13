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
Changes:
    - Inline grid_edge_length and diff a single pair (confirmed all equal)
    - Move the slice of threshold threshold arr outside of loop
    - Unpack lat, lon idx with the loop
