
Notes on optimising the netcdf_to_shape.py

Problem:
Current shapefile generation through matlab is quite slow for large spills and large grids.
Extreme use of RAM and has caused Storm to crash with memory leaks?
Like hours -> days
Some runs have not finished in days, and may not finish.
eg. j0232 sc3 entr prob 10ppb using pcolor ran on ripple in 13 hrs 17 mins or 47,820 seconds
eg. j0232 sc1 surf prob 1gm has not finished in 36 hours (Storm) as shows no encouraging signs.

The usual contour generation can take a long time for large grids with much oil.
Some runs have not finished in days, and may not finish.
The susbsequent shape file generation can take a while if there are many polygons.
The pcolor plot generation is quicker but guarantees many polygons to write as shapes.

Opportunity:
If we had a more efficient method of shapefile generation, we could ensure supply.
We would tax Storm less and move on more fun things like table generation

First Try
=========

I wrote a python module netcdf_to_shape.py which

    - Read the netcdf file for a maxrpob array
    - created individual polygon for each grid cell, the "pixel polygon" method
    - wrote each polygon to disk as it was created

Pros

    - Very low memory use - only the original array

Cons

    - Single threaded, possible slowed by function calls in the inner loop
        + Calculation of vertices for each polygon
        + Dict assignment for each polygon feature

    - Possibly choppy IO use as each polygon is written (not sure if the OS cached writes then flushes periodically)

Results
Something like a 17x speed up in comparison to the pcolor plot and shapefile generation of PlotWinfates
Note no plots were generated for netcdf_to_shape.py, though trial plots have rendered quickly with maplotlib

Next steps
==========

Profiling

Set up an ipython notebook to profile and time func:make_shape_file

Expectations
I expect these areas to be slow
    - func: grid_polygon - simple arithmetic in tight loop
    - write of polygon to fiona collection

Approaches
    - numpy for vectorized calc of grid_polygon
    - cython for to move grid_polygon to c space
    - chunk writing of polygons to collection (n=2048?)






