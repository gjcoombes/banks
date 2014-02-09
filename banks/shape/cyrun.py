#! /usr/bin/env python 

import pstats, cProfile
import cyshape
from cyshape import mk_shp_v4

cProfile.runctx("cyshape.mk_shp_v4('sum', 'entr_prob_100ppb', cyshape.CONFIG)",
    globals(), locals(), "Profile.prof")

s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()
# if __name__ == "__main__":

#     args = ("sum", "entr_prob_100ppb")
#     season = "sum"
#     nc_var = "entr_prob_100ppb"
#     mk_shp_v4(season, nc_var, cyshape.CONFIG)