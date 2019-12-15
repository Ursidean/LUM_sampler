"""
This module is used to calculate the area weighted clumpiness error for a
given input map. It is developed for testing purposes.
"""

from read_map import read_map
import numpy as np
from clumpy_module import clumpiness_index


def area_weighted_clu_error(map1, map2, mask, luc, pas, act, luc_count):
    # Calculate the clumpiness.
    map1_clumpiness = clumpiness_index(map1, mask, luc)
    map2_clumpiness = clumpiness_index(map2, mask, luc)
    clu_error = [0] * luc
    for i in range(0, luc):
        clu_error[i] = abs(map1_clumpiness[i] - map2_clumpiness[i])
    # Extract the clumpiness error of the active classes.
    act_clu_error = [0] * (act)
    for i in range(0, act):
        act_clu_error[i] = clu_error[i + pas]
    # Calculate the area-weighted clumpiness error.
    AWCE = 0
    for i in range(0, act):
        AWCE = AWCE + act_clu_error[i] * luc_count[i + pas]
    # Active class luc_count.
    act_luc_count = 0
    for i in range(pas, pas + act):
        act_luc_count = act_luc_count + luc_count[i]
    AWCE = AWCE / act_luc_count
    return AWCE

# Module test.
"""
base_path = ("C:\\Users\\charl\\OneDrive\\Documents\\P3_eval\\"
             "CLU_test")
smap_path = base_path + "\\Land use map_2000-Jan-01 00_00_00.rst"
amap_path = base_path + "\\madrid_2000.asc"
mask_path = base_path + "\\madrid_mask.asc"
# Read in the maps
smap = read_map(smap_path)
amap = read_map(amap_path)
mask = read_map(mask_path)
# Determine the map properties
map_dimensions = np.shape(amap)
rows = map_dimensions[0]
cols = map_dimensions[1]
luc = 15
pas = 1
fea = 6
act = luc - (pas + fea)

luc_count = [0] * luc
for i in range(0, rows):
    for j in range(0, cols):
        if mask[i, j] > 0:
            luc_count[amap[i, j]] = luc_count[amap[i, j]] + 1


AWCE = area_weighted_clu_error(amap, smap, mask, luc, pas, act, luc_count)
print(AWCE)
"""