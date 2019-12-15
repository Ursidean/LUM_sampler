"""
This module is used to calculate the weighted clumpiness error for
the actively modelled land-use classes for two input maps.
"""

from read_map import read_map
import numpy as np
from clumpy_module import clumpiness_index


def w_clu_error(map1, map2, mask, weights, luc, pas, act):
    # Calculate the clumpiness.
    map1_clumpiness = clumpiness_index(map1, mask, luc)
    map2_clumpiness = clumpiness_index(map2, mask, luc)
    clu_error = [0] * luc
    for i in range(0, luc):
        clu_error[i] = abs(map1_clumpiness[i] - map2_clumpiness[i])
    # Generate the error values.
    weighted_clu_error = [0] * act
    for i in range(0, act):
        weighted_clu_error[i] = weights[i] * clu_error[i + pas]
    # Calculate the weighted clumpiness error.
    WCE = sum(weighted_clu_error)
    return WCE

# Module test.
"""
amap_path = ("C:\\Users\\charl\\Dropbox\\PhD\\Journal Papers\\"
             "3. Integrated Automated Calibration\\Results\\"
             "2 Output_maps\\Data_maps\\madrid_2000.asc")
smap_path = ("C:\\Users\\charl\\Dropbox\\PhD\\Journal Papers\\"
             "3. Integrated Automated Calibration\\Results\\"
             "2 Output_maps\\Seeded\\Seed_cal_map_5.rst")
mask_path = ("C:\\Users\\charl\\Dropbox\\PhD\\Journal Papers\\"
             "3. Integrated Automated Calibration\\Results\\"
             "2 Output_maps\\Data_maps\\madrid_mask.asc")
# Read in the maps
smap = read_map(smap_path)
amap = read_map(amap_path)
mask = read_map(mask_path)
# Determine the map properties
map_dimensions = np.shape(amap)
rows = map_dimensions[0]
cols = map_dimensions[1]
luc = np.max(amap) + 1
pas = 1
fea = 5
act = luc - (pas + fea)

weights = [1/11, 1/11, 1/11, 1/11, 2/11, 2/11, 2/11, 1/11]

WCE = w_clu_error(amap, smap, mask, weights, luc, pas, act)
print(WCE)
"""