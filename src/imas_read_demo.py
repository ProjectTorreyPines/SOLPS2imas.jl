#!/bin/env python3

"""
Demonstrates unpacking data from IMAS output.
Needs an h5 file that came from IMAS and that contains SOLPS output.
"""

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import h5py

filename = '/home/eldond/sparc_detach_ctrl/edge_profiles.h5'

imas_data = h5py.File(filename, 'r')
ep = imas_data['edge_profiles']
ggd_num = 1
space_num = 1
points_ref_num = 1
edges_ref_num = 2
cells_ref_num = 3
subset_num = 5
group = 'electrons'
quantity = 'temperature'
num_te_samples_simple = ep['ggd[]&electrons&temperature[]&values_SHAPE']
num_te_samples = ep[f'ggd[]&{group}&{quantity}[]&values_SHAPE'][ggd_num - 1, subset_num - 1, 0]
te = ep[f'ggd[]&{group}&{quantity}[]&values'][ggd_num - 1, subset_num - 1, :num_te_samples]
num_points = ep['grid_ggd[]&space[]&objects_per_dimension[]&object[]&AOS_SHAPE'][ggd_num - 1, space_num - 1, points_ref_num - 1, 0]
num_cells = ep['grid_ggd[]&space[]&objects_per_dimension[]&object[]&AOS_SHAPE'][ggd_num - 1, space_num - 1, cells_ref_num - 1, 0]
points = ep['grid_ggd[]&space[]&objects_per_dimension[]&object[]&geometry'][ggd_num - 1, space_num - 1, points_ref_num - 1, :num_points, :]
r_points = points[:, 0]
z_points = points[:, 1]
boundary_indices = ep['grid_ggd[]&space[]&objects_per_dimension[]&object[]&nodes'][ggd_num - 1, space_num - 1, cells_ref_num - 1].astype(np.int64)
boundary_idx_offset = np.nanmin(boundary_indices)  # I think they index from 1
polygons = []
for cell in range(num_cells):
    sel = boundary_indices[cell, :4] - boundary_idx_offset
    polygons.append(mpl.patches.Polygon([(r_points[sel][i], z_points[sel][i]) for i in range(4)], closed=True))

collection = mpl.collections.PatchCollection(polygons, linewidth=1, norm=None, cmap='viridis')
collection.set_array(te)
collection.set_clim(vmin=0, vmax=np.nanmax(te))

fig, ax = plt.subplots(1)
ax.plot(r_points, z_points, ls='', marker=',', color='gray')
ax.add_collection(collection)
ax.set_xlabel('R / m')
ax.set_ylabel('Z / m')
plt.show()
