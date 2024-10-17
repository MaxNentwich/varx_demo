#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 04:01:07 2024

@author: max
"""

#%% Can use this to get any feature map from neuromaps and average in parcels from .annot files 

#%% Load T1w/T2w data
import sys
import numpy as np
from neuromaps.datasets import available_annotations
import pandas as pd
import matplotlib.pyplot as plt

for annotation in available_annotations(source='hcps1200'):
    print(annotation)
    
sys.path.append('../src')
import echo_utils
    
# Files with vertex coordinates from freesurfer (included in MNE)
lh_annot_file = '../data/lh.aparc.annot'
rh_annot_file = '../data/rh.aparc.annot'

#%% Load myelination map
myelin_data = pd.read_excel('../data/myelin_lh_parcels_aparc.xlsx')
myelin_lh_parcels = myelin_data.annot_val.values

#%% Plot incoming vs outgoing connections
fig_path = '../results/figures/fig7_%s.png'

varx_dir_aparc = pd.read_excel('../data/varx_dir_aparc.xlsx')
varx_dir = varx_dir_aparc.varx_dir.values
varx_parcel = varx_dir_aparc.parcel_name.values

varx_dir_sort = np.empty(len(myelin_data))

# Order values
for p in range(len(myelin_data)):
    
    parcel = str(myelin_data.parcel_name[p]).replace('b\'', '')
    parcel = parcel.replace('\'', '')
    
    idx_parcel = np.in1d(varx_parcel, parcel)
    
    if np.sum(idx_parcel) == 0:
        varx_dir_sort[p] = np.nan
    else:
        varx_dir_sort[p] = varx_dir[idx_parcel]

c_range = np.max([abs(np.nanmin(varx_dir_sort)), abs(np.nanmax(varx_dir_sort))])

plt.rcParams.update({'font.size': 26})

plt.figure(figsize=(6,5))
echo_utils.plot_MMP(varx_dir_sort, annot_file=lh_annot_file, save_file=fig_path%'direction_desikan_killiany', 
                    minmax=[-c_range, c_range], bp=0, title=r'$R-R^T$',
                    cmap='seismic')

#%% Plot atlas of T1w/T2w values for comparison
c_min = np.min(myelin_lh_parcels[myelin_lh_parcels != 0])
c_max = abs(np.max(myelin_lh_parcels))

plt.figure(figsize=(6,5))
echo_utils.plot_MMP(myelin_lh_parcels, annot_file=lh_annot_file, save_file=fig_path%'myelin_desikan_killiany', 
                    minmax=[c_min, c_max], bp=0, title='T1w/T2w', cmap='Reds')

