#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:22:09 2024

@author: max

Map T1w/T2w ratio on parcels in DK atlas
"""
import numpy as np
import pandas as pd
from neuromaps.datasets import fetch_annotation
from neuromaps import transforms
import nibabel as nib

#%% Myeling maps T1w/T2w

# Load the full high resolution map from neuromaps
myelin = fetch_annotation(source='hcps1200', desc='myelinmap')

# Convert to freesurfer surface
myelin_fsaverage_164k = transforms.fslr_to_fsaverage(myelin, '164k')
lh_fs164k_data = myelin_fsaverage_164k[0].agg_data()
rh_fs164k_data = myelin_fsaverage_164k[1].agg_data()

# Load vertex coordinates from freesurfer (included in MNE)
lh_annot_file = '../data/lh.aparc.annot'
rh_annot_file = '../data/rh.aparc.annot'

# Parse data
lh_annot_data = nib.freesurfer.io.read_annot(lh_annot_file)
rh_annot_data = nib.freesurfer.io.read_annot(rh_annot_file)

lh_voxel_id = lh_annot_data[0]
rh_voxel_id = rh_annot_data[0]

lh_parcel_name = np.asarray(lh_annot_data[2])
rh_parcel_name = np.asarray(rh_annot_data[2])

lh_parcel_id = np.arange(len(lh_parcel_name))
rh_parcel_id = np.arange(len(lh_parcel_name))

lh_parcels = pd.DataFrame({'parcel_id':lh_parcel_id, 'parcel_name':lh_parcel_name})
rh_parcels = pd.DataFrame({'parcel_id':rh_parcel_id, 'parcel_name':rh_parcel_name})

# Map average myelination values to each parcel
myelin_lh_parcels = np.array([np.nan] * len(lh_parcels))

for l in range(len(lh_parcels)):
    myelin_lh_parcels[l] = np.mean(lh_fs164k_data[lh_voxel_id == lh_parcels.parcel_id[l]])
    
myelin_rh_parcels = np.array([np.nan] * len(rh_parcels))

for r in range(len(rh_parcels)):
    myelin_rh_parcels[r] = np.mean(rh_fs164k_data[rh_voxel_id == rh_parcels.parcel_id[r]])

myelin_lh_parcels[np.isnan(myelin_lh_parcels)] = 0
myelin_rh_parcels[np.isnan(myelin_rh_parcels)] = 0

#%% Save data as a table 
myelin_lh_parcels = pd.DataFrame({'parcel_name':lh_parcel_name, 't1wt2w':myelin_lh_parcels})
myelin_rh_parcels = pd.DataFrame({'parcel_name':rh_parcel_name, 't1wt2w':myelin_rh_parcels})

myelin_lh_parcels.to_excel('../data/myelin_lh_parcels_aparc.xlsx')
myelin_rh_parcels.to_excel('../data/myelin_rh_parcels_aparc.xlsx')