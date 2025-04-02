#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 10:23:58 2024

@author: max
"""

import numpy as np
from scipy.io import loadmat
from nilearn import plotting
import matplotlib.pyplot as plt

p_thresh = 0.001

data = loadmat('../results/spatial_plots/fig4_matrices_rest_movie_coords.mat')

varx_R_rest = data['varX_Rvalue_rest']
varx_R_rest = varx_R_rest - (np.eye(varx_R_rest.shape[0]) * np.diag(varx_R_rest))
varx_R_movie = data['varX_Rvalue_dme']
varx_R_movie = varx_R_movie - (np.eye(varx_R_movie.shape[0]) * np.diag(varx_R_movie))

varx_p_rest = data['varX_pvalue_rest']
varx_p_rest = varx_p_rest - (np.eye(varx_p_rest.shape[0]) * np.diag(varx_p_rest))
varx_p_movie = data['varX_pvalue_dme']
varx_p_movie = varx_p_movie - (np.eye(varx_p_movie.shape[0]) * np.diag(varx_p_movie))

coordinates = data['coords']

# Symmetric version
varx_R_rest = np.triu(varx_R_rest) + np.triu(varx_R_rest).T
varx_R_movie = np.triu(varx_R_movie) + np.triu(varx_R_movie).T

varx_p_rest = np.triu(varx_p_rest) + np.triu(varx_p_rest).T
varx_p_movie = np.triu(varx_p_movie) + np.triu(varx_p_movie).T

idx_sig = np.logical_and(varx_p_rest < p_thresh, varx_p_movie < p_thresh)

varx_R_rest[np.invert(idx_sig)] = np.nan
varx_R_movie[np.invert(idx_sig)] = np.nan

plt.rcParams.update({'font.size': 18})

plot_range = data['plot_range'][0]

#%% Rest
plotting.plot_connectome(
    varx_R_rest,
    coordinates,
    node_color='k',
    edge_threshold=0,
    edge_vmin=0, 
    edge_vmax=plot_range[1],
    node_size=2,
    edge_cmap='Reds',
    edge_kwargs={'linewidth':2},
    colorbar=True,
    display_mode='z',
)

plt.gcf().set_figwidth(5)
plt.gcf().set_figheight(5)

plt.savefig('../results/figures/fig4_rest_connectome_plot_R_LFP.png', dpi=300)

#%% Movie
plotting.plot_connectome(
    varx_R_movie,
    coordinates,
    node_color='k',
    edge_threshold=0,
    edge_vmin=0, 
    edge_vmax=plot_range[1],
    node_size=2,
    edge_cmap='Reds',
    edge_kwargs={'linewidth':2},
    colorbar=True,
    display_mode='z',
)

plt.gcf().set_figwidth(5)
plt.gcf().set_figheight(5)

plt.savefig('../results/figures/fig4_movie_connectome_plot_R_LFP.png', dpi=300)

#%% Movie - Rest
plotting.plot_connectome(
    varx_R_movie - varx_R_rest,
    coordinates,
    node_color='k',
    edge_threshold=0,
    edge_vmin=-plot_range[1], 
    edge_vmax=plot_range[1],
    node_size=2,
    edge_cmap='seismic',
    edge_kwargs={'linewidth':2},
    colorbar=True,
    display_mode='z',
)

plt.gcf().set_figwidth(5)
plt.gcf().set_figheight(5)

plt.savefig('../results/figures/fig4_diff_connectome_plot_R_LFP.png', dpi=300)
