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

data = loadmat('../results/spatial_plots/fig4_matrices_rest_movie_coords.mat')

varx_R_rest = data['varX_Rvalue_rest']
varx_R_rest = varx_R_rest - (np.eye(varx_R_rest.shape[0]) * np.diag(varx_R_rest))
varx_R_movie = data['varX_Rvalue_dme']
varx_R_movie = varx_R_movie - (np.eye(varx_R_movie.shape[0]) * np.diag(varx_R_movie))
coordinates = data['coords']

# Symmetric version
varx_R_rest = np.triu(varx_R_rest) + np.triu(varx_R_rest).T
varx_R_movie = np.triu(varx_R_movie) + np.triu(varx_R_movie).T

plt.rcParams.update({'font.size': 26})

#%% Rest
plotting.plot_connectome(
    varx_R_rest,
    coordinates,
    edge_threshold=0.05,
    edge_vmin=0, 
    edge_vmax=0.35,
    node_size=20,
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
    edge_threshold=0.05,
    edge_vmin=0, 
    edge_vmax=0.35,
    node_size=20,
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
    edge_threshold=0.05,
    edge_vmin=0, 
    edge_vmax=0.35,
    node_size=20,
    edge_cmap='Reds',
    edge_kwargs={'linewidth':2},
    colorbar=True,
    display_mode='z',
)

plt.gcf().set_figwidth(5)
plt.gcf().set_figheight(5)

plt.savefig('../results/figures/fig4_diff_connectome_plot_R_LFP.png', dpi=300)
