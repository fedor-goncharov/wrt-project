#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 15:04:09 2020

@author: fedor.goncharov.ol@gmail.com
"""

import numpy as np

def em_emission2d(sinogram_counts, syst_matrix, avg_scattered, max_iterations, eps, init_point):
    ''' 
  two-dimensional EM algorithm of Shepp-Vardi without any regularization
  data is assumed to be given in sinorgam form for parallel beam geometry

  Input 
    sinogram_counts : matrix of type (nshift, nphi)    : sinogram of photon counts 
    syst_matrix     : matrix of type (dim_sin, dim_im) : system matrix for the Radon transform
	avg_scattered   : matrix of type (nshift, nphi)    : average number of scattered photons 
    max_iterations  : integer                          : maximal number of iterations
    eps             : float                            : error level for the stopping rule (max norm is used)
	init_point      : matrix of type (dim_im, dim_im)  : starting point for the iterative algorithm
    
  Output
    image : matrix of size size(syst_matrix, 2)
	'''

    npixels_dim = syst_matrix.shape[0]
    nshift = sinogram_counts.shape[0]
    nphi = sinogram_counts.shape[1]
	
    current_density = np.reshape(init_point, (npixels_dim**2, 1)) # init lambda_0
    scattered = np.reshape(avg_scattered, (nshift*nphi, 1)) # init r
    current_counts = syst_matrix.dot(current_density) + scattered;
    data_counts = np.reshape(sinogram_counts, (nshift*nphi, 1))
    iteration = 0
    err = np.inf	
	
    denominator_a = np.sum(syst_matrix, 1) # vector of denominators in EM-algorithm
	
    while (iteration < max_iterations) and (err > eps):
        current_density_new = np.multiply(current_density, 
                                          np.divide(syst_matrix.T.dot(np.divide(data_counts, current_counts)), 
                                                    denominator_a))
        current_counts = syst_matrix.dot(current_density_new) + scattered;
        err = np.linalg.norm(current_density_new - current_density, 'fro') / np.lingalg.norm(current_density, 'fro')
        iteration += 1
        current_density = current_density_new
		
    return np.reshape(current_density, (npixels_dim, npixels_dim))



