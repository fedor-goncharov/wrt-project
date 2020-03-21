#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 15:04:09 2020

@author: fedor.goncharov.ol@gmail.com
"""

import numpy as np

def EM_algorithm_emission2d(sinogram_counts, syst_matrix, avg_scattered, max_iterations, eps, init_point):
    ''' 
  two-dimensional EM algorithm of Shepp-Vardi without any regularization
  data is assumed to be given in sinorgam form for parallel beam geometry

  Input 
    sinogram_counts : matrix of type (nshift, nphi)    : sinogram of photon counts 
    syst_matrix     : matrix of type (dim_sin=nshift*nphi, dim_im) : system matrix for the Radon transform
	avg_scattered   : matrix of type (nshift, nphi)    : average number of scattered photons 
    max_iterations  : integer                          : maximal number of iterations
    eps             : float                            : error level for the stopping rule (max norm is used)
	init_point      : matrix of type (dim_im, dim_im)  : starting point for the iterative algorithm
    
  Output
    image : matrix of size size(syst_matrix, 2)
	'''

    npixels = np.sqrt(syst_matrix.shape[0]).astype(int)
    nshift = sinogram_counts.shape[0]
    nphi = sinogram_counts.shape[1]
	
    current_density = np.reshape(init_point, (npixels**2, 1)) # init lambda_0
    scattered = np.reshape(avg_scattered, (nshift*nphi, 1)) # init r
    current_counts = syst_matrix.dot(current_density) + scattered;
    data_counts = np.reshape(sinogram_counts, (nshift*nphi, 1))
    iteration = 0
    err = np.inf	
	
    denominator_a = np.sum(syst_matrix, 0, keepdims=True).T # vector of denominators in EM-algorithm
	
    while (iteration < max_iterations) and (err > eps):
        current_density_new = np.multiply(current_density, 
                                          np.divide(syst_matrix.T.dot(np.divide(data_counts, current_counts)), 
                                                    denominator_a))
        current_counts = syst_matrix.dot(current_density_new) + scattered;
        err = np.linalg.norm(current_density_new - current_density, 'fro') / np.linalg.norm(current_density, 'fro')
        iteration += 1
        current_density = current_density_new
        print(f'Iteration {iteration},  relative l2-error : {err} \n')
		
    return np.reshape(current_density, (npixels, npixels))


def EM_algorithm_emsission2d_mlem3(sinogram_counts, syst_matrix, avg_scattered, m_coeffs, 
                                   max_iterations, eps, init_point):
    ''' 
  two-dimensional EM algorithm ML-EM-3 without any regularization
  data is assumed to be given in sinorgam form for parallel beam geometry
  (see lecture notes of J. Fessler : Statistical image reconstruction methods for transmission
  tomography, section 1.9.2)

  Input 
    sinogram_counts : matrix of type (nshift, nphi)    : sinogram of photon counts 
    syst_matrix     : matrix of type (dim_sin=nshift*nphi, dim_im) : system matrix for the Radon transform
	avg_scattered   : matrix of type (nshift, nphi)    : average number of scattered photons 
    m_coeffs        : matrix of type (dim_im, dim_m)   : coefficients for accelerating convergence
    max_iterations  : integer                          : maximal number of iterations
    eps             : float                            : error level for the stopping rule (max norm is used)
	init_point      : matrix of type (dim_im, dim_im)  : starting point for the iterative algorithm
    
  Output
    image : matrix of size size(syst_matrix, 2)
  
    If m_coeffs will not satisfy monote growth condition : r - sys_matrix*m > 0, then None will be returned. 
	'''
    npixels = np.sqrt(syst_matrix.shape[0]).astype(int)
    nshift = sinogram_counts.shape[0]
    nphi = sinogram_counts.shape[1]
	
    current_density = np.reshape(init_point, (npixels**2, 1)) # init lambda_0
    m_coeffs_vector = np.reshape(m_coeffs, (npixels**2, 1))  
    scattered = np.reshape(avg_scattered, (nshift*nphi, 1)) # init r
    
    # check m_coeffs satisfy monotone growth condition
    test = scattered - syst_matrix.dot(m_coeffs_vector)
    if not (test > 0).all():
        print("Coefficients M are too big and do not satisfy monotone growrth condition. Abort\n")
        return None
    
    current_counts = syst_matrix.dot(current_density) + scattered;
    data_counts = np.reshape(sinogram_counts, (nshift*nphi, 1))
    iteration = 0
    err = np.inf	
	
    denominator_a = np.sum(syst_matrix, 0, keepdims=True).T # vector of denominators in EM-algorithm
	
    while (iteration < max_iterations) and (err > eps):
        current_density_new = np.maximum(np.multiply(current_density + m_coeffs_vector, 
                                          np.divide(syst_matrix.T.dot(np.divide(data_counts, current_counts)), 
                                                    denominator_a)) - m_coeffs_vector, 0)
        current_counts = syst_matrix.dot(current_density_new) + scattered;
        err = np.linalg.norm(current_density_new - current_density, 'fro') / np.linalg.norm(current_density, 'fro')
        iteration += 1
        current_density = current_density_new
        print(f'Iteration {iteration}, relative l2-error : {err} \n')
		
    return np.reshape(current_density, (npixels, npixels))




