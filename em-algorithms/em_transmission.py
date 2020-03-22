#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 22:49:04 2020

@author: fedor.goncharov.ol@gmail.com
"""

import numpy as np

def em_transmission_convex_nr1(sinogram_counts, syst_matrix, source_intensity,
                               avg_scattered, max_iterations, eps, init_point):
    """
    Two-dimensional EM-algorithm Convex-NR-1 for X-ray transmission tomography 
    with low counts. A detailed descritption of the algorithm 
    can be found in lecture notes of J. Fessler : Statistical Image 
    Reconstruction methods. Chapter 1.4.5.
    
    Briefly : expectation step - full log-likelyhood for Poissin model
              maximization step - surrogate via a convex trick with free 
              parameter matrix a_{ij}, i = 1,..,nphi*ntheta, j = 1..dim_im*dim_um.
              
              In NR-1 : Weights a_{ij} are computed online, since they depend
              on current estimate of attenuation.
              
              To optimize a surrogate diagonally-scaled gradient ascent 
              is used (this breaks down global monotone-growth condition)
              
              In view of this, the algorithm has local convergence, not global.
    
              
    
    Input
    ------
    sinogram_counts : matrix of type (nshift, nphi) : sinogram of photon counts
    syst_matrix: matrix of type (nshift*nphi, dim_im * dim_im) : system matrix 
    source_intensity : matrix of type (nshift, nphi) : average flux 
        of photons per LOR without the object
    avg_scattered : matrix of type (nshift, nphi) : average flux of 
        scattered photons 
    max_iterations : integer : maximal number of iterations 
    err : float : bound for l2-relative error
    init_point : matrix of type (dim_im, dim_im) : initial guess for attenuation

    Returns
    -------
    matrix of type (dim_im, dim_im) -- reconstructed image 

    """
    
    npixels = np.sqrt(syst_matrix.shape[0]).astype(int)
    nshift = sinogram_counts.shape[0]
    nphi = sinogram_counts.shape[1]
    
    vector_counts = np.reshape(sinogram_counts, (nphi*nshift, 1))
    vector_scattered = np.reshape(avg_scattered, (nphi*nshift,1))
    vector_intensity = np.reshape(source_intensity, (nshift*nphi, 1))
    current_attenuation = np.reshape(init_point, (npixels*npixels, 1)) # vector for attenuation map
    
    iteration = 0
    err = np.inf	
	
    while (iteration < max_iterations) and (err > eps):
        vector_l = syst_matrix.dot(current_attenuation)
        vector_bl = np.multiply(vector_intensity, np.exp((-1)*vector_l))
        
        nominator_arg = np.multiply(np.ones((nshift*nphi, 1)) - np.divide(vector_counts, 
                                                                          vector_bl + vector_scattered), 
                                    vector_bl)
        denominator_arg = np.multiply(vector_l, vector_bl)
        
        
        vector_nominator = syst_matrix.T.dot(nominator_arg)
        vector_denominator = syst_matrix.T.dot(denominator_arg)
        
        current_attenuation_new = current_attenuation + np.multiply(current_attenuation, 
                                               np.divide(vector_nominator, vector_denominator))
        
        err = np.linalg.norm(current_attenuation_new - current_attenuation, 'fro') / np.linalg.norm(current_attenuation, 'fro')
        iteration += 1
        current_attenuation = current_attenuation_new
        
        # clean current variables
        del(current_attenuation_new)
        del(vector_l)
        del(vector_bl)
        del(nominator_arg)
        del(denominator_arg)
        del(vector_nominator)
        del(vector_denominator)
        print(f'Iteration {iteration},  relative l2-error : {err} \n')      
        
    return np.reshape(current_attenuation, (npixels, npixels))


def em_transmission_convex_nr2(sinogram_counts, syst_matrix, transmission_intensity,
                               avg_scattered, max_iterations, eps, init_point):
    """
    Two-dimensional EM-algorithm Convex-NR-2 for X-ray transmission tomography 
    with low counts. A detailed descritption of the algorithm 
    can be found in lecture notes of J. Fessler : Statistical Image 
    Reconstruction methods. Chapter 1.4.5.
    
    Briefly : expectation step - full log-likelyhood for Poissin model
              maximization step - surrogate via a convex trick with free 
              parameter matrix a_{ij}, i = 1,..,nphi*ntheta, j = 1..dim_im*dim_um.
              
              In NR-2 : Weights a_{ij} are precomputed online, since they depend
              only on system matrix.
              
              
              To optimize a surrogate diagonally-scaled gradient ascent 
              is used (this breaks down global monotone-growth condition)
              
              In view of this, the algorithm has local convergence, not global.
              
    
    Input
    ------
    sinogram_counts : matrix of type (nshift, nphi) : sinogram of photon counts
    syst_matrix: matrix of type (nshift*nphi, dim_im * dim_im) : system matrix 
    transmission_intensity : matrix of type (nshift, nphi) : average flux 
        of photons per LOR without the object
    avg_scattered : matrix of type (nshift, nphi) : average flux of 
        scattered photons 
    max_iterations : integer : maximal number of iterations 
    err : float : bound for l2-relative error
    init_point : matrix of type (dim_im, dim_im) : initial guess for attenuation

    Returns
    -------
    matrix of type (dim_im, dim_im) -- reconstructed image
    """
    npixels = np.sqrt(syst_matrix.shape[0]).astype(int)
    nshift = sinogram_counts.shape[0]
    nphi = sinogram_counts.shape[1]
    
    vector_counts = np.reshape(sinogram_counts, (nphi*nshift, 1))
    vector_scattered = np.reshape(avg_scattered, (nphi*nshift,1))
    vector_intensity = np.reshape(transmission_intensity,(nshift*nphi, 1))
    current_attenuation = np.reshape(init_point, (npixels*npixels, 1)) # vector for attenuation map
    
    # precomputation of vector a
    vector_a = syst_matrix.sum(axis=1, keepdims=True)
    
    iteration = 0
    err = np.inf	
	
    while (iteration < max_iterations) and (err > eps):
        vector_l = syst_matrix.dot(current_attenuation)
        vector_bl = np.multiply(vector_intensity, np.exp((-1)*vector_l))
        
        nominator_arg = np.multiply(np.ones((nshift*nphi, 1)) - np.divide(vector_counts, 
                                                                          vector_bl + vector_scattered), 
                                    vector_bl)
        denominator_arg = np.multiply(vector_bl, vector_a)
        
        
        vector_nominator = syst_matrix.T.dot(nominator_arg)
        vector_denominator = syst_matrix.T.dot(denominator_arg)
        
        current_attenuation_new = current_attenuation + np.divide(vector_nominator, vector_denominator)
        
        err = np.linalg.norm(current_attenuation_new - current_attenuation, 'fro') / np.linalg.norm(current_attenuation, 'fro')
        iteration += 1
        current_attenuation = current_attenuation_new
        
        # clean current variables
        del(current_attenuation_new)
        del(vector_l)
        del(vector_bl)
        del(nominator_arg)
        del(denominator_arg)
        del(vector_nominator)
        del(vector_denominator)
        print(f'Iteration {iteration},  relative l2-error : {err} \n')      
        
    return np.reshape(current_attenuation, (npixels, npixels))

