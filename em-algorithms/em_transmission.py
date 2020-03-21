#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 22:49:04 2020

@author: fedor.goncharov.ol@gmail.com
"""

import numpy as np

def em_transmission_convex_nr1():
    """
    Two-dimensional EM-algorithm Convex-NR-1 for X-ray transmission tomography 
    with low counts. A detailed descritption of the algorithm 
    can be found in lecture notes of J. Fessler : Statistical Image 
    Reconstruction methods. Chapter 1.4.5.
    
    Briefly : expectation step - full log-likelyhood for Poissin models
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
    return None

def em_transmission_convex_nr2():
    """
    Two-dimensional EM-algorithm Convex-NR-2 for X-ray transmission tomography 
    with low counts. A detailed descritption of the algorithm 
    can be found in lecture notes of J. Fessler : Statistical Image 
    Reconstruction methods. Chapter 1.4.5.
    
    Briefly : expectation step - full log-likelyhood for Poissin models
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
    
    return None

