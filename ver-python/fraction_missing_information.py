#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 19:25:12 2020

@author: fedor.goncharov.ol@gmail.com
"""

import numpy as np

def fraction_missing_pixel(distribution, syst_matrix, avg_scattered, projector):
    """
    Asymptotic value for the Bayesian fraction of missing information.
    
    Possibly this function is not correct. Some theoretical analysis 
    is necessary.
    
    Parameters
    ----------
    distribution : TYPE
        DESCRIPTION.
    syst_matrix : TYPE
        DESCRIPTION.
    projector : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    npixels = np.sqrt(syst_matrix.shape[1]).astype(int)
    
    projector_vec = np.reshape(projector, (npixels**2, 1))
    distribution_vec = np.reshape(distribution, (npixels**2, 1))
    scattered_vec = np.reshape(avg_scattered, (syst_matrix.shape[0],1))
    
    vec_a = syst_matrix.dot(projector_vec)
    pixel_intensity = np.sum(np.multiply(distribution_vec, projector_vec))
    
    normalization = syst_matrix.dot(distribution_vec) + scattered_vec
    
    multiplier_to_sensitivity = np.ones((syst_matrix.shape[0], 1)) - np.divide(pixel_intensity*vec_a, 
                                                      normalization)
    
    main_nominator = np.sum(np.multiply(vec_a, multiplier_to_sensitivity))
    main_denominator = np.sum(vec_a)
    
    return (1-1./(1. + main_nominator/main_denominator))

def fraction_missing_linear(distribution, syst_matrix, projector):
    """
    

    Parameters
    ----------
    distribution : TYPE
        DESCRIPTION.
    syst_matrix : TYPE
        DESCRIPTION.
    projector : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    return None