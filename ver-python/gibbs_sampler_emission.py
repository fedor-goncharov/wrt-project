#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 20:30:56 2020

@author: fedor.goncharov.ol@gmail.com
"""

import numpy as np

def gibbs_sampler_emission_gamma(sinogram_counts, observation_time, 
                                 syst_matrix, avg_scattered,
                                 gamma_prior_params, 
                                 init_point, 
                                 burn_in_size, sample_size):
    """
    Gibbs sampler for emission tomography with 
    separable Gamma-prior (i.e., without spatial regularizaiton). 
    That is, in each pixel we assume Gamma prior distribution (alpha, beta).
    
    Description
    -------
    Sampling using Gibbs sampler from the posterior distribution p(lambda | y) 
    the for isotope distribution. Auxiliary variable n_ij is introduced, 
    which corresponds to number of photons generated in j-th pixel and 
    reached i-th detector (or i-th LOR). 
    
    Each sample is an image generated in two steps: 
        
        Step 1. Backprojection of observed data n_ij ~ p(n_ij | y, lambda)
        p(n_ij, j | y, lambda) - multinomial
        
        
        Step 2. Sampling of image lambda ~ p(lambda | n_ij, y)
        
    Remark
    -------
    This is a direct Gibbs sampler (without any spatial regularization),
    nor acceleration (see Data Augmentation techniques).
    
    Input
    -------
    sinogram_counts : 
        matrix of type (nshift, nphi) : photon counts along LORs
    
    observation_time : 
        scalar : time for which the given photon counts were generated
    
    syst_matrix : 
        matrix of type (nshift*nphi, dim_im*dim_im) : system matrix of observation 
                                                                data
    avg_scattered : 
        matrix of type (nshift, nphi) : average number of scattered photons per LOR
        
    gamma_prior_params : 
        tuple of type (float, float) : parameters (shape, rate) for Gamma-distribution prior
        
    init_point : 
        matrix of type (dim_im, dim_im) : initial guess for the distribution
        
    burn_in_size : 
        integer : number of iterations for burn-in 
        
    sample_size : 
        integer : number of samples
    
    Returns
    -------
        matrix of type (dim_im, dim_im, sample_size) containing all required samples

    """
    npixels = np.sqrt(syst_matrix.shape[1]).astype(int)
    nshift = sinogram_counts.shape[0]
    nphi = sinogram_counts.shape[1]
    prior_shape = gamma_prior_params[0] # alpha - shape parameter for prior Gamma distribution
    prior_rate = gamma_prior_params[1] # beta - rate parameter for prior Gamma distribution
    
    # vectorize input 
    init_point_vec = np.reshape(init_point, (npixels*npixels, 1))
    sinogram_counts_vec = np.reshape(sinogram_counts, (nshift*nphi, 1))
    avg_scattered_counts_vec = np.reshape(avg_scattered, (nshift*nphi, 1))
    
    # precomputations, reserve output storage 
    sensitivity_array = syst_matrix.sum(axis=0, keepdims=True) # returns a row vector
    output_array = np.zeros((npixels, npixels, sample_size))
    
    # initialize variable in a cycle for Gibbs sampler
    current_density_vec = init_point_vec
    
    for iteration in range(burn_in_size + sample_size):
            print(f'Iteration {iteration}\n')
        # Step 1 : backprojection of observed data : n_ij ~ p(n_ij, y, lambda) (multinomial)
            backprojection_data = np.zeros((nshift*nphi, npixels*npixels))
    
        # assemble probability matrix
            multinomial_prob_matrix_denominator = syst_matrix.dot(current_density_vec) + avg_scattered_counts_vec
            multinomial_prob_matrix_nominator = np.multiply(syst_matrix, current_density_vec.T)            
            multinomial_prob_matrix = np.divide(multinomial_prob_matrix_nominator, multinomial_prob_matrix_denominator)
            
            
            multinomial_prob_matrix_last_column = np.ones((nshift*nphi, 1)) - multinomial_prob_matrix.sum(axis=1, 
                                                                                                  keepdims=True)
            multinomial_prob_matrix = np.append(multinomial_prob_matrix, multinomial_prob_matrix_last_column, axis=1)
    
        # generate backprojected data n_ij, scattered photons
            for i in range(nshift*nphi):
                backprojection_data[i, :] = np.random.multinomial(sinogram_counts_vec[i], multinomial_prob_matrix[i, :])[:-1] 
        # end of Step 1
        
        # Step 2 : generate intensities using backprojection data : lambda ~ p(lambda | n_ij, y) (Gamma)
            array_shape = prior_shape + backprojection_data.sum(axis=0, keepdims=True)[0, :]
            array_scale = 1. / (prior_rate + sensitivity_array*observation_time)
        # generate random image from posterior 
            current_density_vec = np.random.gamma(array_shape, array_scale).T
        # end of Step 2
    
        # save the sample
            if (iteration > burn_in_size-1):
                output_array[:, :, iteration-burn_in_size] = np.reshape(current_density_vec, (npixels, npixels))
    # end of iteration loop 
                
        # clean memory in a loop (not very efficient, better set variables before loop)
            del(backprojection_data)
            del(multinomial_prob_matrix)
            del(multinomial_prob_matrix_nominator)
            del(multinomial_prob_matrix_denominator)
            del(multinomial_prob_matrix_last_column)
            del(array_shape)
            del(array_scale)
                
    return output_array


