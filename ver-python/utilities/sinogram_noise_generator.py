#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 21:39:18 2020

@author: fedor.goncharov.ol@gmail.com
"""

"""
 generators of noise for different tomographical modailities

 X-ray transmission tomography -- transmitted photons (Poisson) + scattered (Poisson)
 PET/SPECT emission tomography -- emitted photons (Poisson) + scattered (Poisson)

"""

import numpy as np


def generate_noise_from_file_xray_transmission(input_filename, output_filename, avg_intensity=1e5, T=1.0, sc_intensity=1e1):
    
    """
    function (1) reads ray/Radon transforms - (2) takes exponentiation exp(-Ra) -
    (3) generates Poisson rdv Po(N0 * t * exp(-Ra) + r * t), where N0 is given in the arguments

     Arguments : 

          input_filename : name of the input file
          output_filename : output file for noisy transforms
          avg_intensity : average intensity of the source (photons / sec)
          T : irradiation time in seconds 
          sc_intensity : intensity of scattered photons (scattered photons / sec)

    Output : 
          binary file containing array with values of type np.int64
  
    Example : 
          generate_nosie_xray_transmission('ray_transforms.bin', 'noise_ray_transforms.bin', avg_intensity=1e5, sc=1e1)
    """    
    # read data
    dtype = np.dtype(np.double)
    data_array = np.fromfile(input_filename, dtype, count=-1, sep="")

    # generate noise
    output_array = np.random.poisson(avg_intensity * T * np.exp((-1.0)*data_array) + sc_intensity*T)
    
    # write to file
    fid = open(output_filename, mode='wb')
    output_array.tofile(fid, sep="")
    fid.close()
    
    
def generate_noise_xray_transmission(sinogram, avg_intensity=1e5, T=1.0, sc_intensity=1e1):
    
    """
    function 
    (1) takes on input ray/Radon transforms (sinogram) 
    (2) takes exponentiation exp(-Ra) 
    (3) generates Poisson rdv Po(N0 * t * exp(-Ra) + r * t), where N0 is given in the arguments
     Arguments : 

          sinogram      : matrix of ray transforms on the plane of type (nshift, nphi)
          avg_intensity : average intensity of the source (photons / sec)
          T : irradiation time in seconds 
          sc_intensity  : intensity of scattered photons (scattered photons / sec)

    Output : 
          matrix of type as sinogram with Poisson noise data 
  
    Example : 
          noisy_sinogram = generate_nosie_xray_transmission(sinogram, avg_intensity=1e5, sc=1e1)
    """
    return np.random.poisson(avg_intensity * T * np.exp(-1.0 * sinogram) + sc_intensity*T)
    
def generate_noise_from_file_emission(input_filename, output_filename, N_max, N_sc):
    """
    function 
    (1) reads weighted ray/Radon transforms 
    (2) finds maximal value P_wf and normalizes data P_wf / max(P_wf) -
    (3) generates Poisson rdvs Po(N_max * P_wf / max(P_wf) + r), where N_max is given in the arguments

    Arguments : 
    
        input_filename  : file with ray transforms
        output_filename : output file for noisy transforms
        N_max : average maximal number of photons per ray
        N_sc : average number of scattered photons

    Output : 
        binary file containing array with values of type np.int64

    Example : 
        generate_noise_from_file_emission('ray_transforms.bin', 'noise_ray_transforms.bin', 1e2, 3)
    """
    
    # read data
    dtype = np.dtype(np.double)
    data_array = np.fromfile(input_filename, dtype, count=-1, sep="")

    # generate noise 
    output_array = np.random.poisson(N_max*data_array / np.max(data_array) + N_sc)

    # write to file
    fid = open(output_filename, mode="wb")
    output_array.tofile(fid, sep="")
    fid.close()

def generate_noise_emission(sinogram, N_max, N_sc):
    
    """
    function 
    (1) takes weighted ray/Radon transforms as input matrix  
    (2) finds maximal value P_wf and normalizes data P_wf / max(P_wf) 
    (3) generates Poisson rdvs Po(N_max * P_wf / max(P_wf) + r), where N_max is given in the arguments

    Arguments : 
    
        sinogram : matrix of type (nshift, nphi) with Radon transfroms on the plane
        N_max : average maximal number of photons per ray
        N_sc : average number of scattered photons

    Output : 
        matrix of the same type as sinogram with Radon transforms with Poisson noise

    Example : 
       noisy_sinogram=generate_noise_emission(sinogram, 1e2, 3)
    """
    return np.random.poisson(N_max*sinogram / np.max(sinogram[:]) + N_sc)
