#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 17:26:28 2020

@author: fedor.goncharov.ol@gmail.com
"""


"""
    Run the script from the folder of the project 
    (in order to have successfull imports)
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "utilities")

from radon_transform_matrix import radon_transform2d_xray_matrix
from sinogram_noise_generator import generate_noise_emission
from gibbs_sampler_emission import gibbs_sampler_emission_gamma

# make an image - spherical layer with radiuses r_out = 0.5, r_in = 0.25
lin = np.linspace(-1., 1., 64)
[XX, YY] = np.meshgrid(lin, lin)
RR = np.sqrt(XX**2 + YY**2)
image = np.zeros((64,64))
image[RR < 0.5] = 1.
image[RR < 0.25] = 0.

# compute matrix for the Radon transform (this make take a while)
rt_system_matrix = radon_transform2d_xray_matrix(64, 64, 64, 1.0)

# compute denoised sinogram and add poisson noise
sinogram_vector = rt_system_matrix.dot(np.reshape(image, (64*64, 1)))
noise_sinogram_vector = generate_noise_emission(sinogram_vector, 20000, 3)
noise_sinogram = np.reshape(noise_sinogram_vector, (64, 64))

avg_scattered = 3*np.ones((64,64))
gamma_prior_params = (1.,1.)
init_point = np.ones((64,64))
burn_in_size = 1000
sample_size = 200

posterior_sample = gibbs_sampler_emission_gamma(noise_sinogram, rt_system_matrix, avg_scattered, gamma_prior_params, init_point, burn_in_size, sample_size)

fig1 = plt.figure()
plt.imshow(posterior_sample.sum(axis=2)/sample_size)

