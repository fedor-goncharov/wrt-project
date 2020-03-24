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
system_matrix = radon_transform2d_xray_matrix(64, 64, 64, 1.0)

sinogram_vector = system_matrix.dot(np.reshape(image, (64*64, 1)))
noise_sinogram_vector = generate_noise_emission(sinogram_vector, 200, 10)
noise_sinogram = np.reshape(noise_sinogram_vector, (64, 64))

avg_scattered = 10*np.ones((64,64))
gamma_prior_params = (1.,1.)
init_point = np.ones((64,64))
burn_in_size = 0
sample_size = 100

posterior_sample = gibbs_sampler_emission_gamma(noise_sinogram, system_matrix, avg_scattered, gamma_prior_params, init_point, burn_in_size, sample_size)

fig1 = plt.figure()
plt.imshow(np.mean(posterior_sample, axis=2))
plt.title("Posterior mean")
plt.colorbar()

fig2 = plt.figure()
plt.imshow(np.std(posterior_sample, axis=2))
plt.title("Posterior standard deviation")
plt.colorbar()

