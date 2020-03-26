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
noise_sinogram_vector = generate_noise_emission(sinogram_vector, 20000000, 10)
noise_sinogram = np.reshape(noise_sinogram_vector, (64, 64))

avg_scattered = 10*np.ones((64,64))
gamma_prior_params = (1.,1.)
init_point = np.ones((64,64))
burn_in_size = 0
sample_size = 1000

posterior_sample = gibbs_sampler_emission_gamma(noise_sinogram, system_matrix, avg_scattered, gamma_prior_params, init_point, burn_in_size, sample_size)

fig1 = plt.figure()
plt.imshow(np.mean(posterior_sample, axis=2))
plt.title("Posterior mean")
plt.colorbar()

fig2 = plt.figure()
plt.imshow(np.std(posterior_sample, axis=2))
plt.title("Posterior standard deviation")
plt.colorbar()

# autocorrelation test

# create image with one "hot pixel"
hot_image = image 
hot_image[32,22] = 10. # create hot pixel

# for comparison we choose also a "cold" pixel
cold_projector = np.zeros((64, 64)) 
cold_projector[32,32] = 1.  # projector for "cold" pixel

hot_projector = np.zeros((64, 64))
hot_projector[32,22] = 1. # projector for "hot" pixel

hot_sinogram_vector = system_matrix.dot(np.reshape(hot_image, (64*64, 1)))
hot_noise_sinogram_vector = generate_noise_emission(hot_sinogram_vector, 2e6, 10)
hot_noise_sinogram = np.reshape(hot_noise_sinogram_vector, (64, 64))

avg_scattered = 10*np.ones((64,64))
gamma_prior_params = (1.,1.)
init_point = np.ones((64,64))
burn_in_size = 0
sample_size = 1e3

hot_posterior_sample = gibbs_sampler_emission_gamma(hot_noise_sinogram, system_matrix, avg_scattered, gamma_prior_params, init_point, burn_in_size, sample_size)

autocorrelation_array_hot = np.array([])
autocorrelation_array_cold = np.array([])

for i in range(sample_size):
    autocorrelation_array_hot = np.append(autocorrelation_array_hot, np.sum(np.multiply(hot_projector, hot_posterior_sample[:, :, i])[:]))
    autocorrelation_array_cold = np.append(autocorrelation_array_cold, np.sum(np.multiply(cold_projector, hot_posterior_sample[:, :, i])[:]))

# autocorrelation function
def acf(x, length=100):
    return np.array([1]+[np.corrcoef(x[:-i], x[i:])[0,1]  \
        for i in range(1, length)])

autocorrelation_lag_array_hot = acf(autocorrelation_array_hot)
autocorrelation_lag_array_cold = acf(autocorrelation_array_cold)

fig3 = plt.plot(autocorrelation_lag_array_hot)
plt.plot(autocorrelation_lag_array_cold)
plt.legend(["Hot pixel", "Cold pixel"])
plt.xlabel("lag")
plt.ylabel("correlation")
plt.title("autocorrelation function")


