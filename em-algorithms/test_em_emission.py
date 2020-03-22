#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 11:56:15 2020

@author: fedor.goncharov.ol@gmail.com
"""

"""
    Run the script from the folder of the project 
    (in order to have successfull imports)
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "../ver-python/utilities")

from radon_transform_matrix import radon_transform2d_xray_matrix
from sinogram_noise_generator import generate_noise_emission
from em_emission import EM_algorithm_emission2d, EM_algorithm_emsission2d_mlem3

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
noise_sinogram_vector = generate_noise_emission(sinogram_vector, 20, 3)
noise_sinogram = np.reshape(noise_sinogram_vector, (64, 64))

# run EM-algorithm
avg_scattered = 3*np.ones((64,64))
max_iterations = 400
relative_err_level = 1e-3
init_point = np.ones((64,64))

reconstruction_em = EM_algorithm_emission2d(noise_sinogram, rt_system_matrix, avg_scattered, max_iterations, 
                                         relative_err_level, init_point)
fig1 = plt.figure()
plt.imshow(reconstruction_em)

# testing EM_emission_algorithm_mlem3

m_coeffs = np.ones((64, 64))
reconstruction_em3 = EM_algorithm_emsission2d_mlem3(noise_sinogram, rt_system_matrix, avg_scattered, m_coeffs, 
                                               max_iterations, relative_err_level, init_point)
fig2 = plt.figure()
plt.imshow(reconstruction_em3)

# plt.close(fig1)
# plt.close(fig2)
