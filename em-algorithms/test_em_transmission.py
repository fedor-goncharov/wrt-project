#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 13:52:27 2020

@author: fedor.goncharov.ol@gmail.com
"""


import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "../ver-python/utilities")

from radon_transform_matrix import radon_transform2d_xray_matrix
from sinogram_noise_generator import generate_noise_xray_transmission
from em_transmission import em_transmission_convex_nr1, em_transmission_convex_nr2

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
ray_transforms_vector = rt_system_matrix.dot(np.reshape(image, (64*64, 1)))
noise_ray_transforms_vector = generate_noise_xray_transmission(ray_transforms_vector, avg_intensity=1e3, T=1.0, 
                                                         sc_intensity=1e1)
noise_ray_transforms = np.reshape(noise_ray_transforms_vector, (64, 64))



# run EM-algorithm
avg_scattered = 1e1*np.ones((64,64))
max_iterations = 100
relative_err_level = 1e-3
init_point = np.ones((64,64))

reconstruction_em_nr1 = em_transmission_convex_nr1(noise_ray_transforms, rt_system_matrix, 1e3*np.ones((64,64)),
                                               avg_scattered, max_iterations, 
                                               relative_err_level, init_point)
# nr1 - algorithm has a tendency to be numerically unstable when estimating attenuation values near zero
# in this example iterations from 1 to 8 give reasonable images, then the process completely diverges

fig1 = plt.figure()
plt.imshow(reconstruction_em_nr1)

# testing EM_emission_algorithm_mlem3

reconstruction_em_nr2 = em_transmission_convex_nr2(noise_ray_transforms, rt_system_matrix, 1e3*np.ones((64,64)), 
                                               avg_scattered, max_iterations, 
                                               relative_err_level, init_point)
fig2 = plt.figure()
plt.imshow(reconstruction_em_nr2)

# plt.close(fig1)
# plt.close(fig2)
