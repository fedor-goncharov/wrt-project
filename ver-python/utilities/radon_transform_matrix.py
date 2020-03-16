#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 21:39:18 2020

@author: fedor.goncharov.ol@gmail.com
"""

"""
  Generators for system matrices for various Radon-type transforms 
  on the plane.
"""

import numpy as np

def radon_transform2d_classical(image, phi, shift, radius):
    """
	  Classical Radon transform along the line on the plane for parallel beam geometry. 
	  
	  Arguments:
	  	  image : matrix of type (npixels, npixels) (only square images are acceptable)
          phi : number of projections (angles)
          shift -- number of shifts per projection
          radius -- radius of the support       
          
      Output: 
          value of the Radon transform R(image)(shift, phi)
	"""
    # if line is outside the radius
    if (shift > radius):
        return 0.
    
    npixels = image.shape[0]
    
    center_point = np.array([[shift * np.cos(phi)], [shift * np.sin(phi)]])
    dir_vector = np.array([[np.cos(phi + np.pi/2.0)], [np.sin(phi + np.pi/2.0)]])
    
    
    line_points = center_point +  np.multiply(dir_vector, np.linspace(-1.5*radius, 1.5*radius, 4*npixels))
    step = 3.0*radius / (4*npixels - 1)
    
    delta = 2.0*radius / npixels
    line_pixels = np.floor((line_points + radius)/delta).astype(int)
    line_pixels = line_pixels[:, ((line_pixels[0] > -1)*(line_pixels[0] < npixels))
                              *((line_pixels[1] > -1)*(line_pixels[1] < npixels))]
    
    
    values = image[line_pixels[1, :], line_pixels[0, :]]
    return np.sum(values) * step
    
    #integral = 0.0
    #for point in line_pixels.T:
        
    #    # find pixel corresponding to point
    #    if ((np.min(point) > -1) and (np.max(point) < npixels)):
    #        integral += image[point[1]][point[0]]
    #   
    #return integral * step

def radon_transform2d_xray(image, nproj, nshift, radius, max_proj_angle=np.pi):
    """
      Compute the Radon transform in X-ray on the plane for an image for 
      parallel beam geometry. 
      
      Arguments: 
          npixels -- number of pixels per dimension
          nproj -- number of projections (angles)
          nshift -- number of shifts per projection
          radius -- radius of the support
          max_angle -- maximal angle of sampling (by default = pi)
          
      Output: 
          matrix of size (nproj*nshift, npixels^2)
          
    """
    npixels = image.shape[0]
    delta = 2.0*radius / npixels
    shifts = np.linspace(-radius + delta/2, radius-delta/2, nshift)
    projections = np.linspace(0., max_proj_angle, nproj, endpoint = False)
    
    sinogram = np.zeros((nshift, nproj))
    for i_shift in range(nshift):
        for i_proj in range(nproj):
            sinogram[i_shift][i_proj] = radon_transform2d_classical(image, projections[i_proj], 
                                                          shifts[i_shift], radius)
    
    return sinogram


def radon_transform2d_xray_matrix(npixels, nproj, nshift, radius, max_proj_angle=np.pi):
    """
      Compute matrix for the Radon transform in X-ray on the plane for 
      parallel beam geometry. 
      
      Arguments: 
          npixels -- number of pixels per dimension
          nproj -- number of projections (angles)
          nshift -- number of shifts per projection
          radius -- radius of the support
          max_angle -- maximal angle of sampling (by default = pi)
          
      Output: 
          matrix of size (nproj*nshift, npixels^2)
          
    """
    delta = 2.0*radius / npixels
    shifts = np.linspace(-radius + delta/2., radius-delta/2., nshift)
    proj_angles = np.linspace(0., max_proj_angle, nproj, endpoint = False)  
    # create output matrix
    matrix = np.zeros((nshift*nproj, npixels**2))
    
    for i_pixel in range(npixels**2):
        i_pixel_x = i_pixel % npixels
        i_pixel_y = i_pixel // npixels
        pixel_x = -radius + delta/2 + i_pixel_x * delta
        pixel_y = -radius + delta/2 + i_pixel_y * delta
        
        # create image with fixed pixel
        image = np.zeros((npixels, npixels))
        image[i_pixel_y, i_pixel_x] = 1.0
        
        #create sinogram
        sinogram = np.zeros((nshift, nproj))
        
        # for each pixel compute sinogram
        for i_shift in range(nshift):
            for i_phi in range(nproj):
                
                shift = shifts[i_shift]
                phi = proj_angles[i_phi]
                if (np.abs((np.dot([pixel_x, pixel_y], [np.cos(phi), np.sin(phi)])) - shift) > delta*np.sqrt(2)):
                    sinogram[i_shift, i_phi] = 0.
                else : 
                    sinogram[i_shift, i_phi] = radon_transform2d_classical(image, phi, shift, radius)
        
        matrix[:, i_pixel] = np.reshape(sinogram, (1, nshift*nproj))
        
    return matrix

def radon_transform2d_matrix_pet():
    """
      Compute matrix for the Radon transform in PET on the plane for 
      parallel beam geometry. 
      
      Arguments: 
          npixels -- number of pixels per dimension
          att_map -- matrix of type (npixels, npixels) for the attenuation map
          nproj -- number of projections (angles)
          nshift -- number of shifts per projection
          radius -- radius of the support
          max_angle -- maximal angle of sampling (by default = 2pi)
          
      Output: 
          matrix of size (nproj*nshift, npixels^2)
     
      Remark: 
          for attenuation map it is assumed that it is supported 
          in the centered ball of given radius
          
    """
    return None

def radon_transform2d_spect(): 
    """
      Compute matrix for the Radon transform on the plane for 
      parallel beam geometry. 
      
      Arguments: 
          npixels -- number of pixels per dimension
          att_map -- matrix of type (npixels, npixels) for the attenuation map
          nproj -- number of projections (angles)
          nshift -- number of shifts per projection
          radius -- radius of the support
          max_angle -- maximal angle of sampling (by default = 2pi)
          
      Output: 
          matrix of size (nproj*nshift, npixels^2)
          
    """
    return None