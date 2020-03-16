import numpy as np

"""
  laplace2dsqrt.py

 Function takes on input a two-dimensional image and returns its
 square root of the laplacian in 2D.

 INPUT PARAMETERS :
     u        : matrix of size ngrid x ngrid
     ngrid    : size of one side of the matrix
     period   : length of the period of an image (geometrical length of one side)

 RETURN VALUE :
     Two-dimensional matrix of size ngrid x ngrid

 COMMENT: related mathematical theory can be found at
 Notes on FFT-based differentiation, Steven G. Johnson, MIT Applied Mathematics
"""

def spectral_laplace2dsqrt(u, ngrid, period=2.0):

  Y = np.fft.fftn(u)
  freq_factor = np.zeros((1, ngrid))
  # spectral derivative of order sqrt(laplacian)

  for k in range(ngrid):
    if (k <= (ngrid / 2)):
      freq_factor[0][k] = k
    else:
      freq_factor[0][k] = k-ngrid

  freq_factor_dim1, freq_factor_dim2 = np.meshgrid(freq_factor, freq_factor);
  U = (2*np.pi / period)*np.sqrt(freq_factor_dim1 * freq_factor_dim1 + freq_factor_dim2 * freq_factor_dim2)

  return np.fft.ifftn(U * Y)
