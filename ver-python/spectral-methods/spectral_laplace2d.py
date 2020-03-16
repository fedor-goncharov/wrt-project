import numpy as np

"""
 laplace2d.py

 Function takes on input an image matrix and returns its spectral Laplace transform in 3D.

 INPUT PARAMETERS :
     u        : matrix of size ngrid x ngrid x ngrid
     ngrid    : number of elements of one side of the matrix
     period   : length of the period of an image (geometrical length of one side)

 RETURN VALUE :
     A two-dimensional matrix of size ngrid x ngrid

 COMMENT: related mathematical theory can be found at
 Notes on FFT-based differentiation, Steven G. Johnson, MIT Applied Mathematics

"""

def spectral_laplace2d(u, ngrid, period=2.0):

  Y = np.fft.fftn(u)
  freq_factor = np.zeros((1, ngrid))

  # spectral derivatives of order 2
  for k in range(ngrid):
    if (k < (ngrid / 2 + 1)):
      freq_factor[0][k] = k
    else:
      freq_factor[0][k] = (k-ngrid)


  u1, u2 = np.meshgrid(freq_factor, freq_factor)
  U = ((2*np.pi / period)**2)*((u1*u1 + u2*u2))

  return np.real(np.fft.ifftn(U * Y))
