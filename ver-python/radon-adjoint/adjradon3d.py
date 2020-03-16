import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

"""
 adjradon3d.py

 Function takes on input 2-plane Radon transforms in 3D of a function and
 computes its adjoint Radon transform. The Radon data are assumed to be given on
 a grid with uniform shifts and uniform-gauss directions on the sphere.

 INPUT PARAMETERS :
     filename : file with Radon transforms
     ngrid    : number of pixels in one direction in the XYZ-domain
     nphi     : number of equatorial angles
     ntheta   : number of latitude angles
     nshift   : number of shifts along each projection
     rsupp    : radius of the support of the test-function (by default = 1.0)

     size     : adjoint is computed in cube of the form [-size, size]^3
     ngrid    : number of points in the grid of the cube in one dimension


RETURN VALUE :
     A three-dimensional matrix is returned of size ngrid x ngrid x ngrid, where
     each element stands for the adjoint  transform at the corresponding points
     of the grid.
"""

def adjradon3d(filename, nphi, ntheta, nshift, rsupp=1.0, size=1.0, ngrid=129) :
  # read (weighted) Radon transforms from file
  radon_data = pd.read_csv(filename, header=None, sep="\n")
  radon_values = radon_data.values[:]
  radon_values = np.reshape(radon_values, (ntheta, nphi, nshift), order="F")

  # set standard variables
  dphi   = 2 * np.pi / nphi;
  dshift = 2 * rsupp / (nshift-1)

  # set grid in Radon space
  # get equatorial angles 'phi' on [0, 2*pi)
  angles_phi = np.arange(nphi) * dphi

  # get latitude angles by the Gauss-Quadrature rule; 'theta' on (0, pi)
  nodes_gauss, weights_gauss =  np.polynomial.legendre.leggauss(ntheta)
  angles_theta = np.arccos(nodes_gauss);
  shifts = np.linspace(-rsupp, rsupp, nshift)

  # create grid XYZ
  lin = np.linspace(-size, size, ngrid)
  XX, YY, ZZ = np.meshgrid(lin, lin, lin)
  XX = np.reshape(XX, (1, ngrid**3), order="F")
  YY = np.reshape(YY, (1, ngrid**3), order="F")
  ZZ = np.reshape(ZZ, (1, ngrid**3), order="F")

  adjoint = np.zeros((1, ngrid**3))

  # points in the slice of size 3 x ngrid^2
  points = np.concatenate((XX, YY, ZZ))

    # integral over sphere for each z
  for i_phi in range(nphi):
    for i_theta in range(ntheta):

        phi = angles_phi[i_phi]
        theta = angles_theta[i_theta]

        # direction on the sphere S^2
        direction = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])

        # compute scalar products with direction for all points in the grid
        interpolate_shifts = np.dot(np.transpose(points), direction)

        # adjoint Radon at all points of the z-slice for direciton (phi, theta)
        add_interpolated = np.interp(interpolate_shifts, shifts, radon_values[i_theta, i_phi, :],  left=0., right=0.)
        #interpolant = interp1d(shifts, radon_values[i_theta, i_phi, :], kind='zero', bounds_error=False, fill_value='extrapolate', assume_sorted=True)
        #add_interpolated = interpolant(interpolate_shifts)

        # add
        adjoint += add_interpolated * dphi * weights_gauss[i_theta];

  # give back array in the fortran ordering
  adjoint = np.reshape(adjoint, (ngrid, ngrid, ngrid), order="F");
  return adjoint
