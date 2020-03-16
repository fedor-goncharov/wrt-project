import numpy as np
import pandas as pd

"""
 adjradon2d.py

 Function takes on input Radon transforms of a function in 2D and computes its
 adjoint Radon transform. The Radon data are assumed to be given on the
 (well-known) parallel beam grid.

 INPUT PARAMETERS :
     filename : file with Radon transforms
     nphi     : number projections
     nshift   : number of shifts per projection
     rsupp    : radius of the support of a test-function (by default = 1.0)

     size     : adjont is computed in a cube of the form [-size,size]^2
     (by default size=1.0)
     ngrid    : number of pixels per dimension in the cubic grid
     (by default ngrid = 129)

 RETURN VALUE :
     A two-dimensional matrix is returned of size ngrid x ngrid where elements stand for the adjoint
     transform at the corresponding geometrical points.
"""

def adjradon2d(filename, nphi, nshift, rsupp=1.0, size=1.0, ngrid=129) :
  # read (weighted) Radon transforms from file
  radon_data = pd.read_csv(filename, header=None, sep="\n")
  radon_values = radon_data.values[:]
  radon_values = np.reshape(radon_values, (nphi, nshift), order="F")

  # set standard variables
  dphi   = 2 * np.pi / nphi;
  dshift = 2 * rsupp / (nshift-1)

  # set grid in Radon space
  # angles 'phi' on [0, 2*pi)
  angles_phi = np.arange(nphi) * dphi
  # shifts on [-rsupp, rsupp]
  shifts = np.linspace(-rsupp, rsupp, nshift)

  # create grid XY
  lin = np.linspace(-size, size, ngrid)
  XX,YY = np.meshgrid(lin, lin)
  XX = np.reshape(XX, (1, ngrid**2), order="F") #reshape into long arrays
  YY = np.reshape(YY, (1, ngrid**2), order="F")


  # main computation
  adjoint = np.zeros((1, ngrid**2))

  # points in the slice of size 2 x L^2
  points = np.concatenate((XX, YY))

    # integral over sphere for each z
  for i_phi in range(nphi):

        phi = angles_phi[i_phi]

        # direction on the circle S^1
        direction = np.array([np.cos(phi), np.sin(phi)])

        # compute scalar products with direction for all points in the grid
        new_shifts = np.dot(np.transpose(points), direction)

        # values of the Radon transforms for new_shifts, fixed angle phi[i_phi]
        add_interpolated = np.interp(new_shifts, shifts, radon_values[i_phi, :],  left=0., right=0.)

        # Riemann term in the int. for adjoint Radon transform
        adjoint += add_interpolated * dphi

  # give back the matrix
  adjoint = np.reshape(adjoint, (ngrid, ngrid), order="F");
  return adjoint
