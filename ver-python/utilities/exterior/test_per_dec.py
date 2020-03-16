import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../radon-adjoint/') # add path to adjoint transform
from adjradon2d import adjradon2d
from per_dec2d import periodic_smooth_decomp2d
from adjradon3d import adjradon3d
from per_dec3d import periodic_smooth_decomp3d

filename = '~/Documents/programming/tomography/wrt-project/reduction_ray_radon/binary/thesis_data/ph1_att_weak_ray_slice.csv'
ngrid = 513
nphi = 128
nshift = 129

adjoint = adjradon2d(filename, nphi, nshift, 1.0, 1.0, ngrid)
periodic_adjoint, smooth_adjoint = periodic_smooth_decomp2d(adjoint)

plt.imshow(periodic_adjoint)
plt.axis('off')
plt.clim(0,1)
plt.savefig('test_periodic2d.png')
plt.clf()

plt.imshow(smooth_adjoint)
plt.axis('off')
plt.clim(0,1)
plt.savefig('test_smooth2d.png')
plt.clf()


filename = '~/Documents/programming/tomography/wrt-project/reduction_ray_radon/binary/thesis_data/ph1_att_weak_radon.csv'
ngrid = 33
ntheta = 128
nphi = 128
nshift = 129

adjoint = adjradon3d(filename, ntheta, nphi, nshift, 1.0, 1.0, ngrid)
periodic_adjoint, smooth_adjoint = periodic_smooth_decomp3d(adjoint)

plt.imshow(periodic_adjoint[:, :, 16])
plt.axis('off')
plt.clim(0,1)
plt.savefig('test_periodic3d.png')
plt.clf()

plt.imshow(smooth_adjoint[:, :, 16])
plt.axis('off')
plt.clim(0,1)
plt.savefig('test_smooth3d.png')
plt.clf()