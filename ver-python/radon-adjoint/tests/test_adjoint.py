import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from adjradon2d import adjradon2d
from adjradon3d import adjradon3d

"""
    Various test of python utilities for Radon transforms.
"""


# test adjoint3d.py
#filename = '~/Documents/programming/tomography/wrt-project/reduction_ray_radon/binary/ph1_att_weak_radon.csv'
#ngrid = 33
#nphi = 128
#ntheta = 128
#nshift = 129

#adjoint = adjradon3d(filename, nphi, ntheta, nshift, 1.0, 1.0, ngrid)

#plt.imshow(adjoint[:, :, 16])
#plt.axis('off')
#plt.clim(0,1)
#plt.savefig('test_adjoint3d.png')
#plt.clf()


# test adjoint2d.py
filename = '~/wrt-link/reduction-algorithm/binary/thesis_data/ph2_att_weak_ray_slice.csv'
ngrid = 129
nphi = 128
nshift = 129
adjoint = adjradon2d(filename, nphi, nshift, 1.0, 1.0, 512)

radon_data = pd.read_csv(filename, header=None, sep="\n")
radon_values = radon_data.values[:]
radon_values = np.reshape(radon_values, (nphi, nshift), order="F")


plt.imshow(adjoint)
plt.axis('off')
plt.savefig('test-results/test-adjradon2d.png')
plt.show()
plt.clf()


# test adjoint3d.py
filename = 'test-source/test-adjradon3d.csv'
ngrid = 129
nphi = 128
ntheta = 128
nshift = 129
adjoint = adjradon3d(filename, nphi, ntheta, nshift, 1.0, 2.0, 32)

plt.imshow(adjoint[:, :, 16])
plt.axis('off')
plt.savefig('test-results/test-adjradon3d.png')
plt.show()
plt.clf()
