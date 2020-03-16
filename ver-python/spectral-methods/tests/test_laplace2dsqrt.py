import numpy as np
import matplotlib.pyplot as plt
from hfilter2d import hfilter2d

import sys
sys.path.append('../radon-adjoint/') # add path to adjoint transform
from adjradon2d import adjradon2d


filename = '~/Documents/programming/tomography/wrt-project/reduction_ray_radon/binary/thesis_data/ph1_att_weak_ray_slice.csv'
ngrid = 513
nphi = 128
nshift = 129

adjoint = adjradon2d(filename, nphi, nshift, 1.0, 1.0, ngrid) / (2*np.pi)
hadjoint = np.real(hfilter2d(adjoint, 513, 2.0))/2


plt.imshow(hadjoint)
plt.clim(0,1)
plt.axis('off')
plt.savefig('test_hfilter2d.png')
plt.clf()

