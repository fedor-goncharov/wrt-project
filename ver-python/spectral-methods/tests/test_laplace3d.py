import numpy as np
import matplotlib.pyplot as plt
from spectral_laplace3d import spectral_laplace3d


# test-function - paraboloid in centered disk of radius 0.5
# outside of disk function is set to 0

# therefore spectral laplacian should be equals 6.0 in disk
# and 0 outside (with a singularity at the border)

ngrid = 128
lin = np.linspace(-1,1, ngrid)
XX, YY, ZZ = np.meshgrid(lin, lin, lin)
RR = XX*XX + YY*YY + ZZ*ZZ
#RR[RR >  (0.25)] = 0

laplace = spectral_laplace3d(RR, ngrid, 2*np.pi)
#laplace[np.fabs(laplace) > 6] = 6

# plot result of the laplacian
plt.imshow(laplace[5:100, 5:100, 64])
plt.colorbar()
#plt.axis('off')
#plt.clim(0, 1)
plt.show()
plt.savefig('test-results/test_laplace3d.png')
plt.clf()
