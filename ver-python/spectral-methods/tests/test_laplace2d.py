import numpy as np
import matplotlib.pyplot as plt
from laplace2d import laplace2d

ngrid = 513

# create a disk in the center 
# laplace transform should fine edges
lin = np.linspace(-1,1, ngrid)
XX, YY = np.meshgrid(lin, lin)
RR = np.sqrt(XX*XX + YY*YY)
RR[RR >  0.5] = 0
RR[RR > 0] = 1

laplace_adjoint = laplace2d(RR, ngrid, 1000.0)

plt.imshow(laplace_adjoint)
plt.axis('off')
plt.colorbar()
plt.clim(0, 1)
plt.savefig('test_laplace2d.png')
plt.clf()