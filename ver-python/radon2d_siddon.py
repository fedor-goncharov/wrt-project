import numba as nb
import numpy as np

@nb.njit(parallel=True)
def radon2d_sidon(image, ntheta, nshift, radius=1.0):
    
    # image size
    npixels = image.shape[0]
    dx = 2.0*radius / npixels
    shifts = np.linspace(-radius + dx/2, radius-dx/2, nshift)
    theta = np.linspace(0, 2*np.pi-2*np.pi/ntheta, ntheta)
    
    proj = np.zeros((ntheta, nshift))
    for i_theta in range(ntheta):
        for i_shift in range(nshift):
            proj[i_theta][i_shift] = siddon_line_projector(image, theta[i_theta], shifts[i_shift], radius)    
    return proj
