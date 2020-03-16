'''
This is a direct extension of result of Moisan, L. to 3D.
'''

'''Periodic & smooth image decomposition
References
----------
Periodic Plus Smooth Image Decomposition
Moisan, L. J Math Imaging Vis (2011) 39: 161.
doi.org/10.1007/s10851-010-0227-1
'''


import numpy as np

def periodic_smooth_decomp3d(I):
    '''Performs periodic-smooth image decomposition
    Parameters
    ----------
    I : np.ndarray
        [M, N] image. will be coerced to a float.
    Returns
    -------
    P : np.ndarray
        [M, N] image, float. periodic portion.
    S : np.ndarray
        [M, N] image, float. smooth portion.
    '''
    u = I.astype(np.float64)
    v = u2v(u)
    v_fft = np.fft.fftn(v)
    s = v2s(v_fft)
    s_f = np.real(np.fft.ifftn(s))
    p = u - s_f # u = p + s
    return p, s_f

def u2v(u):
    '''Converts the image `u` into the image `v`
    Parameters
    ----------
    u : np.ndarray
        [M, N] image
    Returns
    -------
    v : np.ndarray
        [M, N] image, zeroed expect for the outermost rows and cols
    '''
    v = np.zeros(u.shape, dtype=np.float64)

    v[0, :, :] = np.subtract(u[-1, :, :], u[0,  :, :], dtype=np.float64) # x = 0
    v[-1,:, :] = np.subtract(u[0,  :, :], u[-1, :, :], dtype=np.float64) # x = M -1

    v[:,  0, :] += np.subtract(u[:, -1, :], u[:,  0, :], dtype=np.float64) # y = 0
    v[:, -1, :] += np.subtract(u[:,  0, :], u[:, -1, :], dtype=np.float64) # y = N - 1

    v[:, :,  0] += np.subtract(u[:, :, -1], u[:, :,  0], dtype=np.float64) # z = 0
    v[:, :, -1] += np.subtract(u[:, :,  0], u[:, :, -1], dtype=np.float64) # z = K - 1
    return v

def v2s(v_hat):
    '''Computes the maximally smooth component of `u`, `s` from `v`
    s[q, r] = v[q, r] / (2*np.cos( (2*np.pi*q)/M )
        + 2*np.cos( (2*np.pi*r)/N ) - 4)
    Parameters
    ----------
    v_hat : np.ndarray
        [M, N] DFT of v
    '''
    M, N, K = v_hat.shape

    q = np.arange(M).reshape(M, 1, 1).astype(v_hat.dtype)
    r = np.arange(N).reshape(1, N, 1).astype(v_hat.dtype)
    p = np.arange(K).reshape(1, 1, K).astype(v_hat.dtype)

    den = (2*np.cos( np.divide((2*np.pi*q), M) ) \
         + 2*np.cos( np.divide((2*np.pi*r), N) ) \
         + 2*np.cos( np.divide((2*np.pi*p), K) )
         - 6)
    s = np.divide(v_hat, den, out=np.zeros_like(v_hat), where=den!=0)
    s[0, 0, 0] = 0
    return s
