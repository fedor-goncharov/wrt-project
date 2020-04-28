% chang3d.m

function reconstruction = chang3d(ngrid, w0, input_filename, ntheta, nphi, nshift, rsupp=1.0, padd=4.0)

% depends on rtft3d.m, nfft_reconstruct_3d.m
 
% Function takes on input values for weighted Radon transforms along 
% planes in 3D from file (input_filename) and applies classical 3D Chang-type formula.
% The formula containts inversion of Radon transforms and here 
% it is implemented via the Fourier Slice Projection Theorem. 

% NOTE: for usage of Chang-type formula the values for w0 (zero order harmonic 
% of the weight) must be provided.

%% Usage of the script:

%   ngrid             : size of the output image in pixels (ngrid x ngrid x ngrid) : even integer
%   w0                : zero harmonic of the weight : float matrix of size (ngrid x ngrid x ngrid)
%   input_filename    : file with Radon transforms (order [shift, phi, theta]) : string 
%   ntheta            : number of latitude angles in (0, pi) : integer 
%   nphi              : number of projections in azimuth angle (0, 2pi) : integer
%   nshift            : number of shifts in [-1,1] per one direction : integer
%   rsupp             : radius of the support of the test-function (by default=1.0) : float
%   padd              : factor for zero padding each projection for FFT's (by default=4) : float

% OUTPUT: voxel image of size ngrid x ngrid x ngrid
 
  [nodes, values, jweights] = rtft3d(input_filename, ntheta, nphi, nshift, rsupp, padd);
  reconstruction = nfft_reconstruct_3d(ngrid, nodes, values, jweights) ./ w0;
  
endfunction