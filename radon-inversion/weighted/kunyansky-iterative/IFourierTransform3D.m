% Computes inverse Fourier transform of function 'u' 
% given in the unit cube [-rsupp, rsupp]^3; where ngrid is the number of points
% on the unit dimension of the cube including endpoints.

% works for odd grid size, but not for even

function iFu = IFourierTransform3D(u, ngrid, rsupp)
  
  dx = 2 * rsupp / (ngrid - 1);
  nyq_freq = 1 / dx;
  dfreq = nyq_freq / ngrid;
  
  shift_frequencies_dim1 = zeros(ngrid, ngrid, ngrid);  
  shift_frequencies_dim2 = zeros(ngrid, ngrid, ngrid);
  shift_frequencies_dim3 = zeros(ngrid, ngrid, ngrid);
  
  for i = 1 : ngrid
    shift_frequencies_dim1(:, i, :) = (i - 1) * dx;
    shift_frequencies_dim2(i, :, :) = (i - 1) * dx;
    shift_frequencies_dim3(:, :, i) = (i - 1) * dx;
  endfor
  
  shift_argument = shift_frequencies_dim1 .+ shift_frequencies_dim2 ...
                 .+ shift_frequencies_dim3;
  shift_operator = exp((-2)*pi * 1i * rsupp * shift_argument);
  
  iFu = (dx^3) * (ngrid)^3 * ifftn(u);   % ifft on the unit cube (normalization, because IDFT ~ 1/N sum (exp)
  iFu = fftshift(iFu .* shift_operator); % correcting phase shift and centralizing frequency 
  
endfunction