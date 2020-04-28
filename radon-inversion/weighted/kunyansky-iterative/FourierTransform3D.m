% Computes direct Fourier transform of function 'u', which is 
% given on the grid in the unit cube [-rsupp, rsupp]^3; 

% ngrid : the number of equispaced points (with endpoints) in [-rsupp, rsupp]
% rsupp : radius of the domain

function Fu = FourierTransform3D(u, ngrid, rsupp)
  
  dx = 2 * rsupp / (ngrid - 1);
  nyq_freq = 1 / dx;
  dfreq = nyq_freq / ngrid;
  
  shift_frequencies_dim1 = zeros(ngrid, ngrid, ngrid);  
  shift_frequencies_dim2 = zeros(ngrid, ngrid, ngrid);
  shift_frequencies_dim3 = zeros(ngrid, ngrid, ngrid);
  
  for i = 1 : ngrid
    shift_frequencies_dim1(:, i, :) = (i-1) * dx;
    shift_frequencies_dim2(i, :, :) = (i-1) * dx;
    shift_frequencies_dim3(:, :, i) = (i-1) * dx;
  endfor
  shift_argument = shift_frequencies_dim1 .+ shift_frequencies_dim2 ...
                 .+ shift_frequencies_dim3;
  shift_operator = exp(2 * pi * 1i * rsupp *  shift_argument);
  
  Fu = (dx^3) * fftn(u);   % fft on the unit cube
  Fu = fftshift(Fu .* shift_operator);
  
endfunction
