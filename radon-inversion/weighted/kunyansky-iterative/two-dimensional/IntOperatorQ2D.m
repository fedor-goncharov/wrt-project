function Qu = IntOperatorQ2D(ngrid, degree, weights_normalized, harmonics, cutoff, u)

% ngrid : number of points in the domain [-rsupp, rsupp]^3
% degree: degree of iterative inversion (degree >= 1)
% weights_normalized : matrix of normalized weights (ngrid x ngrid x ngrid x (degree+1)^2 - 1)
% harmonics : matrix of normalized harmonics (ngrid x ngrid x ngrid x (degree + 1)^2 - 1)
% cutoff : cuttoff function of the unit ball (ngrid x ngrid x ngrid)
% u : argument function (ngrid x ngrid x ngrid)
   
   U = cutoff .* u; % domain cutoff 
   arg = zeros(ngrid, ngrid, 2*degree); % argument vector
  
  % mutipliciation by normalized weights
   arg = weights_normalized .* U; 
  
  % direct Fourier transforms
  for i = 1 : 2*degree
    arg(:, :, i) = fft2(arg(:, :, i));
    
    harmonic = harmonics(:, :, i);
    %harmonic = flipdim(harmonic, 3);
    harmonic = permute(harmonic, [2 1]);
    
    arg(:, :, i) = arg(:, :, i) .* fftshift(harmonic); % harmonics are permutated
  endfor
  
  % summation of terms
   arg = sum(arg, 3);   % summation of all coefficients with weights in Fourier space
  
  % inverse Fourier transform
   arg = real(ifft2(arg));
  
  % domain cutoff  
   arg = cutoff .* arg;
   
  % return output 
   Qu = arg;
   
endfunction