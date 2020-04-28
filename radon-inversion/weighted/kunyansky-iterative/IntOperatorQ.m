% Application of integral operator Q_{W,D,m} to function 'u', where 
% W is the weight of finite number of harmonics, D is a unit ball in 3D,
% m - highest degree of harmonics expansion; [Goncharov, Novikov, 2016]

function Qu = IntOperatorQ(ngrid, degree, weights_normalized, harmonics, cutoff, u)

% ngrid : number of points in the domain [-rsupp, rsupp]^3
% degree: degree of iterative inversion (degree >= 1)
% weights_normalized : matrix of normalized weights (ngrid x ngrid x ngrid x (degree+1)^2 - 1)
% harmonics : matrix of normalized harmonics (ngrid x ngrid x ngrid x (degree + 1)^2 - 1)
% cutoff : cuttoff function of the unit ball (ngrid x ngrid x ngrid)
% u : argument function (ngrid x ngrid x ngrid)
   
   U = cutoff .* u; % domain cutoff 
   arg = zeros(ngrid, ngrid, ngrid, (degree + 1)*(2*degree + 1) - 1); % argument vector
  
  % mutipliciation by normalized weights
   arg = weights_normalized .* U; 
  
  % direct Fourier transforms
  for i = 1 : ((degree + 1)*(2*degree + 1) - 1)
    arg(:, :, :, i) = fftn(arg(:, :, :, i));
    
    harmonic = harmonics(:, :, :, i);
    harmonic = flipdim(harmonic, 3);
    harmonic = permute(harmonic, [2 1 3]);
    
    arg(:, :, :, i) = arg(:, :, :, i) .* fftshift(harmonic); % harmonics are permutated
  endfor
  
  % summation of terms
   arg = sum(arg, 4);   % summation of all coefficients with weights in Fourier space
  
  % inverse Fourier transform
   arg = real(ifftn(arg));
  
  % domain cutoff  
   arg = cutoff .* arg;
   
  % return output 
   Qu = arg;
   
endfunction