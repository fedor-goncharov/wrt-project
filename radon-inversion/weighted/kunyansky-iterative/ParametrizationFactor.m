% Correction of weights due to change of parametrization 
% in integral and in spherical harmonics decomposition.

% In spect in slice-by-slice approach the new weight W depends 
% on w in the following form : W(x, \theta, phi) = w(x, phi + pi/2), 
% so the factor of (pi/2) should be respected when doing harmonics decompositions.
function factor = ParametrizationFactor(k,m)
  
  factor = exp(i*m * (pi/2));
  
endfunction