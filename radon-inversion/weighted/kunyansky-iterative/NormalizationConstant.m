% Normalization constant of weights expansions when doing 
% transitions from expansions of spect weight in 2D to its expansions
% in 3D in spherical harmonics.

% k : degree of spherical harmonic (positive integer)
% m : order of spherical harmonic [-m, m]

function constant = WeightsNormalization(k, m)
  % depends on :
  % assocLegendreIntegral.m
  
  constant = (2*k +1) / (8*pi) *  assocLegendreIntegral(k, abs(m));
  
endfunction