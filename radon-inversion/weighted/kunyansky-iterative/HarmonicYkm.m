% Computation of Schmidt semi-normalized spherical harmonics

% theta : latitude angle [0, pi] (angle from positive OZ axis) 
% phi   : longitude angle [0, 2pi] (angle from positive OX axis)

function Ykm = HarmonicYkm(k, m, theta, phi)
  % degree of harmonic k (positive integer), order of harmonic m in [-k, k]

  Pk = legendre(k, cos(theta), 'sch');  % associated Legendre polyomials (Schmidt semi-normalized)
  Pkm = reshape(Pk(abs(m)+1, :, :, :), size(phi));
  Ykm = Pkm .* exp(i*m*phi);
  
endfunction