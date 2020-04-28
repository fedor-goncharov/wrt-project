% Generates even non-zero spherical harmonics [1, 2*degree] 
% on a cubic grid [-rsupp,rsupp]^3.

% ngrid : number of points (including ends) on [-rsupp, rsupp]
% rsupp : radius of the domain
% degree: maximal degree of even harmonics (max_degree = 2 * degree)

function SphHarmonics = GenSphHarmonics(ngrid, rsupp, degree)
% depends on HarmonicYkm.m  
  
  linx = linspace(-rsupp, rsupp, ngrid);
  liny = linspace(rsupp, -rsupp, ngrid);  % axis Y must be inverted
  linz = linspace(-rsupp, rsupp, ngrid);
  
  [X,Y,Z] = meshgrid(linx, liny, linx);   % grid in the unit cube in 3D
                                                                               
  [PHI, THETA, ] = cart2sph(X, Y, Z);  % grid in spherical coordinates
  THETA = pi/2 - THETA;                % shift angle from [-pi/2, pi/2] to [0,pi]
  
  % array of even non-zero spherical harmonics on the grid
  SphHarmonics = zeros(ngrid, ngrid, ngrid, (degree + 1)*(2*degree + 1) - 1); 
                                                                                
  current_harmonic = 1;  % counter for even non-zero spherical harmonics
  for k = 1 : degree    
    for m = -2*k : 2*k   % iterate over even harmonics
      SphHarmonics(:, :, :, current_harmonic) = HarmonicYkm(2*k, m, THETA, PHI);
      current_harmonic = current_harmonic + 1;
    endfor
  endfor
  
endfunction
  