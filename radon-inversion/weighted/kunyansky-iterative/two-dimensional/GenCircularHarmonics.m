function CircHarmonics = GenCircularHarmonics(ngrid, rsupp, degree)
% depends on HarmonicYkm.m  
  
  linx = linspace(-rsupp, rsupp, ngrid);
  liny = linspace(rsupp, -rsupp, ngrid);  % axis Y must be inverted
  
  [X,Y] = meshgrid(linx, liny);   % grid in the unit cube in 2D
                                                                               
  [THETA, R] = cart2pol(reshape(X, ngrid^2, 1), reshape(Y, ngrid^2, 1));  % grid in polar coordinates
  %THETA = pi/2 - THETA;        % shift angle from [-pi/2, pi/2] to [0,pi]
  
  
  % array of even non-zero spherical harmonics on the grid
  CircHarmonics = zeros(ngrid, ngrid, 2*degree); 
                           
  index = -2*degree;                           
  for k = 1 : 2*degree    
      index  = index + 2*(k-1);
      if (index == 0)
        index = index + 2;
      endif
      CircHarmonics(:, :, k) = reshape(exp(i*index*THETA), ngrid, ngrid);
  endfor
  
endfunction
