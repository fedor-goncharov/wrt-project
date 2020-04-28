% Characteristic function of a unit ball in the unit cube in 3D

% ngrid : number of uniform points (including ends) on [-1, 1]
% Returns matrix of size ngrid x ngrid x ngrid

function output = GenCutoff(ngrid, rsupp)
    
  lin = linspace(-1, 1, ngrid);
  [X,Y,Z] = meshgrid(lin, lin, lin);
  V = X.^2 + Y.^2 + Z.^2;
  
  V(V > rsupp^2) = 0;
  V(V > 0) = 1;
  
  output = V;
  
endfunction