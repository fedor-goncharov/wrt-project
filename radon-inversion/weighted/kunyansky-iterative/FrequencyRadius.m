% Computation of Nyquist frequency radius if a function is given 
% on a grid [-rsupp, rsupp]^3 with 'ngrid' points per dimension (counting 
% ends).

function rfreq = FrequencyRadius(ngrid, rsupp)
  
  dx = 2 * rsupp / (ngrid - 1);
  dfreq = 1 / (ngrid * dx);
  rfreq = (ngrid-1) * dfreq / 2;
  
endfunction