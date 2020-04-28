% Value of the integral of associated legendre polynomial of type (m,n) 
% along interval [-1,1]. The polynomials are assumed to be Schmidt 
% semi-normalized. 

% The formula is taken from: "The integral of the associated Legendre function" 
% by Donald W. Jepsen, Eugene F. Haugh, and Joseph O. Hirshfielder

% Note that in the paper above no Shortley factor (-1)^m is used in the definition
% of polynomials. It is corrected in the case of (m even, n even).

function Rkm = assocLegendreIntegral(k, m)
  % for zero or positive (k, m), k \geq m
  Rkm = 0;
  if (1 == mod(m + k, 2))                                                       % m + n odd
    Rkm = 0;                                                                    
  elseif (0==mod(k,2) && 0==mod(m,2))                                           % m even, n even
    Rkm = (2*m/k) * (factorial(k/2)^2 * factorial(k+m)) / ...
      (factorial((k-m)/2)^2 * factorial((k + m)/2) * factorial(k+1));
  else                                                                          % m odd, n odd
    Rkm = (m*pi)/(k*2^(2*k + 1)) * (factorial(k+m) * factorial(k+1)) / ...
      (factorial((k+1)/2)^2 * factorial((k-m)/2) * factorial((k+m)/2));
  end
  
  % correction constant for Schmidt semi-normalized polynomials
  schmidt_corrector = 1;
  if (m > 0)
    schmidt_corrector = (-1)^m * sqrt(2 * factorial(k-m) / factorial(k + m));   % Schmidt semi-normalization correction
  endif
  Rkm = Rkm * schmidt_corrector;
endfunction