% Computation of coefficients normalizing the weight arising in SPECT
% in slice-by-slice approach. 

% w_{k,m} = (2k+1) / 4pi * (int_{-1}^{1}P_km(x) dx) * wspect(-m),
% where wspect(m) = int_{0}^{2pi} Wspect exp(im*phi) dphi

function nfactor = WeightsNormalization(k, m)
  % depends on : assocLegendreIntegral
  % assocLegendreIntegral.m
  
  nfactor = (2*k +1) / (4*pi) * assocLegendreIntegral(k, abs(m));
  
endfunction