function output = laplace2dfft(u, ngrid, period=2.0) 
  
% laplace2dfft.m 

% Function applies to an image the laplacian in 2D. The 2d-laplacian is realized via FFT
% (look for "Notes on FFT-based differentiation" by Steven G. Johnson (MIT)). 

% Input parameters : 
%     u         : 2D image of size ngrid x ngrid
%     ngrid     : size of image in pixels in one dimension
%     perioid   : geometrical size of the period of the image

% Return value : 
%     A two-dimensional matrix of size ngrid x ngrid x ngrid
  
% A remark : 

% Since for the use of FFT (on which the implementation is based) the signal is assumed to be 
% periodic, application of the transform to a non-periodic image produces artifacts at the boundaries. 
% A way to evade the latter is perform laplace transform on a bigger image which embedds the initial one.
% To avoid singularities, the embedding should be smoothed by function stepfunction in the are of embedding.   
% Another is to use the method described in article : 
  
  Y = fftn(u); % 3D DFT of the signal 
  frequence_factor = zeros(ngrid, 1);
  for k = 0 : ngrid-1
    if (k <= (ngrid/2))
      frequence_factor(k + 1) = k;
    else 
      frequence_factor(k + 1) = k - ngrid;
    endif
   endfor
   
  freq_factor_dim1 = zeros(ngrid, ngrid);  
  freq_factor_dim2 = zeros(ngrid, ngrid);
  
  [freq_factor_dim1, freq_factor_dim2] = meshgrid(frequence_factor, frequence_factor);
  
  U = -(2*pi / period)^2*(freq_factor_dim1.^2 + freq_factor_dim2.^2);
  
  output = ifftn(U .* Y);
  
endfunction