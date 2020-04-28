% Read a weight-coefficient on a grid [ngrid x ngrid x ngrid] from file.
% It is assumed that function is given in a CSV file with ordering 
% z -> x -> y (starting with a minimal coordinate)

function weight_coeff = ReadWeightCoeff3D(ngrid, filename)
  weight_coeff = zeros(ngrid, ngrid, ngrid);
  wcoeff_array = csvread(filename);
  for i = 1 : ngrid
    for j = 1 : ngrid
      for k = 1 : ngrid
        % cycle [z -> x -> | y | <- x <- z]
       weight_coeff(ngrid - k + 1, j, i) = wcoeff_array((i-1)*ngrid^2 + (j-1)*ngrid + k); 
      endfor
    endfor
  endfor
  % return weight_coeff
endfunction