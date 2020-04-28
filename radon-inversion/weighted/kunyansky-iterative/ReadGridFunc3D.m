% Reading a test-function on a grid [ngrid x ngrid x ngrid] from file

% It is assumed that function is given in a CSV file with ordering 
% x -> y -> z (starting with a minimal index)

function func = ReadGridFunc3D(ngrid, file)
  
  func = zeros(ngrid, ngrid, ngrid);
  func_array = csvread(file);
  for i = 1 : ngrid
    for j = 1 : ngrid
      for k = 1 : ngrid
        % cycle [x -> y -> | z | <- y <- x]
        func(ngrid - j + 1, i, k) = func_array((i-1)*ngrid^2 + (j-1)*ngrid + k); 
      endfor
    endfor
  endfor
  % return grid_function
end