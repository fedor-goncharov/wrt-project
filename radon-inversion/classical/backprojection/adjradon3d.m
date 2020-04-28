function adjoint = adjradon3d(filename, ngrid, nphi, ntheta, nshift, interp_method=1, rsupp=1.0, expand_factor=2.0) 
  
% adjradon3d.m 

% Function takes on input Radon transforms of a certain function and computes its adjoint Radon transform in 3D.  
% The Radon data are assumed to be given on grid uniform in shifts and uniform-gauss in directions on the sphere.

% Input parameters : 
%     filename : file with Radon transforms
%     ngrid    : number of pixels in one direction in the XYZ-domain
%     nphi     : number of longitude angles [0, 2pi]
%     ntheta   : number of latitude angles [0, pi]
%     nshift   : 
%     interp_method : interpolation method, because Radon transforms are given on a discrete grid
%                     an interpolated data is needed for computation of the adjoint integral
%                     possible methdods : 1 is "linear" (default), 0 is "nearest"
%
%     rsupp    : length of the XYZ-domain in one direction (by default = 1.0)
%     expand_factor : compute adjoint transform on the grid on the size [expand_factor * ngrid]^3

% Return value : 
%     A three-dimensional matrix is returned of size (expand_factor*ngrid) per dimension, 
%     where elements stand for the adjoint Radon transform at given pixel.
  
  
  % read (weighted) Radon transforms from file
  rt = csvread(filename);                        % read values of the (weighted) Radon transforms
  rt = rt(:, size(rt)(2));
  rt = reshape(rt, ntheta, nphi, nshift);        % reshape into matrix [theta, phi, shift]  
  
  % set standard variables
  dphi   = 2 * pi / nphi;
  dshift = 2 * rsupp / (nshift-1);
  
  % set grid in Radon space
  % get equatorial angle
  angles_phi = (0 : nphi-1) * (2*pi / nphi);                         % directions 'phi' on circle [0, 2*pi)
  % get latitude angles (Gauss-Quadrature rule)
  [nodes_gauss, weights_theta] = lgwt(ntheta, -rsupp, rsupp);        % Gaussian nodes on [-1,1] and weights
  angles_theta = acos(nodes_gauss);                                  % directions 'theta' on (0, pi)
  clear nodes_gauss;
  shifts = linspace(-rsupp, rsupp, nshift)';
 
  % create grid XYZ
  lin = (-(expand_factor/2)*ngrid : (expand_factor/2)*ngrid) * dshift;   % adjoint may be computed on bigger grid to avoid border effects in derivations
  L = length(lin);
  [XX,YY] = meshgrid(lin, lin);
  XX = reshape(XX, L^2, 1);
  YY = reshape(YY, L^2, 1);
  
  % create grid shift-loop
  [FF, SS] = meshgrid(angles_phi, shifts);
  FFi = repmat(angles_phi, L^2, 1);
  
  % initi full grid
  adjoint_full_grid = zeros(L^3, 1);
  
  for i_z= 1 : L  % 2D for each slice
    
    points = [XX YY ones(L^2,1)*lin(i_z)];                   % size of the grid L x L
    adjoint_slice_grid = zeros(L^2, 1);
    
    printf("i_z=%d\n", i_z);
    fflush(stdout); 
    
    % integral over sphere for each z
    for i_theta = 1 : ntheta
      
        % get current angle
        theta = angles_theta(i_theta);
      
        % compute scalar products with direction for all points in the grid
        direction = [sin(theta)*cos(angles_phi); sin(theta)*sin(angles_phi); cos(theta) * ones(1, nphi)]; % size 3 x nphi
        SSi = points * direction;  % temporary shifts, matrix product is of size L^2 x nphi
      
        % values for the grid FF, SS (see above)
        VV = permute(rt(i_theta, :, :), [3 2 1]);
         
        II = qinterp2(FF, SS, VV, SSi, FFi, 2, 0);
      
        % add the term 
        adjoint_slice_grid += sum(II, 2) * dphi * weights_theta(i_theta);
    endfor
    
    % store computed adjoint for one z-slice 
    adjoint_full_grid(1 + (i_z-1)*L^2 : i_z * L^2) = adjoint_slice_grid;
  endfor 
  
  adjoint_full_grid = reshape(adjoint_full_grid, L, L, L);
endfunction