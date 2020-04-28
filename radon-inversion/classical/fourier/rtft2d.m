% rtft2d.m

function [nodes, values, jacobian_weights] = rtft2d(filename, nphi, nshift, rsupp=1.0, padding = 4.0)
% depends on cartprod.m. centered_frequencies.m

% This script reads data given by ray transforms in 2D from file and performs 
% 1D Fourier transform along shift argument. Returns values (values) of the 
% 1D Fourier integral and respective nodes (nodes) with volumes corresponding to nodes (jacobian_weights). 

% Usage of the script:

%   filename          : file with Radon transforms [sigma, phi]
%   nphi              : number of projections in azimuth angle
%   nshift            : number of rays per one direction, shift [-1, 1]
%   rsupp             : radius of the support of the test-function
%   padding           : multiplying factor for appending data with zeros
 
% It is assumed that angles 'phi' and 'shifts' are uniformly spaced 
% in their intervals (0,2pi], [-1,1], respectively.
  
% OUTPUT : nodes - nodes in 3D spaces, where Fourier Transform of the function is known
%          values - values of FT of the function
%          jacobian_weights - volume of cells for each node in fourier space  
  
  fd = fopen(filename, "rb");
  rt_matrix = reshape(fread(fd, 'double'), nphi, nshift);
  fclose(fd);
  
  dphi = 2*pi/ nphi;
  dshift = 2 * rsupp / (nshift-1);
  padding_size = floor(padding * nshift);          % total number of points with zero padding (should be odd)
  ntotal = nshift + 2 * padding_size;
  dfreq   = 1 / (ntotal * dshift);                 % discretization step in frequency domain
  
  phi_vec = vec((0 : nphi-1) * dphi);              % directions 'phi' on circle [0, 2*pi)
  
  frequencies = vec((0 : ntotal-1) * dfreq);
  frequencies_centered = vec(centered_frequencies(ntotal, dfreq));
  
  % Nyquist sinc filter 
  half_nyq_frequency = ntotal * dfreq / 2;
  sinc_filter = (sinc(frequencies_centered ./ half_nyq_frequency)).^2;
  
  nodes = [];            % nodes in frequency domain
  values = [];           % values of Fourier integral in nodes 
  jacobian_weights = []; % jacobian multipliers for each value
  
  printf("Generating nodes in Fourier space...\n");
  fflush(stdout);
  
  for i_phi = 1 : nphi
      % current angle phi
      phi = phi_vec(i_phi);
      
      % append nodes
      direction = [cos(phi) sin(phi)];
      nodes = [nodes; frequencies_centered * direction];
      
      % append jacobian
      jacobian_weight = 0.5 * dfreq * dphi * abs(frequencies_centered);
      
      % application of nyquist sinc filter 
      jacobian_weight = jacobian_weight .* sinc_filter;
      
      % append jacobian weights of cells
      jacobian_weights = [jacobian_weights; jacobian_weight];
      
      % append 1D Fourier transform of RT using FFT along fixed direction
      rt_at_direction = vec(rt_matrix(i_phi, :));                                      % data on fixed direction
      fft_vec = fftshift(fft(ifftshift(padarray(rt_at_direction, padding_size))));    % 1D Fourier integral using 1D FFT
      ft1d_vec = dshift * fft_vec;    % centralizing frequencies and jacobian correction
      
      values = [values; vec(ft1d_vec)];       % 'values' of the Fourier transform of test function at 'nodes'
  end
  
  printf("Done. %d primary nodes created.\n", size(values, 1));
  printf("Appending exterior nodes...");
  fflush(stdout);
  
  % stabilization - append nodes outside of the ball of radius of Nyquist/2 frequency
  append_nodes = cartprod(frequencies_centered, frequencies_centered); % all points in Fourier domain
  append_nodes = [append_nodes, sqrt(sum(append_nodes.^2, 2))];                            % append norms column to the right
  append_nodes = append_nodes(append_nodes(:, 3) > ((ntotal/2 + 1) * dfreq), :);           % choose nodes outside the ball
  append_nodes = append_nodes(:, [1 2]);                                                   % remove norm column
  nodes = [nodes; append_nodes];
  size_append = size(append_nodes, 1);  % remember the number of nodes that have been appended
  clear append_nodes;
  
  printf("Done. %d nodes appended.\n", size_append);
  fflush(stdout);
  
  values = [values; zeros(size_append, 1)];                               % append zero values
  jacobian_weights = [jacobian_weights; ones(size_append, 1)*(dfreq^2)];  % append square volumes
    
  printf("Done.\n");
  fflush(stdout); 
   
end 
