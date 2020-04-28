% Reconstruction of a test function in 3D from its Fouirer transforms

% nfft_reconstruct_3d.m
% depends NFFT library (see Chemnitz-TU for NFFT)
  
% Script performs reconstruction of a function from its Fourier transforms at 'nodes' 
% with values 'values'. Mathematically it works as a discretized version of (inverse) Fourier integral 
% in 3D. Result is given as a 3D matrix of size : ngrid x ngrid x ngrid, which, in turn, is a 
% grid on [-1,1)x[-1,1)x[-1,1). 

% Usage of the script
% ngrid : number of points (as intervals) in [-1.0, 1.0)
% nodes : matrix of size (number_of_nodes x 3), where '3' stands for [x, y, z], these are the 
%         points in space where Fourier transform is known
% values : vector of size (number_of_nodes) values of Fourier transform of a function in 'nodes'
% jacobian_weights : volumes at 'nodes' in the Riemann summ of discretized Fourier integral

% Example : >> test_function = nfft_reconstruct_3d(128, nodes, values, jacobian_weights);
function test_function = nfft_reconstruct_3d(ngrid, nodes, values, jacobian_weights)

% Jacobian weightening
  summands = jacobian_weights .* values; 
  
% Inversion using NFFT
  
  %normalization of 'nodes' so that they belong to torus [-1/2, 1/2)^3
  deltax = 2.0 /  ngrid;  % normalization on Nyquist frequency
  nodes_normalized = nodes * deltax;
  nodes_normalized = nodes_normalized';

  % Initialisation of plan
  n=2^(ceil(log(ngrid)/log(2))+1);
  plan = nfft_init_guru(3, ngrid, ngrid, ngrid, size(nodes_normalized,2), n, n, n, 8, 
  NFFT_OMP_BLOCKWISE_ADJOINT, FFTW_ESTIMATE);
  
  % set nodes in plan
  nfft_set_x(plan, nodes_normalized);
  
  % precomputations
  nfft_precompute_psi(plan);
  
  % set Fourier coefficients
  nfft_set_f(plan, summands);
  
  % inversion
  nfft_adjoint(plan);
  
  % return test-function as a 3d matrix
  test_function = reshape(real(nfft_get_f_hat(plan)), ngrid, ngrid, ngrid);
  test_function = permute(test_function, [3 2 1]);
  test_function = permute(test_function, [2 1 3]);
  
  nfft_finalize(plan);
end