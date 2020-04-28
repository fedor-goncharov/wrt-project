% Reconstruction of a test function in 3D from weighted Radon transforms 
% using iterative approach of F. Goncharov, R. Novikov (2017)

function [test_function, sigma] = IterativeWRadonInv3D(general_params, init_params, 
  weights_params, stopping_params)
% IterativeWRadonInversion3D.m

% depends on ReadGridFunc3D.m,
%            ReadWeightCoeff3D.m, 
%            StopRuleCheck.m, 
%            IntOperatorQ.m,
%            GenCutoffFunc.m,
%            NormalizationConstant.m,
%            ParametrizationFactor.m

% general_params
%     ngrid -- integer, number points on the grid per dimension (on [-1,1], including the ends) 
%     filename -- file with the initial test-function R^{-1}R_W(F)
% init_data_params
%     file_name_init -- file with the data of the test-function on the grid [ngrid x ngrid x ngrid]
% weight_data_params:
%     degree -- integer, highest order of harmonics used
%     file_name_array -- array of files where coefficients of the weight are stored
%                        coefficients are tensors [ngrid x ngrid x ngrid]
% stopping_params
%     flag -- FLAG_NITERATIONS, FLAG_PRECISION, FLAG_COMPOSITE
%     flag_values{1,2} -- integer, double
%         FLAG_NITERATIONS -- algorithm stops after 'stopping_iteration' iterations
%                             stopping_iteration = flag_values{1} 
%         FLAG_PRECISION   -- algorithm stops reaching given precision |u_n - u_(n-1)| / | u_n-1 | < stopping_precision
%                             stopping_precision = flag_values{2}
%         FLAG_COMPOSITE   -- algorithm stops either after 'stopping_iteration' iterations either 
%                             reaching given precision |u_n - u_(n-1)| / | u_n-1 | < stopping_precision
%                             stopping_iteration = flag_values{1} 
%                             stoppping_precision = flag_values{2}
%
%     for a 'good' choice of number of harmonics the algorithm converges geometrically to quasi-optimum
%     set eps as precision to quasi-solution - it implies the exact number of iterations you need to do
%     to reach the quasi-optimum 

% ----------------------------read initial function ----------------------------
  ngrid = general_params{1};      % size of the grid
  rsupp = general_params{2};
  
  % initial point R^{-1}R_W(F)
  F0_filename = general_params{3};
  F0 = ReadGgridFunc3D(ngrid, F0_filename);                                     

  % ---------------------------- read weight coeffs ----------------------------
  degree = weights_params{1};                                                   % degree ( >= 1) of iterative approximation
  weights_filenames = weights_params{2};                                        % files with weights
  
  weights = zeros(ngrid, ngrid, ngrid, (degree + 1)*(2*degree + 1));            % init array of weight coeffs.
  for i = 1 : (degree + 1)*(2*degree + 1)                                       
    weights(:, :, :, i) = ReadWeightCoeff3D(ngrid, weights_filenames(i));       % read weights (corresponding to harmonics)                                                                     
  endfor
% ------------------------------------------------------------------------------


% --------- normalization of weights (w_(2k,m) / w_(0,0)) * correction_constant-
% --------------------------- precomputation of all harmonics ------------------
  
  
  % weight of order {0,0}
  w0 = weights(:, :, :, 1) / (2*pi);
  
  % weight coefficients without (0,0) term
  Nw_normalized = (degree + 1)*(2*degree + 1) - 1; % number of normalized weights
  weights_normalized = zeros(ngrid, ngrid, ngrid, Nw_normalized);
  for i = 1 : Nw_normalized
    weights_normalized(:, :, :, i) = weights(:, :, :, i+1);
  endfor
  
  % weights normalization constants
  current_weight = 1;
  for k = 1 : degree
    for m = 1 : 4*k + 1
      weights_normalized(:, :, :, current_weight) = weights_normalized(:, :, :, current_weight) ...
                                                    * NormalizationConstant(2*k, m - 2*k - 1) * ParametrizationFactor(2*k, m - 2*k - 1);
      current_weight = current_weight + 1;
    endfor
  endfor
  
  % w0 normalization of other weights
  for i = 1 : Nw_normalized
    weights_normalized(:, :, :, i) = weights_normalized(:, :, :, i) ./ w0;      % w0 normalization
  endfor
  
  % cutoff function of the unit ball
  cutoff = GenCutoffFunc(ngrid, rsupp);
  
% -------------- spherical harmonics [2, 2 * degree] on [ngrid x ngrid x ngrid] 
% -------------- create harmonics in Fourier space -----------------------------

  rfreq = FrequencyRadius(ngrid, rsupp);
  harmonics = GenSphHarmonics(ngrid, rfreq, degree);

%--------------- Initialization of functions before iterations -----------------
 
% initialization of a starting point for two layer scheme

  u_new = zeros(ngrid, ngrid, ngrid);  %  next layer
  u_old = F0;                          %  current layer
  u0 = F0;
  
% ------------------------ Cycle of iterations ---------------------------------

  iteration = 0;  % iteration counter
  do 
    u_new = u0 - IntOperatorQ(ngrid, degree, weights_normalized, harmonics, cutoff, rsupp, rfreq, u_old);
    u_old = u_new;
    iteration = iteration + 1;
  
    cycle_params = cell(2);
    cycle_params{1} = iteration; 
    cycle_params{2} = norm(u_new - u_old, 2) / norm(u_old, 2);                  % current precision from the new iteration
  while (StopRuleCheck(stoppping_params, cycle_params) == false)
# ------------------------ end of cycle computations ---------------------------

  test_function = u_old ./ w0;   % return value of the test-function ngrid x ngrid x ngrid
  
  % computation of sigma
  sigma = 0;
  for i = 1 : Nw_normalized
    array = reshape(weights_normalized(:, :, :, i), ngrid^3);
    sigma += norm(array, 'inf');
  endfor
  
endfunction 

