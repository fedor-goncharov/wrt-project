% two-dimensional iterative reconstructions

wm2_real = wm2_strong_real_slice;
wm2_imag = wm2_strong_imag_slice;

% w0, order 0
w0 = w0_att_strong_slice;   % it should be normalized

% weights of order 2
w2_real = w2_strong_real_slice;
w2_imag = w2_strong_imag_slice;

% full complex weights
wm2 = wm2_real + 1i * wm2_imag;
w2  = w2_real + 1i * w2_imag;

% one array of weights
% note: don't forget correction because weight W(phi) = w(phi + pi/2) ???

ngrid = 128;
order = 1;

weights = zeros(ngrid, ngrid, 2*order + 1);
weights(:, :, 1) = w0; % w0, normalization weight
weights(:, :, 2) = w2 / (2*pi);
weights(:, :, 3) = wm2 / (2*pi);

w00 = w0 / (2*pi);

% normalization of weights by w0
weights_normalized = zeros(ngrid, ngrid, 2*order);
for i = 1 : 2*order
  weights_normalized(:, :, i) = weights(:, :, i + 1) ./ w00;
endfor
  
% weights normalization constants
%for k = 1 : 2*order
%    weights_normalized(:, :, k) = weights_normalized(:, :, k);
%endfor


% create harmonics in Fourier space
rsupp = 1.0;
dx = 2 * rsupp / (ngrid - 1);
dfreq = 1 / (ngrid * dx);
rfreq = (ngrid-1) * dfreq / 2;
harmonics = GenCircularHarmonics(ngrid, rfreq, order);


% cutoff function
cutoff = cutoff_slice;

% iterations
u0 = reconstruction2d_ph1_att_strong_nweak;
u_new = zeros(ngrid, ngrid);  %  next function layer
u_old = u0;   % current function layer


% iteration
for iteration = 1 : 4
  u_new = u0 - IntOperatorQ2D(ngrid, order, weights_normalized, harmonics, cutoff, u_old);
  u_old = u_new;
endfor

f0 = reconstruction2d_ph1_att_strong;
error = norm(((u_new - f0) ./ w00), 'fro') / norm(f0 ./ w00, 'fro');
error


s = 0;
for i = 1 : 2*order
  array = reshape(weights_normalized(:, :, i) .* cutoff, ngrid^2, 1); 
  printf("Normalized weight (%d) : %f\n", i, norm(array, 'inf')); 
  s += norm(array, 'inf');  
endfor
