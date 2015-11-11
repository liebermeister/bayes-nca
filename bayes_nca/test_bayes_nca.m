addpath(genpath('~/projekte/metabolic_NCA/'))

% Small example
W_data  = [1 0; 1 1; 0 1];
W_signs = [1 0; 1 -1; 0 1];
A_true  = [1 0; 1 -1; 0 1];
B_true  = [sin(2*pi*[0:0.02:1]); cos(2*pi*[0:0.02:1])];

% Large example
W_signs = rand(50,10) .* [rand(50,10) >0.9];
W_data  = double(W_signs~=0);
A_true  = W_signs;
B_true  = [sin([1:5]'*2*pi*[0:0.02:1]); cos([1:5]'*2*pi*[0:0.02:1])];

Y_true       = A_true * B_true;

Y_mean       = Y_true;
Y_std        = 0.1*ones(size(Y_true));
A_prior_mean = W_signs;% zeros(size(A_true)); % A_true; % 
A_prior_std  = 1 * double(A_true~=0);
B_prior_mean = zeros(size(B_true)); % B_true; % 
B_prior_std  = 1*ones(size(B_true));

options = struct('n_it_max', 100, 'n_accept', 80, 'method', 'gibbs sampling', 'graphics_flag', 1); 

result  = bayes_nca(W_data, W_signs, Y_mean, Y_std, A_prior_mean, A_prior_std, B_prior_mean, B_prior_std, options);

figure(1); 
bayes_nca_graphics(A_prior_mean, B_prior_mean, Y_mean, A_true, B_true, Y_true);

figure(2); 
bayes_nca_graphics(A_prior_mean, B_prior_mean, Y_mean, result.A_post_mean, result.B_post_mean, result.Ypred);
