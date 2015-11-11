function [result,options] = bayes_nca(W_data,W_signs,Y_mean,Y_std,A_prior_mean,A_prior_std,B_prior_mean,B_prior_std,options)

% result = bayes_nca(W_data,W_signs,Y_mean,Y_std,A_prior_mean,A_prior_std,B_prior_mean,B_prior_std,options)
%
% Bayesian NCA; use POSTERIOR MEAN of A or B in every iteration step 
%
%Inputs
%  W_data, W_signs          : matrices of known connections and known connection signs
%  Y_mean, Y_std            : matrices of data values (means and standard deviations)
%  A_prior_mean, A_prior_std: matrices of A prior values (means and standard deviations)
%  B_prior_mean, B_prior_std: matrices of B prior values (means and standard deviations)
%  options                  : for defaults, see bayes_nca_default_options.m
%
% The arguments W_data and Y_mean are mandatory. All others can be replaced by '[]'
%
%Outputs: result, options

eval(default('verbose','0','W_signs','[]','Y_std','[]','A_prior_mean','[]','A_prior_std','[]','B_prior_mean','[]','B_prior_std','[]','options','struct'));

% Default values (if no values have been provided)

if isempty(W_signs),      W_signs      = zeros(size(W_data)); end
if isempty(A_prior_mean), A_prior_mean = zeros(size(W_data)); end
if isempty(A_prior_std),  A_prior_std  = double(W_data~=0);   end
if isempty(B_prior_mean), B_prior_mean = zeros(size(W_data,2),size(Y_mean,2)); end
if isempty(B_prior_std),  B_prior_std  = ones(size(W_data,2),size(Y_mean,2)); end
if isempty(Y_std),        Y_std        = 0.5 * mean(abs(Y_mean(:))) * abs(ones(size(Y_mean))); end

options = join_struct(bayes_nca_default_options(A_prior_mean, B_prior_mean),options);


% -------------------------------------------------------
% If more than one repetition is requested, call this function multiple times and stop


if options.repeat > 1,
  display(sprintf('Running bayesNCA %d times with different random seeds. Only the best result is reported'));
  best_ssr = inf;
  for it = 1:options.repeat,
    fprintf('\nRun %d/%d\n',it,options.repeat);
    my_options        = options;
    my_options.repeat = 1;
    my_options.seed   = rand;
    my_result = bayes_nca(W_data,W_signs,Y_mean,Y_std,A_prior_mean,A_prior_std,B_prior_mean,B_prior_std,my_options);
    if my_result.ssr_Y + my_result.ssr_B < best_ssr,
      result   = my_result;
      best_ssr = my_result.ssr_Y + my_result.ssr_B;
    end
  end
  return
end


% -------------------------------------------------------
% Initialise

[W_data,W_signs,Y_mean,Y_std,A_prior_mean,A_prior_std,B_prior_mean,B_prior_std,options] = bayes_nca_prepare(W_data,W_signs,Y_mean,Y_std,A_prior_mean,A_prior_std,B_prior_mean,B_prior_std,options);

if isfinite(options.seed), randn('state',options.seed); end

if isfield(options,'A_init'), 
  A = options.A_init; 
else,
  A = W_data .* randn(size(W_data));
end

if isfield(options,'B_init'), 
  B = options.B_init; 
else,
  B = randn(size(B_prior_mean));
end

ssr_Y_list = [];

if options.graphics_flag,
  bayes_nca_graphics(A_prior_mean, B_prior_mean, Y_mean, A, B, A*B, ssr_Y_list); drawnow
end


% -------------------------------------------------------
% Run Bayesian NCA

switch options.method
  
  case 'posterior mode',
    fprintf(' Iteration ');

    for it = 1:options.n_it_max,
      fprintf('%d/%d ',it,options.n_it_max)

      B = bayes_nca_B(A, Y_mean, Y_std, B_prior_mean, B_prior_std,options.verbose);

      % algorithm can get stuck when there are zero values in B
      % replace such values by small nonzero (random) values 
      ind_zero = find( B(:) ==0); 
      B(ind_zero) = 0.001 * randn(length(ind_zero),1); 

      A = bayes_nca_A(B, W_data, W_signs, Y_mean, Y_std, A_prior_mean, A_prior_std, options.A_max_abs, options.A_min_abs);

      % algorithm can get stuck when there are (unintended) zero values in A
      % replace such values by small nonzero (random) values 
      ind_zero = find( A(find(W_data~=0))==0);  
      A(ind_zero) = 0.01 * W_data(ind_zero) + 0.05 * W_signs(ind_zero); 

      if options.graphics_flag,
	ssr_Y_list = [ssr_Y_list, nansum(nansum([[A*B-Y_mean]./Y_std].^2))];
        bayes_nca_graphics(A_prior_mean, B_prior_mean, Y_mean, A, B, A*B, ssr_Y_list); drawnow
      end

    end

    % Final iteration
    [B,B_std] = bayes_nca_B(A, Y_mean, Y_std, B_prior_mean, B_prior_std,options.verbose);
    [A,A_std] = bayes_nca_A(B, W_data, W_signs, Y_mean, Y_std, A_prior_mean, A_prior_std, options.A_max_abs, options.A_min_abs);
    
    result.A_post_mean = A;
    result.B_post_mean = B;
    result.A_post_std  = A_std;
    result.B_post_std  = B_std;
    result.Y_pred      = A * B;
    
  case 'gibbs sampling',
 
    A_list      = [];
    B_list      = [];
    Y_pred_list = [];
    A_post_mean = A;

    fprintf(' Iteration ');

    for it = 1:options.n_it_max,
      fprintf('%d/%d ',it,options.n_it_max)

      [B_post_mean, B_post_std] = bayes_nca_B(A_post_mean, Y_mean, Y_std, B_prior_mean, B_prior_std,options.verbose);
      B = B_post_mean + B_post_std .* randn(size(B_post_std));

      [A_post_mean, A_post_std] = bayes_nca_A(B_post_mean, W_data, W_signs, Y_mean, Y_std, A_prior_mean, A_prior_std, options.A_max_abs, options.A_min_abs);
      % since sampling with sign constraints would be hard .. use approximation:
      % sample using true posterior mode (satisfying sign constraints) 
      % and covaraince matrix (not satisfying sign constraints) 
      % and then set all entries with wrong sign to 0 
      A = A_post_mean + A_post_std .* randn(size(A_post_std));
      A(find(A.*W_signs<0)) = 0;

      Y_pred = A_post_mean * B_post_mean;
      A_list(:,:,it) = A; 
      B_list(:,:,it) = B; 
      Y_pred_list(:,:,it) = Y_pred; 

      if options.graphics_flag,
	ssr_Y_list = [ssr_Y_list, nansum(nansum([[A*B-Y_mean]./Y_std].^2))];
        bayes_nca_graphics(A_prior_mean, B_prior_mean, Y_mean, A, B, A*B, ssr_Y_list);  drawnow
      end

    end

    [A_sample_mean, A_sample_std] = tensor_mean_std(A_list,options.n_accept);
    [B_sample_mean, B_sample_std] = tensor_mean_std(B_list,options.n_accept);
    %%[Y_sample_mean, Y_sample_std] = tensor_mean_std(Y_pred_list,options.n_accept);
    
    result.A_post_mean    = A_sample_mean;
    result.A_post_std     = A_sample_std;
    result.B_post_mean    = B_sample_mean;
    result.B_post_std     = B_sample_std;
    result.Y_pred         = A_sample_mean * B_sample_mean;  
    result.A_list         = A_list;
    result.B_list         = B_list;
    result.Y_pred_list    = Y_pred_list;  

  otherwise, error('unknown method');
    
end

% --------------------------------------------------
% Collect sum of squared residuals

result.ssr_Y = nansum(nansum([[result.Y_pred-Y_mean]./Y_std].^2));
result.ssr_A = nansum(nansum([[result.A_post_mean - A_prior_mean]./A_prior_std].^2));
result.ssr_B = nansum(nansum([[result.B_post_mean - B_prior_mean]./B_prior_std].^2));
