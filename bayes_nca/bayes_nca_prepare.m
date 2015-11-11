function [W_data,W_signs,Y_mean,Y_std,A_prior_mean,A_prior_std,B_prior_mean,B_prior_std,options] = bayes_nca_prepare(W_data,W_signs,Y_mean,Y_std,A_prior_mean,A_prior_std,B_prior_mean,B_prior_std,options)

% [W_data,W_signs,Y_mean,Y_std,A_prior_mean,A_prior_std,B_prior_mean,B_prior_std,options] = bayes_nca_prepare(W_data,W_signs,Y_mean,Y_std,A_prior_mean,A_prior_std,B_prior_mean,B_prior_std,options)
%
% complete data structure 'options'
% remove missing values in data Y and priors A and B

eval(default('options', '[]', 'W_signs', '[]', 'Y_std', '[]', 'A_prior_std', '[]', 'B_prior_std', '[]'));

if isempty(options ),       options     = struct; end 
if isempty(W_signs),        W_signs     = 0*W_data;   end
if isempty(Y_std  ),        Y_std       = nan*Y_mean; end
if isempty(A_prior_std  ),  A_prior_std = nan*A_prior_mean; end
if isempty(B_prior_std  ),  B_prior_std = nan*B_prior_mean; end

options_default = bayes_nca_default_options(A_prior_mean,B_prior_mean);
options = join_struct(options_default,options);

Y_mean(isnan(Y_mean))             = options.nan_Y_mean;
A_prior_mean(isnan(A_prior_mean)) = options.nan_A_prior_mean;
B_prior_mean(isnan(B_prior_mean)) = options.nan_B_prior_mean;
Y_std(isnan(Y_std))               = options.nan_Y_std;
A_prior_std(isnan(A_prior_std))   = options.nan_A_prior_std;
B_prior_std(isnan(B_prior_std))   = options.nan_B_prior_std;
