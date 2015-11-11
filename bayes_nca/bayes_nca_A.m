% [A_post_mean, A_post_std] = bayes_nca_A(B, W_data, W_signs, Y_mean, Y_std, A_prior_mean, A_prior_std, A_max_abs, A_min_abs);
%
% Bayesian network component analysis, updating step for matrix A
%
% Outputs: A_post_mean, A_post_std (posterior mean and std. dev. for matrix elements A)
%
% No missing values (nan) are allowed in input matrices
%
% see bayes_nca

function [A_post_mean, A_post_std] = bayes_nca_A(B, W_data, W_signs, Y_mean, Y_std, A_prior_mean, A_prior_std, A_max_abs, A_min_abs);

nr          = size(A_prior_mean,1);
A_post_mean = zeros(size(A_prior_mean));
A_post_std  = zeros(size(A_prior_mean));

for it = 1:nr,
  x_mean    = Y_mean(it,:);
  ind_rel   = find(W_data(it,:));
  nc        = length(ind_rel);
  a_mean    = A_prior_mean(it,ind_rel);
  a_sign    = W_signs(it,ind_rel);
  dg_a_sign = diag(a_sign);
  ind_sign  = find(a_sign);
  B_rel     = B(ind_rel,:);
  Cx_inv    = diag(1./Y_std(it,:).^2);
  Ca_inv    = diag(1./A_prior_std(it,ind_rel).^2);
  M_ineq    = - dg_a_sign(ind_sign,:);
  b_ineq    = - A_min_abs * ones(length(ind_sign),1);
  M_ineq    = [M_ineq; -eye(nc); eye(nc)];
  b_ineq    = [b_ineq; A_max_abs * ones(nc,1); A_max_abs * ones(nc,1)];
  M         = B_rel * Cx_inv * B_rel' +  Ca_inv;
  M         = 0.5 * [M+M']; 
  f         = - [B_rel * Cx_inv * x_mean' +  Ca_inv * a_mean'];

  %if exist('cplexqp','file'),
  %  opt = cplexoptimset('Display','off'); 
  %  A_post_mean(it,ind_rel) = cplexqp(M, f, M_ineq, b_ineq,[],[],[],[],[],opt );
  %else
    opt = optimset('LargeScale','off','Algorithm','active-set','Display','off'); 
    A_post_mean(it,ind_rel) = quadprog(M, f, M_ineq, b_ineq,[],[],[],[],[],opt );
  %end

  if nargout > 1
    A_post_std(it, ind_rel) = sqrt(diag(inv(M)));
  end

end
