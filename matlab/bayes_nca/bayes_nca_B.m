% [B_post_mean, B_post_std] = bayes_nca_B(A, W_data, W_signs, Y_mean, Y_std, B_prior_mean, B_prior_std, verbose);
%
% Bayesian network component analysis: updating of matrix B;
%
% Output B_post_mean, B_post_std (posterior mean and std. dev. for matrix elements A)
%
% No missing values are allowed in input matrices
%
% see bayes_nca

function [B_post_mean, B_post_std] = bayes_nca_B(A, Y_mean, Y_std, B_prior_mean, B_prior_std, verbose);

nc          = size(B_prior_mean,2);
B_post_mean = zeros(size(B_prior_mean));
B_post_std  = zeros(size(B_prior_mean));

for it = 1:nc,

  if verbose, if mod(it,10)==0, fprintf('%d/%d ',it,nc); end; end

  y_mean  = Y_mean(:,it);
  b_mean  = B_prior_mean(:,it);
  Cy_inv  = diag(1./Y_std(:,it).^2);
  Cb_inv  = diag(1./B_prior_std(:,it).^2);
  B_post_mean(:,it) = [A' * Cy_inv * A + Cb_inv] \ [A' * Cy_inv * y_mean + Cb_inv * b_mean];

  if nargout > 1,
    B_post_std(:,it)  = sqrt(diag(inv([A' * Cy_inv * A + Cb_inv])));
  end

end
