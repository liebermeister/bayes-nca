function [M_mean,M_std] = tensor_mean_std(M,n_accept)

% [M_mean, M_std] = tensor_mean_std(M,n_accept)
%
% Compute mean and standard deviations over a 3-dim array
% 
% M is a tensor containing many sampled versions of a matrix; the third index is the sample index
% Use only the last 'n_accept' entries

ss           = size(M);
if n_accept > ss(3), n_accept  = ss(3); end

M_sum        = zeros(ss(1),ss(2));
M_square_sum = M_sum;

for it = ss(3)-n_accept+1:ss(3),
  M_sum        = M_sum        + squeeze(M(:,:,it)); 
  M_square_sum = M_square_sum + squeeze(M(:,:,it).^2); 
end

M_mean = M_sum/n_accept;
M_std  = sqrt(M_square_sum/n_accept);
