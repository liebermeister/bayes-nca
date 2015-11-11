% Functions for Bayesian Network Component Analysis
%
% The NCA problem is defined by the structure of the connection matrix 
% (conections in W_data, signs in W_signs), mean values and standard
% deviations for the data matrix X, and priors (mean values and standard
% deviations) for the connection matrice A and the factor matrix B
%
% NCA is computed iteratively as proposed by Liao et al, but in each step,
% the posterior of A or B is computed. 
% In the next iteration step, one can either use the posterior mean (just 
% like in usual NCA) or sample from the posterior (Gibbs sampling).
%
% Main function: bayes_nca.m
% Test example:  test_bayes_nca.m
