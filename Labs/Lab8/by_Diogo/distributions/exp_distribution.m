function [res, t] = exp_distribution(n, seed, r)
%EXP_DISTRIBUTION Generates samples from an exponential distribution
%   using the inverse transform sampling method.
%
%   INPUTS:
%       n    - number of random samples to generate
%       seed - initial seed
%       r    - rate parameter of the exponential distribution
%
%   OUTPUTS:
%       res  - vector of generated exponential random variables
%       t    - execution time of the generation process
%
%   The exponential distribution is generated using:
%       X = -log(U) / r, where U ~ Uniform(0,1)

tic

F_inv =@(theta, x) (-log(x)/theta);
[res, ~] = inversion_method(n, F_inv, seed, r);

t = toc;
end