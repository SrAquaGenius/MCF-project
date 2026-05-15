function [res, t, s_out] = exp_distribution(n, seed, r)
%EXP_DISTRIBUTION Generates samples from an exponential distribution
%   using the inverse transform sampling method.
%
%   INPUTS:
%       n     - number of random samples to generate
%       seed  - initial seed
%       r     - rate parameter of the exponential distribution
%
%   OUTPUTS:
%       res   - vector of generated exponential random variables
%       t     - execution time
%       s_out - next seed for reproducibility
%
%   The exponential distribution is generated using:
%       X = -log(U) / r, where U ~ Uniform(0,1)

if (nargout >= 2); tic; end

F_inv =@(theta, x) (-log(x)/theta);
if (nargout == 3)
    [res, ~, s_out] = inversion_method(n, F_inv, seed, r);
elseif (nargout <= 2)
    res = inversion_method(n, F_inv, seed, r);
end

if (nargout >= 2); t = toc; end
end