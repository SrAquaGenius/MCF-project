function [res, t, s_out] = nrm_distribution(n, seed, mu, sigma)
%NRM_DISTRIBUTION Generates samples from a normal distribution using the
%   Box-Muller sampling method.
%
%   INPUTS:
%       n     - number of random samples to generate
%       seed  - initial seed
%       mu    - mean value of the normal distribution
%       sigma - standard deviation of the normal distribution
%
%   OUTPUTS:
%       res   - vector of generated exponential random variables
%       t     - execution time
%       s_out - next seed for reproducibility

if (nargout >= 2); tic; end

if (nargout == 3)
    [z, ~, s_out] = box_muller_method(n, seed);
elseif (nargout <= 2)
    z = box_muller_method(n, seed);
end

res = mu + sigma * z;

if (nargout >= 2); t = toc; end
end