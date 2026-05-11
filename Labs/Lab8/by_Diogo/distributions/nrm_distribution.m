function [res, t] = nrm_distribution(n, seed, mu, sigma)
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
%       t     - execution time of the generation process

tic

[x, y, ~] = box_muller_method(n, seed);
z = [x(:); y(:)];
res = mu + sigma * z;

t = toc;
end