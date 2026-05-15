function [res, t, s_out] = ray_distribution(n, seed, s)
%RaY_DISTRIBUTION Generates samples from an rayleigh distribution
%   using the inverse transform sampling method.
%
%   INPUTS:
%       n     - number of random samples to generate
%       seed  - initial seed
%       s     - scale parameter of the rayleigh distribution
%
%   OUTPUTS:
%       res   - vector of generated rayleigh random variables
%       t     - execution time
%       s_out - next seed for reproducibility
%
%   The rayleigh distribution is generated using:
%       X = r * sqrt(-2*log(U)), where U ~ Uniform(0,1)

if (nargout >= 2); tic; end

F_inv =@(theta, x) (theta*sqrt(-2*log(x)));

if (nargout == 3)
    [res, ~, s_out] = inversion_method(n, F_inv, seed, s);
elseif (nargout <= 2)
    res = inversion_method(n, F_inv, seed, s);
end

if (nargout >= 2); t = toc; end
end