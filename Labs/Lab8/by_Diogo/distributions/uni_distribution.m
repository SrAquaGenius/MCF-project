function [res, t, s_out] = uni_distribution(n, seed, a, b)
%UNI_DISTRIBUTION Generates samples from an uniform distribution
%   using the inverse transform sampling method.
%
%   INPUTS:
%       n    - number of random samples to generate
%       seed - initial seed
%       a    - first value of the interval
%       b    - last value of the interval
%
%   OUTPUTS:
%       res  - vector of generated uniform random variables
%       t    - execution time of the generation process
%       s_out - next seed for reproducibility
%
%   The uniform distribution is generated using:
%       X = (b - a)*U + a, where U ~ Uniform(0,1)

if (nargout >= 2); tic; end

F_inv =@(a, b, x) ((b-a)*x + a);
if (nargout == 3)
    [res, ~, s_out] = inversion_method(n, F_inv, seed, a, b);
elseif (nargout <= 2)
    res = inversion_method(n, F_inv, seed, a, b);
end

if (nargout >= 2); t = toc; end
end