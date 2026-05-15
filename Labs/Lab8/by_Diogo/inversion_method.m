function [res, t, s_out] = inversion_method(n, f, seed, varargin)
%INVERSION_METHOD Generate the sample using for the given inverted function
%
%   INPUTS:
%       n     - number of samples
%       f     - inverse CDF function handle
%       seed  - seed for LCG
%       args  - parameters for f
%
%   OUTPUT:
%       res   - generated samples
%       t     - execution time
%       s_out - seed out for reproducibility

if (nargout >= 2); tic; end

u = linear_congruential_generator(n, seed);
res = f(varargin{:}, u);

if (nargout >= 2); t = toc; end
if (nargout == 3); s_out = u(end); end
end