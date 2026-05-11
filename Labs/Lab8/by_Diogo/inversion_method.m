function [res, t] = inversion_method(n, f, seed, varargin)
%INVERSION_METHOD Generate the sample using for the given inverted function
%
%   INPUTS:
%       n    - number of samples
%       f    - inverse CDF function handle
%       seed - seed for LCG
%       args - parameters for f
%
%   OUTPUT:
%       res  - generated samples
%       t    - execution time

tic

[u, ~] = linear_congruential_generator(n, seed);
res = f(varargin{:}, u);

t = toc;
end