function [res, t] = ray_distribution(n, seed, s)
%RaY_DISTRIBUTION Generates samples from an rayleigh distribution
%   using the inverse transform sampling method.
%
%   INPUTS:
%       n    - number of random samples to generate
%       seed - initial seed
%       s    - scale parameter of the rayleigh distribution
%
%   OUTPUTS:
%       res  - vector of generated rayleigh random variables
%       t    - execution time of the generation process
%
%   The rayleigh distribution is generated using:
%       X = r * sqrt(-2*log(U)), where U ~ Uniform(0,1)

tic

F_inv =@(theta, x) (theta*sqrt(-2*log(x)));
[res, ~] = inversion_method(n, F_inv, seed, s);

t = toc;
end