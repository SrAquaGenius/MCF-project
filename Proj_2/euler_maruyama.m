function [X] = euler_maruyama(a, b, x0, T, N, dW)
%EULER_MARUYAMA Numerical approximation of a scalar stochastic differential
% equation (SDE) using the Euler-Maruyama method.
%
% On the interval [0,T], solves
%    dX(t) = a(t,X(t))dt + b(t,X(t))dB(t),
%    X(0) = x0.
%
% The Euler-Maruyama scheme is
%   X_{n+1} = X_n + a(t_n,X_n)h + b(t_n,X_n)dW_n.
%
% Strong convergence order: 1/2
% Weak convergence order:   1
%
% INPUTS:
%   a  - drift coefficient, function handle @(t,x)
%   b  - diffusion coefficient, function handle @(t,x)
%   x0 - initial condition
%   T  - final time
%   N  - number of time steps
%   dW - vector of Brownian increments, typically
%        dW(n) ~ N(0,h), h=T/N
%
% OUTPUT:
%   X  - approximation of X(t_n), n=0,...,N

if N <= 0
    error('N must be positive.')
end

if length(dW) ~= N
    error('dW must contain exactly N Brownian increments.')
end

h = T/N; % Uniform time-step size

% Brownian increments dW(n) = B(t_{n+1}) - B(t_n), with dW(n) ~ N(0,h).

X = zeros(N + 1, 1);
X(1) = x0;

for n = 1:N
    tn = (n - 1)*h;
    X(n + 1) = X(n) + a(tn, X(n))*h + b(tn, X(n))*dW(n);
end
end
