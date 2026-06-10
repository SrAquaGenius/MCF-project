function [X] = milstein(a, b, dbdx, x0, T, N, dW)
%MILSTEIN Numerical solution of a scalar stochastic differential equation
% (SDE) using the Milstein method.
%
% On the interval [0,T], solves:
%    dX(t) = a(t,X(t))dt + b(t,X(t))dB(t),
%    X(0) = x0.
%
% The Milstein scheme is
%    X_{n+1} = X_n + a(t_n,X_n)h + b(t_n,X_n)dW_n
%            + 0.5*b(t_n,X_n)*b_x(t_n,X_n)*(dW_n^2-h).
%
% Strong convergence order: 1
% Weak convergence order:   1
%
% INPUTS:
%   a    - drift coefficient, function handle @(t,x)
%   b    - diffusion coefficient, function handle @(t,x)
%   dbdx - partial derivative of b with respect to x (b_x(t,x) = ∂b/∂x)
%   x0   - initial condition
%   T    - final time
%   N    - number of time steps
%   dW   - Brownian increments, dW(n) ~ N(0,h)
%
% OUTPUT:
%   X    - approximation of X(t_n), n=0,...,N

if N <= 0
    error('N must be positive.')
end

if length(dW) ~= N
    error('dW must contain exactly N Brownian increments.')
end

h = T/N;  % Uniform time-step size

X = zeros(N + 1, 1);
X(1) = x0;

for n = 1:N
    t_n = (n - 1)*h;

    a_n = a(t_n, X(n));
    b_n = b(t_n, X(n));
    db_n = dbdx(t_n, X(n));

    X(n + 1) = X(n) + a_n*h + b_n*dW(n) + 0.5*b_n*db_n*((dW(n))^2 - h);
end
end
