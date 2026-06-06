function [X] = euler_maruyama(a, b, x0, T, N, dW)
% Euler-Maruyama method for dX = a(t,X)dt + b(t,X)dB.

h = T/N;
t = linspace(0, T, N + 1).'; % define the time grid as a column vector.

%the Brownian motion increments of h, then it follows a Normal distribution: N(0,h).

X = zeros(N + 1, 1);
X(1) = x0;

for n = 1:N
    X(n + 1) = X(n) + a(t(n), X(n))*h + b(t(n), X(n))*dW(n);
end
end
