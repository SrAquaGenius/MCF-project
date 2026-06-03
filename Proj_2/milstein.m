function [t, X] = milstein(a, b, dbdx, x0, T, N, dW)
%MILSTEIN Milstein method for scalar SDE dX = a(t,X)dt + b(t,X)dB + 0.5*b(t,X)*dbdX*((dB)^2 - h), where dbdX is the derivative of b(t,X) with respect to X.

h = T/N;
t = linspace(0, T, N + 1).';


X = zeros(N + 1, 1);
X(1) = x0;

for n = 1:N
    bn = b(t(n), X(n));
    X(n + 1) = X(n) + a(t(n), X(n))*h + bn*dW(n) + 0.5*bn*dbdx(t(n), X(n))*((dW(n))^2 - h);
end
end
