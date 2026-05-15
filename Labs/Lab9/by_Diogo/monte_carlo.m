function [res, t] = monte_carlo(n, s, g)
%MONTE_CARLO Summary of this function goes here

tic;

res = (sum(g(rand(1, n))))/n;

t = toc;
end