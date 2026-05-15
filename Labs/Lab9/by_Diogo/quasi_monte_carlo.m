function [res, t] = quasi_monte_carlo(n, s, g)
%QUASI_MONTE_CARLO Summary of this function goes here

tic;

haltonset

t = toc;
end