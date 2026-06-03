function [u, last_seed] = uniform_lcg(seed, n)
%UNIFORM_LCG Realizations of U([0,1]) generated with the LCG.

[u, ~, last_seed] = linear_congruential_generator(seed, n);
end
