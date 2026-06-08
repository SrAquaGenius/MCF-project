function [u, last_seed] = uniform_lcg(seed, n)
%UNIFORM_LCG Realizations of U([0,1]) generated with the LCG.

if (nargout == 2)
    [u, last_seed] = linear_congruential_generator(seed, n);
else
    u = linear_congruential_generator(seed, n);
end
end
