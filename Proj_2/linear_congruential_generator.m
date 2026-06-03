function [res, t, last_seed] = linear_congruential_generator(seed, n)
%LINEAR_CONGRUENTIAL_GENERATOR Park-Miller LCG on (0,1).
%   [res,t,last_seed] = linear_congruential_generator(seed,n) returns n
%   pseudo-random values in (0,1), the elapsed time, and the last integer
%   state. Use last_seed to continue the same stream.

tic

M = 2^31 - 1;
a = 16807;
b = 0;



m = zeros(n, 1);
m(1) = seed;

for i = 2:n
    m(i) = mod(a*m(i-1) + b, M);
end

res = m/M;
last_seed = m(end);
t = toc;
end
