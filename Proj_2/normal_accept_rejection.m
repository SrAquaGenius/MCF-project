function [z, acceptance_rate, last_seed] = normal_accept_rejection(seed, n)
%NORMAL_ACCEPT_REJECTION Standard normal samples using acceptance-rejection.
%   Uses an exponential proposal for the half-normal distribution, the
%   right/positive of a normal distribution.
%      Y ~ Exp(1)
%   , accept with prob exp(-(Y-1)^2/2), then assigns a random sign.

M = 2^31 - 1;
a = 16807;

seed = mod(seed, M);
if seed == 0
    error('For this LCG, seed must not be a multiple of M.');
end

z = zeros(n, 1);
k = 1;
tries = 0;

while k <= n
    seed = mod(a*seed, M);
    u1 = seed/M; % u1 is used to generate the exponential candidate
    seed = mod(a*seed, M);
    u2 = seed/M; % u2 is used to decide whether to accept the candidate
    seed = mod(a*seed, M); 
    u3 = seed/M; % u3 is used to define the sign, as the normal is symmetric

    y = -log(u1);
    tries = tries + 1;

    if u2 <= exp(-0.5*(y - 1)^2) % Condition of acceptance.
        if u3 < 0.5 % Defining the sign.
            z(k) = y;
        else
            z(k) = -y;
        end
        k = k + 1;
    end
end

acceptance_rate = n/tries;
if (nargout == 3); last_seed = seed; end

end