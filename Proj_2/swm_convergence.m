function [results, orders] = swm_convergence(S0, mu, sigma, T, h_values, nSim, seed, chunkSize)
% we simulate in smaller groups (chunks), since we have 1e6 simulations for differents values of h.
%   Uses vectorized chunks to avoid storing all paths.

if nargin < 8
    chunkSize = 10000;
end

nH = numel(h_values);
strongEM = zeros(nH, 1);
strongMil = zeros(nH, 1);
weakEM = zeros(nH, 1);
weakMil = zeros(nH, 1);

for ih = 1:nH
    h = h_values(ih);
    N = round(T/h);
    h = T/N;

    sumAbsEM = 0;
    sumAbsMil = 0;
    sumDiffEM = 0;
    sumDiffMil = 0;
    done = 0;

    while done < nSim
        m = min(chunkSize, nSim - done);

        [Z, ~ ,seed] = normal_accept_rejection(seed, m*N);
        Z = reshape(Z, m, N);

        S_EM = S0*ones(m, 1);
        S_Mil = S0*ones(m, 1);
        B = zeros(m, 1);

        for n = 1:N
            dW = sqrt(h)*Z(:, n);

            B = B + dW;

            S_EM = S_EM + mu*S_EM*h + sigma*S_EM.*dW;

            S_Mil = S_Mil + mu*S_Mil*h + sigma*S_Mil.*dW ...
                + 0.5*sigma^2*S_Mil.*(dW.^2 - h);
        end

        S_exact = S0*exp((mu - 0.5*sigma^2)*T + sigma*B);

        sumAbsEM = sumAbsEM + sum(abs(S_EM - S_exact));
        sumAbsMil = sumAbsMil + sum(abs(S_Mil - S_exact));
        sumDiffEM = sumDiffEM + sum(S_EM - S_exact);
        sumDiffMil = sumDiffMil + sum(S_Mil - S_exact);

        done = done + m;
    end

    strongEM(ih) = sumAbsEM/nSim;
    strongMil(ih) = sumAbsMil/nSim;
    weakEM(ih) = abs(sumDiffEM/nSim);
    weakMil(ih) = abs(sumDiffMil/nSim);
end

results = table(h_values(:), strongEM, strongMil, weakEM, weakMil, ...
    'VariableNames', {'h', 'StrongEM', 'StrongMil', 'WeakEM', 'WeakMil'});

h_values = h_values(:);

pStrongEM_aux    = polyfit(log(h_values(:)), log(strongEM), 1);
pStrongEM        = pStrongEM_aux(1);
pStrongEM_local  = mean(log(strongEM(1:end-1)./strongEM(2:end)) ./ ...
                  log(h_values(1:end-1)./h_values(2:end)));
pStrongMil_aux   = polyfit(log(h_values(:)), log(strongMil), 1);
pStrongMil       = pStrongMil_aux(1);
pStrongMil_local = mean(log(strongMil(1:end-1)./strongMil(2:end)) ./ ...
                  log(h_values(1:end-1)./h_values(2:end)));
pWeakEM_aux      = polyfit(log(h_values(:)), log(weakEM), 1);
pWeakEM          = pWeakEM_aux(1);
pWeakEM_local    = mean(log(weakEM(1:end-1)./weakEM(2:end)) ./ ...
                  log(h_values(1:end-1)./h_values(2:end)));
pWeakMil_aux     = polyfit(log(h_values(:)), log(weakMil), 1);
pWeakMil         = pWeakMil_aux(1);
pWeakMil_local   = mean(log(weakMil(1:end-1)./weakMil(2:end)) ./ ...
                  log(h_values(1:end-1)./h_values(2:end)));

orders = table(...
    pStrongEM, pStrongEM_local, pStrongMil, pStrongMil_local, ...
    pWeakEM, pWeakEM_local, pWeakMil, pWeakMil_local, ...
    'VariableNames', {...
    'pStrongEM', 'pStrongEM_local', 'pStrongMil', 'pStrongMil_local', ...
    'pWeakEM', 'pWeakEM_local', 'pWeakMil', 'pWeakMil_local'});
end
