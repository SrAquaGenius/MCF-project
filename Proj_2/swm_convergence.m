function results = swm_convergence(S0, mu, sigma, T, h_values, nSim, seed, chunkSize)
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

        [Z, seed] = normal_accept_rejection(seed, m*N);
        Z = reshape(Z, m, N);

        SEM = S0*ones(m, 1);
        SMil = S0*ones(m, 1);
        B = zeros(m, 1);

        for n = 1:N
            dW = sqrt(h)*Z(:, n);

            B = B + dW;

            SEM = SEM + mu*SEM*h + sigma*SEM.*dW;

            SMil = SMil + mu*SMil*h + sigma*SMil.*dW ...
                + 0.5*sigma^2*SMil.*(dW.^2 - h);
        end

        Sexact = S0*exp((mu - 0.5*sigma^2)*T + sigma*B);

        sumAbsEM = sumAbsEM + sum(abs(SEM - Sexact));
        sumAbsMil = sumAbsMil + sum(abs(SMil - Sexact));
        sumDiffEM = sumDiffEM + sum(SEM - Sexact);
        sumDiffMil = sumDiffMil + sum(SMil - Sexact);

        done = done + m;
    end

    strongEM(ih) = sumAbsEM/nSim;
    strongMil(ih) = sumAbsMil/nSim;
    weakEM(ih) = abs(sumDiffEM/nSim);
    weakMil(ih) = abs(sumDiffMil/nSim);
end

orderStrongEM = polyfit(log(h_values(:)), log(strongEM), 1);
orderStrongMil = polyfit(log(h_values(:)), log(strongMil), 1);
orderWeakEM = polyfit(log(h_values(:)), log(weakEM), 1);
orderWeakMil = polyfit(log(h_values(:)), log(weakMil), 1);

results = table(h_values(:), strongEM, strongMil, weakEM, weakMil, ...
    'VariableNames', {'h', 'StrongEM', 'StrongMilstein', 'WeakEM', 'WeakMilstein'});

fprintf('Estimated strong order EM: %.3f\n', orderStrongEM(1));
fprintf('Estimated strong order Milstein: %.3f\n', orderStrongMil(1));
fprintf('Estimated weak order EM: %.3f\n', orderWeakEM(1));
fprintf('Estimated weak order Milstein: %.3f\n', orderWeakMil(1));
end
