clear
clc
close all

%% 1. Random numbers  
seed = 1;
% assuming always the 10^6 samples
[u, seed] = uniform_lcg(seed, 10000);
[z, ~, accRate] = normal_accept_rejection(seed, 10000);

figure
histogram(u, 40, 'Normalization', 'pdf')
title('U([0,1]) with LCG')

figure
histogram(z, 60, 'Normalization', 'pdf')
hold on
xx = linspace(-4, 4, 400);
plot(xx, exp(-0.5*xx.^2)/sqrt(2*pi), 'LineWidth', 1.5)
title(sprintf('N(0,1) by acceptance-rejection, rate %.3f', accRate))

figure
pointsHaltonDemo = halton_nodes(1000, [2 3]);
plot(pointsHaltonDemo(:, 1), pointsHaltonDemo(:, 2), '.')
axis square
title('Halton nodes in [0,1]^2')

%% 2. Area curve alpha -> |A_alpha|
Narea = 1e6;
alpha = linspace(-1, 1, 1000).';
seed = 1;
[u2, ~] = uniform_lcg(seed, 2*Narea);
pointsMC = reshape(u2, Narea, 2);
pointsQMC = halton_nodes(Narea, [2 3]);

areaMC = estimate_area_curve(pointsMC, alpha);
areaQMC = estimate_area_curve(pointsQMC, alpha);

figure
plot(alpha, areaMC, 'b.', 'MarkerSize', 6)
hold on
plot(alpha, areaQMC, 'r.', 'MarkerSize', 3)
legend('MC', 'QMC', 'Location', 'northwest')
xlabel('\alpha')
ylabel('|A_\alpha|')
title('Estimated area curve')
grid on

figure
plot(alpha, areaQMC - areaMC, 'k', 'LineWidth', 1)
xlabel('\alpha')
ylabel('|A_\alpha|_{QMC} - |A_\alpha|_{MC}')
title('Difference between QMC and MC area estimates')
grid on



%possible graphic that evaluates the absolute error between the real values
%of the function and the values obtained for each method
%% 3 and 4(a). Euler-Maruyama and Milstein for SDE
S0 = 1; %The convergence orders should not change. The paths and errors would scale with S0, but the theoretical order of convergence remains the same.
T = 1;
mu = 0.6;
sigma = 0.25;
h = 0.001;
N = round(T/h);
seed = 1;

[Z, seed] = normal_accept_rejection(seed, N);
dW = sqrt(h)*Z;
t = linspace(0, T, N + 1).';
B = [0; cumsum(dW)]; %cumulative sums of the increments of the brownian motion to build a brownian motion (create the full brownian path).

drift = @(t, x) mu*x;
diffusion = @(t, x) sigma*x;
diffusionDx = @(t, x) sigma;

[~, SEM] = euler_maruyama(drift, diffusion, S0, T, N, dW);
[~, SMil] = milstein(drift, diffusion, diffusionDx, S0, T, N, dW);
Sexact = S0*exp((mu - 0.5*sigma^2)*t + sigma*B);

figure
plot(t, Sexact, 'k', 'LineWidth', 1.5)
hold on
plot(t, SEM, '--', 'LineWidth', 1.2)
plot(t, SMil, ':', 'LineWidth', 1.5)
legend('Exact', 'Euler-Maruyama', 'Milstein', 'Location', 'northwest')
xlabel('t')
ylabel('S(t)')
title('Exact and numerical paths')
grid on

%here, we are only using one realization of the brownian motion.

errorEM = abs(Sexact - SEM);
errorMil = abs(Sexact - SMil);

figure
plot(t, errorEM, 'b', 'LineWidth', 1.2)
hold on
plot(t, errorMil, 'r', 'LineWidth', 1.2)
legend('Euler-Maruyama error', 'Milstein error', 'Location', 'northwest')
xlabel('t')
ylabel('Absolute error')
title('Absolute errors along one Brownian path')
grid on

%% 4(b). Convergence study
hValues = 0.005*(1/2).^(0:3);
nSim = 1e6;


convResults = swm_convergence(S0, mu, sigma, T, hValues, nSim, 2026, 10000);
disp(convResults)

figure
loglog(convResults.h, convResults.StrongEM, 'o-', ...
    convResults.h, convResults.StrongMilstein, 's-', ...
    convResults.h, convResults.WeakEM, 'o--', ...
    convResults.h, convResults.WeakMilstein, 's--')
grid on
xlabel('h')
ylabel('error')
legend('Strong EM', 'Strong Milstein', 'Weak EM', 'Weak Milstein', ...
    'Location', 'northwest')
title('Strong and weak convergence')
