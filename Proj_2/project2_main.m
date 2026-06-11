clear, clc, close all

%% 1. Random numbers  
seed = 1;          % we choose this value
n_samples = 10000; % assumed 10^6 samples

[u, seed] = uniform_lcg(seed, n_samples);
[z, accRate] = normal_accept_rejection(seed, n_samples);

% Figure for the plot of the uniform distribution
figure
histogram(u, 40, 'Normalization', 'pdf')
hold on
yline(1, 'LineWidth', 1.5, 'Color','r')
title('U([0,1]) with LCG')

% Figure for the plot of the normal distribution
figure
histogram(z, 60, 'Normalization', 'pdf')
hold on
xx = linspace(-4, 4, 400);
plot(xx, exp(-0.5*xx.^2)/sqrt(2*pi), 'LineWidth', 1.5)
title(sprintf('N(0,1) by acceptance-rejection, rate %.3f', accRate))

% Figure for the plot of the halton points
figure
n_points = 1000;
pointsHaltonDemo = halton_nodes(n_points, [2 3]);
plot(pointsHaltonDemo(:, 1), pointsHaltonDemo(:, 2), '.')
axis square
title('Halton nodes in [0,1]^2')

%% 2. Area curve alpha -> |A_alpha|
n_samples = 1e6;
alpha = linspace(-1, 1, 1000).';
seed = 1;                                   % We choose this value

f =@ (x,y) sin(10*(x.^2 - sin(3*y)));       % Evaluating Curve

u = uniform_lcg(seed, 2*n_samples);         % Uniform distribution
pointsMC = reshape(u, n_samples, 2);        % Points of Monte Carlo
pointsQMC = halton_nodes(n_samples, [2 3]); % Points of Quasi Monte Carlo

fMC = f(pointsMC(:,1),  pointsMC(:,2));     % Monte Carlo after f
fQMC = f(pointsQMC(:,1),  pointsQMC(:,2));  % Quasi Monte Carlo after f

areaMC  = mean(fMC  < alpha', 1);           % Area for Monte Carlo
areaQMC = mean(fQMC < alpha', 1);           % Area for Quasi Monte Carlo

% Figure for the plot of the function, on 3D
[x, y] = meshgrid(linspace(0,1,200), linspace(0,1,200));
z = f(x,y);
figure
surf(x, y, z)
shading interp
xlabel('x')
ylabel('y')
zlabel('f(x,y)')
title('Surface plot of f(x,y)')

% Figure for the plot of the estimated areas
figure
plot(alpha, areaMC, 'b.', 'MarkerSize', 6)
hold on
plot(alpha, areaQMC, 'r.', 'MarkerSize', 3)
legend('MC', 'QMC', 'Location', 'northwest')
xlabel('\alpha')
ylabel('|A_\alpha|')
title('Estimated area curve')
grid on

% Figure for the plot of the difference between estimated areas
figure
plot(alpha, areaQMC - areaMC, 'k', 'LineWidth', 1)
xlabel('\alpha')
ylabel('|A_\alpha|_{QMC} - |A_\alpha|_{MC}')
title('Difference between QMC and MC area estimates')
grid on

%% 3 and 4(a). Euler-Maruyama and Mil for SDE
%The convergence orders should not change. The paths and errors would scale
% with S0, but the theoretical order of convergence remains the same.

T = 1;
mu = 0.6;
sigma = 0.25;
h = 0.001;
N = round(T/h);

seed = 1;
S_0 = 1;

t = linspace(0, T, N + 1).';

% Get a single realization of the brownian motion
Z = normal_accept_rejection(seed, N);
if length(Z) ~= N
    error('Mismatch between Brownian increments and N.')
end
dW = sqrt(h)*Z;
% Cumulative sums of the increments to create the full brownian path.
B = [0; cumsum(dW)];

% Get the drift and the diffusion and applying it to the methods
drift = @(t, x) mu*x;
diffusion = @(t, x) sigma*x;
diffusionDx = @(t, x) sigma;

[SEM] = euler_maruyama(drift, diffusion, S_0, T, N, dW);
[SMil] = milstein(drift, diffusion, diffusionDx, S_0, T, N, dW);

% Exact solution (geometric brownian motion)
S_exact = S_0*exp((mu - 0.5*sigma^2)*t + sigma*B);

% Figure of the plot of the Brownian paths trajectories
figure
plot(t, S_exact, 'k', 'LineWidth', 1.5)
hold on
plot(t, SEM, '--', 'LineWidth', 1.2)
plot(t, SMil, ':', 'LineWidth', 1.5)
legend('Exact', 'Euler-Maruyama', 'Mil', 'Location', 'northwest')
xlabel('t')
ylabel('S(t)')
title('Exact and numerical paths')
grid on

% Figure of the plot of the errors of the Brownian paths trajectories
% Here, we are only using one realization of the brownian motion.
errorEM = abs(S_exact - SEM);
errorMil = abs(S_exact - SMil);

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

T = 1;
mu = 0.6;
sigma = 0.25;

S_0 = 1;

i_vals = (0:3);
h_vals = 0.005*(1/2).^i_vals;
n_sim = 1e6; % Number of simulations

[results, orders] = swm_convergence(S_0, mu, sigma, T, h_vals, n_sim, 1, 10000);
disp(results)
disp(orders)

figure
loglog(results.h, results.StrongEM, 'o-', ...
    results.h, results.StrongMil, 's-', ...
    results.h, results.WeakEM, 'o--', ...
    results.h, results.WeakMil, 's--')
grid on
xlabel('h')
ylabel('error')
legend('Strong EM', 'Strong Milstein', 'Weak EM', 'Weak Milstein', ...
    'Location', 'northwest')
title('Strong and weak convergence')
