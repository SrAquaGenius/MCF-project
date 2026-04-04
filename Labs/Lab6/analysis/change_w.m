%% Analysis of PSOR parameters
clc; clearvars; close all;

% Functions
g1 = @(x)(15/16 + 3*x/8 - 25.*x.^2./16);
g2 = @(x)(-0.1 + (x-0.05).^2 .* (sin(pi.*x)).^2);

g = g2;                   % choose function

% Base parameters
N = 200;
p = 10000;
x = linspace(-1,1,N);

n_w_values = 30;          % Number of w on the interval
eps_w = 1e-3;             % "Epsilon" of w, a > 0 small value
w_values = linspace(eps_w, 2-eps_w, n_w_values);

% Storage
residuals_data = zeros(size(w_values));
iterations_data = zeros(size(w_values));
energy_data = zeros(size(w_values));

%% Loop
tic;
for i = 1:length(w_values)
    w = w_values(i);

    % Ideal: PSOR devolve info extra
    [v, it, res, energy] = PSORextra(g(x(2:end-1)), N-2, p, w);

    iterations_data(i) = it;
    residuals_data(i) = res(end);
    energy_data(i) = energy;
end
toc

%% Plots

figure;
plot(w_values, residuals_data, '-o');
xlabel('\omega'); ylabel('Final residual');
title('Convergence quality');
grid on;

figure;
plot(w_values, iterations_data, '-o');
xlabel('\omega'); ylabel('Iterations');
title('Speed of convergence');
grid on;

figure;
plot(w_values, energy_data, '-o');
xlabel('\omega'); ylabel('Energy');
title('Solution smoothness');
grid on;
