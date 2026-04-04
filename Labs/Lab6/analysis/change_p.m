%% Analysis of PSOR parameter p
clc; clearvars; close all;

% Functions
g1 = @(x)(15/16 + 3*x/8 - 25.*x.^2./16);
g2 = @(x)(-0.1 + (x-0.05).^2 .* (sin(pi.*x)).^2);

g = g2;                   % choose function

% Base parameters
N = 200;
w = 1.8;
x = linspace(-1,1,N);

% Range of p values
p_values = round(linspace(10, 2000, 30));   % podes ajustar isto

% Storage
residuals_data = zeros(size(p_values));
iterations_data = zeros(size(p_values));
energy_data = zeros(size(p_values));

%% Loop
tic;
for i = 1:length(p_values)
    p = p_values(i);

    [v, it, res, energy] = PSORextra(g(x(2:end-1)), N-2, p, w);

    iterations_data(i) = it;
    residuals_data(i) = res(end);
    energy_data(i) = energy;
end
toc

%% Plots

figure;
plot(p_values, residuals_data, '-o');
xlabel('p'); ylabel('Final residual');
title('Convergence quality vs p');
grid on;

figure;
plot(p_values, iterations_data, '-o');
xlabel('p'); ylabel('Iterations used');
title('Effective iterations vs p');
grid on;

figure;
plot(p_values, energy_data, '-o');
xlabel('p'); ylabel('Energy');
title('Solution smoothness vs p');
grid on;