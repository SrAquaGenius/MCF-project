%% Initialize Values
clc; clearvars; close all;

tic
seed = 0.1; % Must be a value in [0,1]
n = 100000; % Size of pseudo random numbers array
mu = 1;
sigma = 1.5;

% Generate samples using Box-Muller
[z, t] = nrm_distribution(n, seed, mu, sigma);
fprintf('The execution time was %.6f seconds\n', t);

% Plot histogram vs theoretical normal density
figure;

histogram(z, 100, "Normalization", "pdf");
hold on;

x_axis = linspace(min(z), max(z), 1000);
y_axis = (1/(sigma*sqrt(2*pi))) * exp(-((x_axis - mu).^2) / (2*sigma^2));

plot(x_axis, y_axis, 'r', 'LineWidth', 2);

xlabel('x');
ylabel('Probability Density');
title('Box-Muller Method - Standard Normal Distribution');

legend('Generated samples', 'N(0,1) PDF');
grid on;

fprintf('The total time elapsed was %.4f seconds\n', toc);
