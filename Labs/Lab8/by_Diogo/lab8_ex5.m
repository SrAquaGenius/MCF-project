%% Initialize Values
clc; clearvars; close all;

tic
seed = 0.1; % Must be a value in [0,1]
n = 100000; % Size of pseudo random numbers array

% Generate samples using Box-Muller
[x, y, t] = box_muller_method(n, seed);
fprintf('The execution time was %.6f seconds\n', t);

% Plot histogram vs theoretical normal density
z = [x(:); y(:)];
figure;

histogram(x, 100, "Normalization", "pdf");
hold on;

x_axis = linspace(min(x), max(x), 1000);
y_axis = (1/sqrt(2*pi)) * exp(-x_axis.^2 / 2);

plot(x_axis, y_axis, 'r', 'LineWidth', 2);

xlabel('x');
ylabel('Probability Density');
title('Box-Muller Method - Standard Normal Distribution');

legend('Generated samples', 'N(0,1) PDF');
grid on;

fprintf('The total time elapsed was %.4f seconds\n', toc);
