%% Initialize Values
clc; clearvars; close all;

tic
seed = 0.5; % Must be a value in [0,1]
n = 100000; % Size of pseudo random numbers array

% Apply Linear Congruential Generator
[u, exec_time] = linear_congruential_generator(n, seed);
fprintf('The Linear Congruential Generator took %.4f seconds\n', exec_time);

% Plot the results
figure;
histogram(u, 50, "Normalization", "pdf");
hold on;

% Theoretical uniform density on [0,1]
x = linspace(0, 1, n);
y = ones(size(x));

plot(x, y, 'r', 'LineWidth', 2);

xlabel('x');
ylabel('Probability Density');
title('Linear Congruential Generator - Uniform Distribution');

legend('Generated samples', 'Uniform PDF');
grid on;

fprintf('The total time elapsed was %.4f seconds\n', toc);
