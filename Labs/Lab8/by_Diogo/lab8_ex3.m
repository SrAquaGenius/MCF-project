%% Initialize Values
clc; clearvars; close all;

tic
seed = 0.5; % Must be a value in [0,1]
n = 100000; % Size of pseudo random numbers array

X =@(theta, x) (theta*exp(-theta*x));
Y =@(theta, y) ((y.*exp(-(y.^2)/(2*theta^2)))/(theta^2));
X_inv =@(theta, x) (-log(x)/theta);
Y_inv =@(theta, y) (theta*sqrt(-2*log(y)));

% Select testing function
f = X;
F_inv = X_inv;
theta = 1;

% Apply the inversion method
[realization, exec_time] = inversion_method(n, F_inv, seed, theta);
fprintf('The Invertion Method took %.4f seconds\n', exec_time);

% Plot both histogram of samples and theoretical density
n_sub_int = 100;
histogram(realization, n_sub_int, "Normalization", "pdf");
hold on;

x_axis = linspace(0, max(realization));
y_axis = f(theta, x_axis);
plot(x_axis, y_axis, 'r', 'LineWidth',2);

xlabel('Number of samples');
ylabel('Probability Density');
title('Histogram and Theoretical Density');

legend('Histogram of Samples', ...
       'Theoretical Distribution');

fprintf('The total time elapsed was %.4f seconds\n', toc);
