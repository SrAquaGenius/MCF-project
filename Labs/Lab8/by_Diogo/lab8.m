%% Initialize Values
clc; clearvars; close all;

tic
seed = 0.5; % Must be a value in [0,1]
n = 100000; % Size of pseudo random numbers array

X =@(theta, x) (theta*exp(-theta*x));
Y =@(theta, y) ((y.*exp(-(y.^2)/(2*theta^2)))/(theta^2));
X_inv =@(theta, x) (-log(x)/theta);
Y_inv =@(theta, y) (theta*sqrt(-2*log(y)));

%% 1) Apply Linear Congruential Generator
[u, exec_time] = linear_congruential_generator(n, seed)
toc

%% 2) Apply the Inversion Method
F = Y_inv;
theta = 0.1;
[f, exec_time] = inversion_method(n, F, theta, seed)
toc

%% 3) Plot the results on an histogram
f = Y;
F_inv = Y_inv;
theta = 1;

[realization, ] = inversion_method(n, F_inv, theta, seed);
n_sub_int = 100;
histogram(realization, n_sub_int, "Normalization", "pdf");
hold on;

x_axis = linspace(0, max(realization));
y_axis = f(theta, x_axis);
plot(x_axis, y_axis, 'r', 'LineWidth',2);

toc
