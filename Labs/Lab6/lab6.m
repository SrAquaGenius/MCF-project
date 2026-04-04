%% Initialize Values
clc; clearvars; close all; tic
addpath('Labs/Lab6/')

g1 =@ (x)(15/16+3*x/8-25.*x.^2./16);
g2 =@ (x)(-0.1+(x-0.05).^2.*(sin(pi.*x)).^2);

N = 300;
p = 20000; % Precision Value = number of iterations
w = 1.8;   % Relaxation omega value in ]0, 2[
x = linspace(-1,1,N);

g = g2;    % function selector

%% Apply it for a choosen g function
u = PSOR(g(x(2:(end-1))), N-2, p, w);
[xCrit, yCrit] = findCriticalPoints(x, u, g(x));
plotPSOR(g, x, u, xCrit, yCrit);
toc
