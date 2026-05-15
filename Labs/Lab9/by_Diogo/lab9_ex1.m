%% Initialize Values
clc; clearvars; close all;

tic
seed = 0.1; % Must be a value in [0,1]

g =@(x)(exp(x));
I = integral(g, 0, 1);

N = [100, 500, 1000, 5000, 10000, 50000, 100000];
size_n = length(N);
m = zeros(size_n, 1);
e = zeros(size_n, 1);
t = zeros(size_n, 1);

for i = 1:size_n
    [m(i), t(i)] = monte_carlo(N(i), seed, g);
    e(i) = abs(I - m(i));
end

T = table(N', m, e, t, ...
    'VariableNames', {'N Samples', 'Monte Carlo', 'Error', 'Elapsed Time'});

fig = figure('Name', 'Monte Carlo Results');

uitable(fig, ...
    'Data', T, ...
    'ColumnName', T.Properties.VariableNames, ...
    'Units', 'Normalized', ...
    'Position', [0 0 1 1]);

fprintf('The total time elapsed was %.4f seconds\n', toc);
