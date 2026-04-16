%% Initialize Values
clc, clearvars
addpath('Proj_1/code')

% Change the below values as the exercise request!

type = "Eu"; % or "Am"
op = "Call";  % or "Put"

Ss = 15;     % S* stock price boundary
T = 1;       % Maturity time
K = 10;      % Strike price
r = 0.06;    % Risk free interest rate
sigma = 0.3; % Market volatility

%% Crank-Nicolson Method
[tempo_execucao_cn, U] = crank_nicolson(type, op, Ss, T, K, r, sigma, 100, 200);
%disp(U);

%% Runge-Kutta Method
[tempo_execucao_rk, U] = runge_kutta_4(type, op, Ss, T, K, r, sigma, 100, 200);
%disp(U);
