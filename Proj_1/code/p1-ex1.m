%% Initialize Values
clc, clearvars

% Change the below values as the exercise request!

type = "Eu"; % or "Am"
op = "Call";  % or "Call"

Ss = 15;     % S* stock price boundary
T = 1;       % Maturity time
K = 10;      % Strike price
r = 0.06;    % Risk free interest rate
sigma = 0.3; % Market volatility

%% Crank-Nicolson Method
t_exec = crank_nicolson(type, op, Ss, T, K, r, sigma, 100, 200)

%% Runge-Kutta Method
t_exec = runge_kutta_4(type, op, Ss, T, K, r, sigma, 100, 200)
