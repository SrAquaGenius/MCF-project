%% Initialize Values
clc, clearvars
addpath('Proj_1/code')

% Change the below values as the exercise request!

type = "Eu"; % "Eu" or "Am"
op = "Call";  % "Call" or "Put"

Ss = 300;     % S* stock price boundary
T = 3;       % Maturity time
K = 90;      % Strike price
r = -0.01535;    % Risk free interest rate
sigma = 0.0916; % Market volatility
Ns=3000;      % Numbers of intervals in space
Nt=1000;      % Numbers of intervals in time

%% Crank-Nicolson Method
[tempo_execucao_cn, U] = crank_nicolson(type, op, Ss, T, K, r, sigma, Ns, Nt);
tempo_execucao_cn
%disp(U);

%% Runge-Kutta Method
[tempo_execucao_rk, U] = runge_kutta_4(type, op, Ss, T, K, r, sigma, Ns, Nt);
tempo_execucao_rk
%disp(U);
