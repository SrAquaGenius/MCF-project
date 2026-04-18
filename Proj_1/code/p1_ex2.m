%% Initialize Values
clc, clearvars
%addpath('Proj_1/code')

% Change the below values as the exercise request!

type = "AM"; % or "Eu"
op = "Put";  % or "Call"

Ss = 15;     % S* stock price boundary
T = 1;       % Maturity time
K = 10;      % Strike price
r = 0.06;    % Risk free interest rate
sigma = 0.3; % Market volatility
Ns=1500;      % Numbers of intervals in space
Nt=100;      % Numbers of intervals in time
%% Crank-Nicholson & PSOR
%define only in the continuation region. 
%find the critical points.
%in the stopping region, we have it defined as the payoff function.
[tempo_execucao_cn_psor] = CN_PSOR(type, op, Ss, T, K, r, sigma, Nt, Ns);
tempo_execucao_cn_psor