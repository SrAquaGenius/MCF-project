clc, clearvars

x = linspace(0,20,100);
K = 10;

Payoff_Call = max(x-K,0);
Payoff_Put  = max(K-x,0);

figure
plot(x, Payoff_Call, 'LineWidth', 2)
xlabel('S')
ylabel('Payoff')
title('Call Option Payoff')
grid on

figure
plot(x, Payoff_Put, 'LineWidth', 2)
xlabel('S')
ylabel('Payoff')
title('Put Option Payoff')
grid on