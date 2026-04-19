%%

function [t_exec] = CN_PSOR(type, op, Ss, T, K, r, sigma, Nt, Ns)
tic
s_0 = 0; % Inicial value on space
t_0 = 0; % Inicial value on time

ds = (Ss-s_0)/Ns; % difference between points on space
dt = (T-t_0)/Nt; % difference between points on time
s = s_0:ds:Ss;    % points of space
t = t_0:dt:T;    % points of time

S = linspace(0, Ss, Ns+1)';


omega = 1.2; %can be defined between 0 and 2. Particularly 1 and 2?
tol = 1e-6; 
max_iter = 5000;

i = (1:Ns-1)';
    

alpha = 0.25*dt*(sigma^2*i.^2 - r*i);
beta  = -0.5*dt*(sigma^2*i.^2 + r);
gamma = 0.25*dt*(sigma^2*i.^2 + r*i);

% Payoff for put
V = max(K - S, 0);

for n = Nt:-1:1
    
 
 
    
    % Matrices
    A = diag(1 - beta) ...
      + diag(-alpha(2:end), -1) ...
      + diag(-gamma(1:end-1), 1);

    B = diag(1 + beta) ...
      + diag(alpha(2:end), -1) ...
      + diag(gamma(1:end-1), 1);
    
    % RHS
    b = B * V(2:Ns);
    
    % Boundary conditions
    b(1) = b(1) + alpha(1)*(K*exp(-r*(T - (n-1)*dt)));
    
    % Obstacle
    g = max(K - S(2:Ns), 0);
    
    % Solve LCP
    V_inner = PSOR_proj(A, b, g, omega, tol, max_iter);
    
    % Update solution
    V = [K*exp(-r*(T - (n-1)*dt)); V_inner; 0];
end

%At maturity:
V(:, end) = payoff(type, op, S, K);

%early exercise constraint in each time of the loop
V(:, n) = max(V(:, n), payoff(type, op, S, K));


g_full = max(K - S, 0);

[xCrit, ~] = findCriticalPoints(S, V, g_full);

if ~isempty(xCrit)
    Sstar(n) = xCrit(1);   % take first crossing (correct one)
else
    Sstar(n) = NaN;
end

figure;
plot(time, Sstar, 'LineWidth',2);
xlabel('Time');
ylabel('S^*(t)');
title('Exercise Boundary');
grid on;




t_exec = toc;
end
