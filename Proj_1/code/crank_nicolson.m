function [t_exec, U] = crank_nicolson(type, op, S, T, K, r, sigma, Ns, Nt)
% Crank-Nicolson:
% - type: Type of Option: "Am" American or "Eu" European
% - op: Operation type: "Put" or "Call"
% - S: S* stock price boundary
% - T: Maturity time
% - K: Strike price
% - r: Risk free interest rate
% - sigma: Market Volatility
% - Ns: Number of subdivisions on space
% - Nt: Number of subdivisions on time

tic

s_0 = 0; % Inicial value on space
s_S = S; % Final value on space
t_0 = 0; % Inicial value on time
t_T = T; % Final value on time = Maturity time

hs = (S-s_0)/Ns; % difference between points on space
ht = (T-t_0)/Nt; % difference between points on time
s = s_0:hs:S;    % points of space
t = t_0:ht:T;    % points of time

% Initialize matrices for storing results (lines -> time; columns -> space)
U = zeros(Nt+1, Ns+1); % Matrix that stores the final results
U(end,:) = payoff(type, op, s, K);                 % Condition t = t_0
U(:,1) = bc(type, op, "left", s_0, t, K, r, T);    % Condition s = s_0
U(:,end) = bc(type, op, "right", s_S, t, K, r, T); % Condition s = S*

% Define Aproximation Formula
i = 2:Ns;      % Indexes on space
S_i = s(i);    % Values of space for each index

% Coeficientes do método de Crank-Nicolson
alpha = sigma^2/2.*ht.*S_i.^2;
beta = r*ht/2;

a_i = -alpha/2 + beta.*S_i/2;
b_i = 1 + alpha + beta;
c_i = -alpha/2 - beta.*S_i/2;
d_i = 1 - alpha - beta;

% Matrix A (lado esquerdo do sistema)
A = diag(b_i) + diag(a_i(2:end),-1) + diag(c_i(1:end-1),1);

% Matrix B (lado direito do sistema)
B = diag(d_i) + diag(-a_i(2:end),-1) + diag(-c_i(1:end-1),1);


% A = B*u => u_new = A \ (B*u_old)
for j = Nt:-1:1
    % vetor da solução no tempo atual (pontos interiores)
    U_now = U(j+1,2:Ns)';

    % lado direito do sistema
    rhs = B * U_now;

    % adicionar efeitos das fronteiras
    rhs(1) = rhs(1) - a_i(1) * U(j+1,1) - a_i(1) * U(j,1);
    rhs(end) = rhs(end) - c_i(end) * U(j+1,end) - c_i(end) * U(j,end);

    % resolver sistema linear e guardar solução
    U(j,2:Ns) = (A \ rhs)';
end

plot_surface(s, t, U, type, op);
setup_time_selector(s, t, U);

t_exec = toc;
end