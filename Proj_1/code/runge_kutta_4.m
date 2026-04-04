function [t_exec] = runge_kutta_4(type, op, S, T, K, r, sigma, Ns, Nt)
% Runge-Kutta 4 applied to Black-Scholes by the method of lines

clc, tic

s_0 = 0; % Inicial value on space
s_S = S; % Final value on space
t_0 = 0; % Inicial value on time
t_T = T; % Final value on time = Maturity time

hs = (S-s_0)/Ns; % difference between points on space
ht = (T-t_0)/Nt; % difference between points on time
s = s_0:hs:S;    % points of space
t = t_0:ht:T;    % points of time

% Initialize matrix that stores the final results
U = zeros(Nt+1, Ns+1);
U(end,:) = payoff(type, op, s, K);              % Condition t = t_0
U(:,1)   = bc(type, op, "left",  s_0, t, K, r); % Condition s = s_0
U(:,end) = bc(type, op, "right", s_S, t, K, r); % Condition s = S*

% Construção da matriz espacial (igual ao CN mas sem ht)
i = 2:Ns;
S_i = s(i);

alpha = sigma^2/2 .* S_i.^2;
beta  = r .* S_i;

a_i = alpha/(hs^2) - beta/(2*hs);
b_i = -2*alpha/(hs^2) - r;
c_i = alpha/(hs^2) + beta/(2*hs);

A = diag(b_i) + diag(a_i(2:end),-1) + diag(c_i(1:end-1),1);

% Função do sistema: dV/dt = A*V + b(t)
F = @(V, t_idx) A*V + boundary_term(t_idx);

% Integração no tempo (BACKWARD!)
for n = Nt:-1:1
    
    Vn = U(n+1,2:Ns)';

    % RK4
    k1 = F(Vn, n+1);
    k2 = F(Vn + ht/2*k1, n+1);
    k3 = F(Vn + ht/2*k2, n+1);
    k4 = F(Vn + ht*k3,   n+1);

    Vprev = Vn + ht/6 * (k1 + 2*k2 + 2*k3 + k4);

    % American option (projeção)
    if type == "Am"
        payoff_vec = payoff(type, op, s(2:Ns), K)';
        Vprev = max(Vprev, payoff_vec);
    end

    U(n,2:Ns) = Vprev';
end

% Plotting
[Sg,Tg] = meshgrid(s,t);

figure
surf(Sg,Tg,U,'FaceLighting','phong','EdgeColor','none','FaceColor','interp');
xlabel('Stock price S')
ylabel('Time t')
zlabel('Option value V')
title('Black-Scholes RK4 Surface')
colorbar

t_exec = toc;

% ===== função auxiliar =====
function b = boundary_term(t_idx)
    b = zeros(Ns-1,1);
    
    % contribuição das fronteiras
    b(1)   = a_i(1)   * U(t_idx,1);
    b(end) = c_i(end) * U(t_idx,end);
end

end