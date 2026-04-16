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
U(:,1)   = bc(type, op, "left",  s_0, t, K, r, T); % Condition s = s_0
U(:,end) = bc(type, op, "right", s_S, t, K, r, T); % Condition s = S*

% Construção da matriz espacial (igual ao CN mas sem ht)
i = 2:Ns;
S_i = s(i);

alpha = sigma^2/2 .* S_i.^2;
beta  = r .* S_i;

a_i = alpha - beta / 2;
b_i = -(sigma^2 .* S_i.^2) - r;
c_i = alpha + beta / 2;

A = diag(b_i) + diag(a_i(2:end),-1) + diag(c_i(1:end-1),1);

%% Boundary contribution (generalized)
boundary_vec = @(t_idx) compute_boundary(U, t_idx, a_i, c_i, Ns);

% Função do sistema: dV/dt = A*V + b(t)
F = @(V, t_idx) A*V + boundary_term(t_idx);

% Integração no tempo (BACKWARD!)
for n = Nt:-1:1
    
    Vn = U(n+1,2:Ns)';

    % Time interpolation helpers
    theta = @(a,b,alpha) (1-alpha)*a + alpha*b;

    % Boundary values at different RK stages
    b1 = boundary_vec(n+1);
    b2 = theta(boundary_vec(n+1), boundary_vec(n), 0.5);
    b3 = b2;
    b4 = boundary_vec(n);

    % RK4 stages (fully consistent)
    k1 = A*Vn + b1;
    k2 = A*(Vn + (ht/2)*k1) + b2;
    k3 = A*(Vn + (ht/2)*k2) + b3;
    k4 = A*(Vn + ht*k3)     + b4;

    Vprev = Vn + (ht/6) * (k1 + 2*k2 + 2*k3 + k4);

    %American option (projeção)
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



% Auxiliar function
function b = compute_boundary(U, t_idx, a_i, c_i, Ns)
    b = zeros(Ns-1,1);
    b(1)   = a_i(1)   * U(t_idx,1);
    b(end) = c_i(end) * U(t_idx,end);
end
end