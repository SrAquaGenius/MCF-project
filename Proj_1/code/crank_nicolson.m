function [t_exec] = crank_nicolson(type, op, S, T, K, r, sigma, Ns, Nt)
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

clc, tic

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
U(1,:) = payoff(type, op, s, K);                % Condition t = t_0
U(:,1) = bc(type, op, "left", s_0, t, K, r);    % Ccondition s = s_0
U(:,end) = bc(type, op, "right", s_S, t, K, r); % Ccondition s = S*


% Define Aproximation Formula

i = 2:Ns;      % Indexes on space
S_i = s(i);    % Values of space for each index

alpha = sigma^2/2.*ht.*S_i.^2;
beta = r*ht/2;

a_i = -alpha/2 + beta.*S_i/2;
b_i = 1 + alpha + beta;
c_i = -alpha/2 - beta.*S_i/2;
d_i = 1 - alpha - beta;

% Matrix A (Left side of the system)
A = diag(b_i) + diag(a_i(2:end),-1) + diag(c_i(1:end-1),1);

% Matrix B (Right side of the system)
B = diag(d_i) + diag(-a_i(2:end),-1) + diag(-c_i(1:end-1),1);

for j = 1:Nt
    % vetor da solução no tempo atual (pontos interiores)
    U_now = U(j,2:Ns)';

    % lado direito do sistema
    rhs = B * U_now;

    % adicionar efeitos das fronteiras
    rhs(1)   = rhs(1)   - a_i(1)   * U(j+1,1);
    rhs(end) = rhs(end) - c_i(end) * U(j+1,end);

    U_next = A \ rhs;      % resolver sistema linear
    U(j+1,2:Ns) = U_next'; % guardar solução
end

% Crate surface for visualizing the result
[S,T]=meshgrid(s,t);
figure
surf(S,T,U,'FaceLighting','phong','EdgeColor','none','FaceColor','interp');
xlabel('Stock price S')
ylabel('Time t')
zlabel('Option value V')
title('Black-Scholes Call Price Surface')
colorbar
hold on

% Selecionar 5 valores de t (0%, 25%, 50%, 75%, 100%)
idx_list = round(linspace(1, length(t), 5));

% Linha destacada inicial
idx = idx_list(1);
hLine = plot3(s, t(idx)*ones(size(s)), U(idx,:), ...
              'k', 'LineWidth', 2);

% Botões apenas para esses 5 valores
for k = 1:length(idx_list)
    i = idx_list(k);
    uicontrol('Style','pushbutton',...
        'String',['t = ' sprintf('%.2f', t(i))],...
        'Position',[10 30*k 100 25],...
        'Callback', @(~,~) updateLine(i));
end

% Atualização
function updateLine(i)
    set(hLine, ...
        'YData', t(i)*ones(size(s)), ...
        'ZData', U(i,:));
end

t_exec = toc;
end