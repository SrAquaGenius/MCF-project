function [t_exec, result] = crank_nicolson(type, op, S, T, K, r, sigma, Ns, Nt)
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
U(end,:) = payoff(type, op, s, K);                 % Condition t = t_0
U(:,1) = bc(type, op, "left", s_0, t, K, r, T);    % Condition s = s_0
U(:,end) = bc(type, op, "right", s_S, t, K, r, T); % Condition s = S*


% Define Aproximation Formula
i = 2:Ns;      % Indexes on space
S_i = s(i);    % Values of space for each index

% Coeficiente do método de Crank-Nicolson
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

%disp(U)

% Crate surface for visualizing the result
[S,T]=meshgrid(s,t);
figure
surf(S,T,U,'FaceLighting','phong','EdgeColor','none','FaceColor','interp');
if strcmpi(op,'put')
    view(45, 25)    % eixo do espaço à direita
else
    view(-45, 25)   % eixo do espaço à esquerda
end
xlabel('Stock price S')
ylabel('Time t')
zlabel('Option value V')

if strcmpi(op,'put'), Op = 'Put'; else Op = 'Call'; end
if strcmpi(type,'am'), Type = 'American'; else Type = 'European'; end
title(['Black-Scholes ' Op ' Price Surface on an ' Type ' Option'])
colorbar
hold on

% ------------------ Criar grelha na malha de superfície ------------------
t_vals = linspace(t(1), t(end), 10);
s_vals = linspace(s(1), s(end), 6);

% Linhas paralelas ao eixo S (tempo fixo)
for tt = t_vals

    [~, idx_t] = min(abs(t - tt)); % encontrar índice mais próximo

    plot3(s, ...
          t(idx_t)*ones(size(s)), ...
          U(idx_t,:), ...
          'k-', 'LineWidth', 0.8);
end

% Linhas paralelas ao eixo t (espaço fixo)
for ss = s_vals
    
    [~, idx_s] = min(abs(s - ss)); % encontrar índice mais próximo
    
    plot3(s(idx_s)*ones(size(t)), ...
          t, ...
          U(:,idx_s), ...
          'k-', 'LineWidth', 0.8);
end

% -------------------- Selecionador de linhas temporais -------------------
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
result = U;
end