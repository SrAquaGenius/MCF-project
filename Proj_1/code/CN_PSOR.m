function [t_exec] = CN_PSOR(type, op, Ss, T, K, r, sigma, Nt, Ns)
tic

s_0 = 0;
t_0 = 0;
dt  = (T  - t_0) / Nt;
ds  = (Ss - s_0) / Ns;

%  Grids: ALWAYS column vectors 
t = (t_0 : dt : T)';            % (Nt+1) x 1
S = (s_0 : ds : Ss)';           % (Ns+1) x 1

m   = Ns - 1;                   % number of interior nodes
idx = (1:m)';                   % column vector of interior indices

%  CN coefficients (column vectors, length m) 
alpha = 0.25*dt*(sigma^2*idx.^2 - r*idx);
beta  =  0.5*dt*(sigma^2*idx.^2 + r);      % note: positive; enters as (1+beta)
gamma = 0.25*dt*(sigma^2*idx.^2 + r*idx);

%  System matrices (m x m)
A = diag(1 + beta)          ...
  - diag(gamma(1:m-1), +1)  ...
  - diag(alpha(2:m),   -1);

B = diag(1 - beta)          ...
  + diag(gamma(1:m-1), +1)  ...
  + diag(alpha(2:m),   -1);

% ── PSOR parameters 
omega    = 2 / (1 + sin(pi / (m+1)));
tol      = 1e-6;
max_iter = 5000;

% ── Initial condition V(S, T) = payoff, column vector (Ns+1) x 1 
V = max(K - S, 0);

Sstar       = NaN(Nt, 1);
V_old_inner = V(2:Ns);          % (m x 1) warm-start for PSOR
V_surf          = zeros(Ns+1, Nt+1); % store full surface
V_surf(:, Nt+1) = V;                 % payoff at t=T

% Backward time-stepping 
for n = Nt:-1:1

    % RHS: (m x m) * (m x 1) = (m x 1)
    b = B * V(2:Ns);

    % Lower BC: American put at S=0 is always worth K (exercise immediately,
    % no need to discount — the holder receives K now, not at maturity).
    % Both CN levels use K (it is constant in time for an American option).
    b(1) = b(1) + alpha(1) * (K + K);   % = 2*alpha(1)*K

    % Obstacle: (m x 1) column, computed inline — no function call risk
    g = max(K - S(2:Ns), 0);

    % Solve the LCP
    V_new_inner = PSOR_proj(A, b, g, omega, tol, max_iter, V_old_inner);

    % Rebuild full solution
    V = [K; V_new_inner; 0];
    V_old_inner = V_new_inner;
    V_surf(:, n)    = V;   

    % Early-exercise boundary
    [xCrit, ~] = findCriticalPoints(S, V, max(K - S, 0));
    if ~isempty(xCrit)
        Sstar(n) = xCrit(1);
    end
end

t_exec = toc;

%%  Figure 1: V(S, t=0) in the continuation region 
g_final  = max(K - S, 0);
cont_idx = V > g_final + 1e-10;

figure('Name','Option Value at t=0','NumberTitle','off');
hold on;
plot(S, V,       'b-',  'LineWidth', 2,   'DisplayName', 'V(S,0) CN-PSOR');
plot(S, g_final, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Intrinsic g(S)');

S_cont = S(cont_idx);  V_cont = V(cont_idx);
if ~isempty(S_cont)
    fill([S_cont; flipud(S_cont)], ...
         [V_cont; flipud(g_final(cont_idx))], ...
         [0.8 0.9 1.0], 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
         'DisplayName', 'Continuation region');
end
if ~isnan(Sstar(1))
    xline(Sstar(1), 'k:', 'LineWidth', 1.5, ...
          'DisplayName', sprintf('S^*(0) \\approx %.3f', Sstar(1)));
end
xlabel('Stock price S', 'FontSize', 13);
ylabel('Option value V', 'FontSize', 13);
title(sprintf('American Put V(S,0): K=%g, r=%g, \\sigma=%g, T=%g', ...
              K, r, sigma, T), 'FontSize', 13);
legend('Location','northeast','FontSize',11);
grid on; box on; xlim([0 Ss]); ylim([0 K*1.1]);
hold off;

%% ── Figure 2: Continuation region 
t_bnd   = t(1:Nt);
valid   = ~isnan(Sstar);
t_valid = t_bnd(valid);
S_valid = Sstar(valid);
 
figure('Name','Continuation Region','NumberTitle','off');
hold on;
 
% Fill continuation region with horizontal strips coloured by time t
% Each strip spans from S*(t) to Ss at height t — colour = t value
nv = numel(t_valid);
for k = 1:nv-1
    t1 = t_valid(k);   t2 = t_valid(k+1);
    s1 = S_valid(k);   s2 = S_valid(k+1);
    patch([s1 Ss Ss s2], [t1 t1 t2 t2], (t1+t2)/2, 'EdgeColor','none');
end
colormap(gca, parula);
clim([0 T]);
 
windowSize = max(1, round(nv/50));
S_sm = movmean(S_valid, windowSize);
plot(S_sm, t_valid, 'k-', 'LineWidth', 2);
 
xlabel('S', 'FontSize', 13);
ylabel('t', 'FontSize', 13);
title('Continuation region', 'FontWeight', 'bold', 'FontSize', 13);
set(gca, 'Color', 'w');
grid on; box on;
xlim([7 Ss]); ylim([0 T]);
hold off;

end