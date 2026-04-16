function plot_surface(s, t, U, type, op)

[S,T] = meshgrid(s,t);

figure
surf(S,T,U,'FaceLighting','phong','EdgeColor','none','FaceColor','interp');

if strcmpi(op,'put')
    view(45, 25)
else
    view(-45, 25)
end

xlabel('Stock price S')
ylabel('Time t')
zlabel('Option value V')

if strcmpi(op,'put'), Op = 'Put'; else Op = 'Call'; end
if strcmpi(type,'am'), Type = 'American'; else Type = 'European'; end

title(['Black-Scholes ' Op ' Price Surface on an ' Type ' Option'])
colorbar
hold on

% Grelha visual
t_vals = linspace(t(1), t(end), 10);
s_vals = linspace(s(1), s(end), 6);

for tt = t_vals
    [~, idx_t] = min(abs(t - tt));
    plot3(s, t(idx_t)*ones(size(s)), U(idx_t,:), 'k-', 'LineWidth', 0.8);
end

for ss = s_vals
    [~, idx_s] = min(abs(s - ss));
    plot3(s(idx_s)*ones(size(t)), t, U(:,idx_s), 'k-', 'LineWidth', 0.8);
end
end
