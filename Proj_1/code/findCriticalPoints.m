function [xCrit, yCrit] = findCriticalPoints(x, v, g)
% FINDCRITICALPOINTS encontra pontos onde v' ≈ g'
%
% Entrada:
%   x - vetor de pontos da malha
%   v - solução (u + g)
%   g - função obstáculo avaliada nos pontos x
%
% Saída:
%   idx  - índices dos pontos críticos na malha
%   xCrit - posições dos pontos críticos

    dx = x(2) - x(1);
    dv = gradient(v, dx);
    dg = gradient(g, dx);
    h = dv(:) - dg(:);
    
    % Detect indexes where the zero crossing occurs
    tol = 1e-6;
    idx = find(h(1:end-1).*h(2:end) <= 0 & ...
          (abs(h(1:end-1)) > tol | abs(h(2:end)) > tol));
    n_idx = length(idx);
    
    % Aproximate critical point
    xCrit = zeros(n_idx,1);
    yCrit = zeros(n_idx,1);

    for k = 1:n_idx
        i = idx(k);
        x1 = x(i);   x2 = x(i+1);
        h1 = h(i);   h2 = h(i+1);
        x0 = x1 - h1*(x2 - x1)/(h2 - h1);
        xCrit(k) = x0;
        v1 = v(i); v2 = v(i+1);
        yCrit(k) = v1 + (v2 - v1)*(x0 - x1)/(x2 - x1);
    end
end