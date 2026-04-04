function [v, it, residuals, energy, violation] = PSORextra(g, n, p, w)
% PSORextra - PSOR com métricas para análise
%
% OUTPUTS:
% v           -> solução final (com fronteiras)
% it          -> número de iterações realizadas
% residuals   -> histórico ||u^{k+1} - u^k||_inf
% energy      -> energia final da solução
% violation   -> violação do obstáculo

    %% Inicialização
    u = zeros(n,1);
    A = 2*diag(ones(n,1)) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);
    b = -A*g(:);

    stop = 1e-6;
    residuals = zeros(p,1);

    %% Iterações PSOR
    for k = 1:p
        u_old = u;

        for i = 1:n
            sum1 = A(i,1:i-1) * u(1:i-1);
            sum2 = A(i,i+1:n) * u_old(i+1:n);

            u_GS = (b(i) - sum1 - sum2) / A(i,i);

            % Relaxação
            u_relax = u_old(i) + w * (u_GS - u_old(i));

            % Projeção (obstacle problem)
            u(i) = max(0, u_relax);
        end

        residuals(k) = norm(u - u_old, inf);
        if residuals(k) < stop
            break
        end
    end

    % Iterações efetivas
    it = k;
    residuals = residuals(1:it);

    %% Reconstrução da solução completa
    v_inner = u + g(:);
    v = [0; v_inner(:); 0];

    % Métricas adicionais
    energy = sum(diff(u).^2);
end