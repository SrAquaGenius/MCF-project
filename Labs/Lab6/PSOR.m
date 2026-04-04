function [v] = PSOR(g, n, p, w)
% PSOR returns the PSOR method of a certain funtion

    u = zeros(n,1); % Result vector
    A = 2*diag(ones(n,1)) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);
    b = -A*g(:);

    stop = 1e-6;    % Stop parameter

    %K = zeros(p,n);

    for k = 1:p
        u_old = u;
        
        for i = 1:n
            sum1 = A(i,1:i-1) * u(1:i-1);
            sum2 = A(i,i+1:n) * u_old(i+1:n);
            u_GS = (b(i) - sum1 - sum2) / A(i,i);
            u_relax = u_old(i) + w * (u_GS - u_old(i));
            u(i) = max(0, u_relax);
        end

        %K(k,:) = u'
        
        % Stop Criterium
        if norm(u - u_old, inf) < stop
            break
        end
    end

    v_inner = u + g(:);
    v = [0; v_inner(:); 0];
end