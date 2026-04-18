%I don´t understand if it will be necessary to define critical points or
%inPSOR we can define the boundarios with the conditions of the system.
function x = PSOR_proj(A, b, g, w, tol, max_iter)

n = length(b);
x = g; % initial guess satisfies constraint
b_tilde = b - A*g(:);
for k = 1:max_iter
    x_old = x;
    
    for i = 1:n
        sigma1 = A(i,1:i-1)*x(1:i-1);
        sigma2 = A(i,i+1:n)*x_old(i+1:n);
        
        x_PSOR = (b_tilde(i) - sigma1 - sigma2)/A(i,i);
        
        % Relaxation
        x_relax = x_old(i) + w*(x_PSOR - x_old(i));
        
        % Projection (obstacle condition)
        x(i) = max(g(i), x_relax);
    end
    
    %convergence check
    if norm(x - x_old, inf) < tol
        break;
    end
end

end