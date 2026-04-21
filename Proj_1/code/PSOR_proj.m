function [x, iters] = PSOR_proj(A, b, g, omega, tol, max_iter, x_init)
% PSOR_PROJ  PSOR for the LCP defined by the 3 constraints

% Create column vectors
b = b(:);
g = g(:);
n = numel(b);

% Initial guess
if nargin < 7 || isempty(x_init)
    x = max(g, b ./ diag(A));   
else
    x = max(g, x_init(:));
end

% Extract diagonals once 
%   dl(i) = A(i+1, i)  sub-diagonal   (length n-1)
%   d(i)  = A(i,   i)  main diagonal  (length n)
%   du(i) = A(i, i+1)  super-diagonal (length n-1)
d  = diag(A,  0);
dl = diag(A, -1);
du = diag(A, +1);

iters = max_iter;

for k = 1:max_iter

    x_old = x;

    for i = 1:n
        % Gauss-Seidel sum: lower uses updated x, upper uses x_old
        lo = 0;  if i > 1, lo = dl(i-1) * x(i-1);     end
        hi = 0;  if i < n, hi = du(i)   * x_old(i+1); end

        x_gs = (b(i) - lo - hi) / d(i);

        % The projection onto x >= g
        x(i) = max(g(i), x_old(i) + omega*(x_gs - x_old(i)));
    end

    % Complementarity residual: it should be zero at solution
    r_vec    = A*x - b;
    comp_res = min(x - g, r_vec);
    if norm(comp_res, inf) < tol
        iters = k;
        break;
    end
end
end