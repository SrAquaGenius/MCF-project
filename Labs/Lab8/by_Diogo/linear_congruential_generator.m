function [res, t] = linear_congruential_generator(n, seed)
%LINEAR_CONGRUENTIAL_GENERATOR Generate Pseudo-Random numbers from U([0,1])

if (nargout == 2); tic; end

M = 2^31 - 1;
a = 16807;
b = 0;

m = zeros(1, n);
m_aux = seed;

for i = 1:n
    m_aux = mod(a*m_aux+b, M);
    m(i) = m_aux;
end

res = m/M;

if (nargout == 2); t = toc; end
end