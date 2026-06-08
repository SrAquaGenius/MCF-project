function p = halton_nodes(N, bases)
%HALTON_NODES Halton nodes in [0,1]^d.
%   p = halton_nodes(N) returns N points in [0,1]^2 using bases [2 3].

if nargin < 2
    bases = [2 3];
end

d = numel(bases); % Number of coordinates
% (Each base requires a different prime number to generate Halton nodes.)
p = zeros(N, d); % Matrix with number of points by the dimension.

for j = 1:d
    for i = 1:N
        % The Van der Corput sequence is the 1D building block.
        p(i, j) = van_der_corput(i, bases(j));
    end
end
end
