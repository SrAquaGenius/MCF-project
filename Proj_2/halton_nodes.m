function p = halton_nodes(N, bases)
%HALTON_NODES Halton nodes in [0,1]^d.
%   p = halton_nodes(N) returns N points in [0,1]^2 using bases [2 3].

if nargin < 2
    bases = [2 3];
end

d = numel(bases); %the number of coordinates corresponds to the number of different prime bases to generate Halton nodes.
p = zeros(N, d); %matrix with the number of rows equal to the number of points tht it's pretended to generate, and as number of columns the dimension of the coordinates of the nodes.

for j = 1:d
    for i = 1:N
        p(i, j) = van_der_corput(i, bases(j)); % The Van der Corput sequence is the 1D building block of the Halton sequence.
    end
end
end

function x = van_der_corput(i, base) %This function builds the Van der Corput number by taking the digits of the integer i in a given base, reversing them, and putting them after the decimal point. Actually, the definition of a Van der Crompt sequence.
x = 0;
f = 1/base;

while i > 0
    digit = mod(i, base); % it gets the current last digit.
    x = x + digit*f; %value of the i value in a base with a decimal ofc.
    i = floor(i/base); % it removes the last digit assessed.
    f = f/base; % it moves one place further after the decimal point.
end
end
