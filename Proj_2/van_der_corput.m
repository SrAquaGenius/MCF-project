function x = van_der_corput(i, base)
%VAN_DER_CORPUT Builds Van der Corput number by taking the digits of the
% integer i in a given base, reversing, and putting them after the decimal point.

x = 0;
f = 1/base;

while i > 0
    digit = mod(i, base); % it gets the current last digit.
    x = x + digit * f;    % value of the i value in a base with a decimal ofc.
    i = floor(i / base);  % it removes the last digit assessed.
    f = f / base;         % it moves one place further after the decimal point.
end
end