function [res,t] = inversion_method(n, f, arg, seed)
%INVERSION_METHOD Generate the sample using for the given inverted function

tic
[u, ] = linear_congruential_generator(n, seed);
res = f(arg, u);
t = toc;
end