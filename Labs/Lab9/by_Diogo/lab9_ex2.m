%% Initialize Values
clc; clearvars; close all;

tic
seed = 0.1; % Must be a value in [0,1]

D = 1:3;
%N = [100, 1000, 10000, 100000];
n = 1000;

for d = D

    U = zeros(d, n);
    

    for i = 1:d
        U(i, :) = rand(1, n);
    end

    S = sum(U.^2<=1);





end



fprintf('The total time elapsed was %.4f seconds\n', toc);
