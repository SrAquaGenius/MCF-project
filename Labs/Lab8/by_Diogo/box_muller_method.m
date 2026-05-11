function [x, y, t] = box_muller_method(n, seed)
%BOX_MULLER_METHOD Generates standard normal variables using Box-Muller

tic

[r, ~] = ray_distribution(n, seed, 1);
[u, ~] = uni_distribution(n, seed + 0.37, 0, 2*pi);

x = r.*cos(u);
y = r.*sin(u);

t = toc;
end