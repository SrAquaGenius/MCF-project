function [z, t, s_out] = box_muller_method(n, seed)
%BOX_MULLER_METHOD Generates standard normal variables using Box-Muller

if (nargout >= 2); tic; end

[r, ~, s_in] = ray_distribution(n, seed, 1);

if (nargout == 3)
    [u, ~, s_out] = uni_distribution(n, s_in, 0, 2*pi);
elseif (nargout <= 2)
    u = uni_distribution(n, s_in, 0, 2*pi);
end

x = r.*cos(u);
y = r.*sin(u);
z = [x(:); y(:)];

if (nargout >= 2); t = toc; end
end