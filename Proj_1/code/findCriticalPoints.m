function [xCrit, vCrit] = findCriticalPoints(S, v, g)
% FINDCRITICALPOINTS: Find the points where v' ≈ g'
%
% Input:
%   x - Points vector of the grid
%   v - Solution (u + g)
%   g - Obstacle function assessed in the points x
%
% Output:
%   idx  - Index of the critical points of the grid
%   xCrit - Position of the critical points
S = S(:);
V = v(:);
g = g(:);
 
h   = V - g;            % positive in continuation, zero in stopping
tol = 1e-10;
 
% Detect signal changes (zero crossings of h)
cross = find(h(1:end-1) .* h(2:end) <= 0 & ...
             (abs(h(1:end-1)) > tol | abs(h(2:end)) > tol));
 
nc    = numel(cross);
xCrit = zeros(nc, 1);
vCrit = zeros(nc, 1);
 
for k = 1:nc
    i  = cross(k);
    s1 = S(i);   s2 = S(i+1);
    h1 = h(i);   h2 = h(i+1);
 
    % Linear interpolation for the zero of h
    s0 = s1 - h1 * (s2 - s1) / (h2 - h1);
 
    xCrit(k) = s0;
    vCrit(k) = V(i) + (V(i+1) - V(i)) * (s0 - s1) / (s2 - s1);
end
end