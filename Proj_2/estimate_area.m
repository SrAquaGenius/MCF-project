function areas = estimate_area(points, f, alpha)
%ESTIMATE_AREA Estimate |A_alpha| for f(x,y) < alpha.

x = points(:, 1);
y = points(:, 2);
values = f(x,y);

values = sort(values);
N = numel(values); % Number of values defined above for the function f

areas = zeros(size(alpha)); % Matrix of zeros with the same size as alpha. 
for k = 1:numel(alpha) % Loop for one of the 1000 alphas defined.
    areas(k) = upper_count(values, alpha(k))/N; % This counts how many function values are strictly smaller than each index of alpha. Then it' divided by N to get a proportion of the values that hold the condition. 
end
end

function c = upper_count(sorted_values, threshold)
% Number of entries strictly smaller than threshold, it will be equal to c.
lo = 1;
hi = numel(sorted_values);
c = 0;

while lo <= hi
    mid = floor((lo + hi)/2);
    if sorted_values(mid) < threshold
        c = mid;
        lo = mid + 1;
    else
        hi = mid - 1;
    end
end
end
