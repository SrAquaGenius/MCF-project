function areas = estimate_area_curve(points, alpha)
%ESTIMATE_AREA_CURVE Estimate |A_alpha| for f(x,y) < alpha.

x = points(:, 1);
y = points(:, 2);
values = sin(10*(x.^2 - sin(3*y)));

values = sort(values);
N = numel(values); %number of values defined above for the function f

areas = zeros(size(alpha)); %matrix of zeros with the same size as alpha. 
for k = 1:numel(alpha) %loop for one of the 1000 alphas defined.
    areas(k) = upper_count(values, alpha(k))/N; %This counts how many function values are strictly smaller than each index of alpha. Then it' divided by N to get a proportion of the values that hold the condition. 
end
end

function c = upper_count(sorted_values, threshold)
% Number of entries strictly smaller than threshold, it will correspond to
% c.
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
