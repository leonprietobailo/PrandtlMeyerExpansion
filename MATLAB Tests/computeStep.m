function [dXi] = computeStep(M, theta, ho, Courant, dy)

anglePlus = [];
angleMin = [];
for n=1:size(M,1)
    mu = asin(1 / M(n, ho));
    anglePlus(end + 1) = abs(tan(theta + mu));
    angleMin(end + 1) = abs(tan(theta - mu));
end
angleMAX = max([max(anglePlus) max(angleMin)]);

dXi = Courant * dy / angleMAX; % Hard to check

end

