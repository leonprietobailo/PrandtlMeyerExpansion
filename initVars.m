function [dEtadX, dEta, h, dy] = initVars(j, xi, E, H, theta)
for i = 1:j
    if xi < E
        ys = 0;
        h = H;
    else
        h = H + (xi - E) * tan(theta);
        ys = -(xi - E) * tan(theta);
    end
    dy = h / (j - 1);
    y(i) = dy * (j-1);
    eta(i) = (y(i) - ys) / h;
    dEta = 1 / (j - 1);
    
    if xi < E
        dEtadX(i) = 0;
    else
        dEtadX(i) = (1 - eta(i)) * tan(theta) / h;
    end
end
end
