function [dEtadX, dXi, dEta, h] = initVars(j, xi, E, H, theta, ve, ho, M, Courant)
        if xi <= E
            ys = 0;
            h = H;
        else
            h = H + (xi - E) * tan(theta);
            ys = -(xi - E) * tan(theta);
        end
        dy = h / j;
        y = dy * (ve-1);
        eta = (y - ys) / h;
        dEta = 1 / j;

        if xi < E
            dEtadX = 0;
        else
            dEtadX = (1 - eta) * tan(theta) / h;
        end
        mu = asin(1/M(ve, ho));
        dXi = Courant * dy / max([tan(theta - mu) tan(theta + mu)]); % Hard to check
end
