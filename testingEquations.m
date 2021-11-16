i = 10;
j = 10;
gamma = 1.4;
xMax = 35;
E = 10;
hMin = 40;
theta = 5.352*pi/180;
R = 8.314463;
Cy = 0.1;


base = zeros(j, i);
base(:,1) = ones(j,1);


u = 0.678e3 * base;
v = 0 * base;
ro = 0.123e1 * base;
p = 0.101e6 * base;
T = 0.286e3 * base;
M = 0.200e1 * base;

F1 = ro.*u;
F2 = ro.*u.^2 + p;
F3 = ro.*u.*v;
F4 = gamma / (gamma - 1) * p.*u + ro.*u.*(u.^2 + v.^2) / 2;

G1 = ro.*F3./F1;
G2 = F3;
G3 = ro.*(F3./F1).^2+F2-F1.^2./ro;
G4 = gamma / (gamma - 1) * (F2 - F1.^2./ro)*F3./F1 + ro.*F3./F1/2*((F1./ro).^2+(F3./F1)^2);

F1pre = zeros(j, i);
F2pre = zeros(j, i);
F3pre = zeros(j, i);
F4pre = zeros(j, i);

G1pre = zeros(j, i);
G2pre = zeros(j, i);
G3pre = zeros(j, i);
G4pre = zeros(j, i);

dF1de = zeros(j, i);
dF2de = zeros(j, i);
dF3de = zeros(j, i);
dF4de = zeros(j, i);

dF1dePre = zeros(j, i);
dF2dePre = zeros(j, i);
dF3dePre = zeros(j, i);
dF4dePre = zeros(j, i);

for ho=1:i-1
    for ve=1:j-1
        dx = xMax / i;
        x = dx * ho;
        if x < E
            ys = 0;
            h = hMin;
        else
            h = hMin + (x - E) * tan(theta);
            ys = -(x - E) * tan(theta);
        end
        dy = h / j;
        y = dy * ve;
        eta = (y - ys) / h;
        dEta = 1 / j;

        if x < E
            dEtadX = 0;
        else
            dEtadX = (1 - eta) * tan(theta) / h;
        end
        dF1de(ho, ve) = dEtadX * (F1(ho, ve) - F1(ho, ve + 1)/ dEta) + 1 / h * ((G1(ho, ve) - G1(ho, ve+1))/dEta);
        dF2de(ho, ve) = dEtadX * (F2(ho, ve) - F2(ho, ve + 1)/ dEta) + 1 / h * ((G2(ho, ve) - G2(ho, ve+1))/dEta);
        dF3de(ho, ve) = dEtadX * (F3(ho, ve) - F3(ho, ve + 1)/ dEta) + 1 / h * ((G3(ho, ve) - G3(ho, ve+1))/dEta);
        dF4de(ho, ve) = dEtadX * (F4(ho, ve) - F4(ho, ve + 1)/ dEta) + 1 / h * ((G4(ho, ve) - G4(ho, ve+1))/dEta);
        mu = asin(1/M(ho, ve));
        dE = C * dEta / max([tan(theta - mu) tan(theta + mu)]);

        % Artificial viscosity:
        SF1(ho,ve) = Cy * abs(p(ho, ve + 1) - 2*p(ho,ve) + p(ho, ve - 1)) / (p(ho, ve + 1) + 2*p(ho,ve) + p(ho, ve - 1)) * (F1(ho, ve + 1) - 2*F1(ho,ve) + F1(ho,ve-1));
        SF2(ho,ve) = Cy * abs(p(ho, ve + 1) - 2*p(ho,ve) + p(ho, ve - 1)) / (p(ho, ve + 1) + 2*p(ho,ve) + p(ho, ve - 1)) * (F2(ho, ve + 1) - 2*F2(ho,ve) + F2(ho,ve-1));
        SF3(ho,ve) = Cy * abs(p(ho, ve + 1) - 2*p(ho,ve) + p(ho, ve - 1)) / (p(ho, ve + 1) + 2*p(ho,ve) + p(ho, ve - 1)) * (F3(ho, ve + 1) - 2*F3(ho,ve) + F3(ho,ve-1));
        SF4(ho,ve) = Cy * abs(p(ho, ve + 1) - 2*p(ho,ve) + p(ho, ve - 1)) / (p(ho, ve + 1) + 2*p(ho,ve) + p(ho, ve - 1)) * (F4(ho, ve + 1) - 2*F4(ho,ve) + F4(ho,ve-1));


        % Seguramente sobre estructura matriz.
        F1pre(ho + 1, ve) = F1(ho, ve) + dF1de * dE + SF1(ho, ve);
        F2pre(ho + 1, ve) = F2(ho, ve) + dF2de * dE + SF2(ho, ve);
        F3pre(ho + 1, ve) = F3(ho, ve) + dF3de * dE + SF3(ho, ve);
        F4pre(ho + 1, ve) = F4(ho, ve) + dF4de * dE + SF4(ho, ve);

        % Predicted values:
        A = F3pre(ho + 1, ve) ^2 / 2 / F1pre(ho + 1, ve) - F4pre(ho + 1, ve);
        B = gamma / (gamma - 1) * F1pre(ho + 1, ve) * F2pre(ho + 1, ve);
        C = - (gamma + 1) / 2 / (gamma - 1) * F1pre(ho + 1, ve)^3;
        predictedRo = (-B + sqrt(B^2 - 4 * A * C)) / 2 / A;

        % Predicted G
        G1pre(ho + 1, ve) = predictedRo * F3pre(ho + 1, ve) / F1pre(ho + 1, ve);
        G2pre(ho + 1, ve) = F3pre;
        G3pre(ho + 1, ve) = predictedRo * (F3pre(ho + 1, ve) / F1pre(ho + 1, ve))^2 + F2pre(ho + 1, ve) - F1pre(ho + 1, ve)^2 / predictedRo;
        G4pre(ho + 1, ve) = gamma / (gamma - 1) * (F2pre(ho + 1, ve) - F1pre(ho + 1, ve)^2 / predictedRo) * F3pre(ho + 1, ve) / F1pre(ho + 1, ve) + predictedRo / 2 * F3pre(ho + 1, ve) / F1pre(ho + 1, ve) * ((F1pre(ho + 1, ve)/predictedRo)^2 + (F3pre(ho + 1, ve) / F1pre(ho + 1, ve))^2);

        if(ve ~= 1) % No boundary
            % Predicted dFde
            dF1dePre(ho + 1, ve) = dEtadX * (F1(ho+1, ve-1) - F1(ho + 1, ve)/ dEta) + 1 / h * ((G1(ho+1, ve-1) - G1(ho + 1, ve))/dEta);
            dF2dePre(ho + 1, ve) = dEtadX * (F2(ho+1, ve-1) - F2(ho + 1, ve)/ dEta) + 1 / h * ((G2(ho+1, ve-1) - G2(ho + 1, ve))/dEta);
            dF3dePre(ho + 1, ve) = dEtadX * (F3(ho+1, ve-1) - F3(ho + 1, ve)/ dEta) + 1 / h * ((G3(ho+1, ve-1) - G3(ho + 1, ve))/dEta);
            dF4dePre(ho + 1, ve) = dEtadX * (F4(ho+1, ve-1) - F4(ho + 1, ve)/ dEta) + 1 / h * ((G4(ho+1, ve-1) - G4(ho + 1, ve))/dEta);

            % Averaging
            dF1deAvg = (dF1de(ho, ve) + dF1dePre(ho + 1, ve)) / 2;
            dF2deAvg = (dF2de(ho, ve) + dF2dePre(ho + 1, ve)) / 2;
            dF3deAvg = (dF3de(ho, ve) + dF3dePre(ho + 1, ve)) / 2;
            dF4deAvg = (dF4de(ho, ve) + dF4dePre(ho + 1, ve)) / 2;

            % Predicted pressure
            pPre(ho, ve) = F2pre - F1pre^2 / predictedRo;

            % Predicted SF
            SF1pre(ho + 1,ve) = Cy * abs(pPre(ho, ve + 1) - 2*pPre(ho,ve) + pPre(ho, ve - 1)) / (pPre(ho, ve + 1) + 2*pPre(ho,ve) + pPre(ho, ve - 1)) * (F1pre(ho, ve + 1) - 2*F1pre(ho,ve) + F1pre(ho,ve-1));
            SF2pre(ho + 1,ve) = Cy * abs(pPre(ho, ve + 1) - 2*pPre(ho,ve) + pPre(ho, ve - 1)) / (pPre(ho, ve + 1) + 2*pPre(ho,ve) + pPre(ho, ve - 1)) * (F2pre(ho, ve + 1) - 2*F2pre(ho,ve) + F2pre(ho,ve-1));
            SF3pre(ho + 1,ve) = Cy * abs(pPre(ho, ve + 1) - 2*pPre(ho,ve) + pPre(ho, ve - 1)) / (pPre(ho, ve + 1) + 2*pPre(ho,ve) + pPre(ho, ve - 1)) * (F3pre(ho, ve + 1) - 2*F3pre(ho,ve) + F3pre(ho,ve-1));
            SF4pre(ho + 1,ve) = Cy * abs(pPre(ho, ve + 1) - 2*pPre(ho,ve) + pPre(ho, ve - 1)) / (pPre(ho, ve + 1) + 2*pPre(ho,ve) + pPre(ho, ve - 1)) * (F4pre(ho, ve + 1) - 2*F4pre(ho,ve) + F4pre(ho,ve-1));

            % New F
            F1(ho + 1, ve) = F1(ho, ve) + dF1deAvg * dE + SF1pre(ho+1,ve);
            F2(ho + 1, ve) = F2(ho, ve) + dF2deAvg * dE + SF2pre(ho+1,ve);
            F3(ho + 1, ve) = F3(ho, ve) + dF3deAvg * dE + SF3pre(ho+1,ve);
            F4(ho + 1, ve) = F4(ho, ve) + dF4deAvg * dE + SF4pre(ho+1,ve);

            % Magnitudes fisicas
            A = F3(ho, ve) ^2 / 2 / F1(ho, ve) - F4(ho, ve);
            B = gamma / (gamma - 1) * F1(ho, ve) * F2(ho, ve);
            C = - (gamma + 1) / 2 / (gamma - 1) * F1(ho, ve)^3;
            ro(ho,ve) = (-B + sqrt(B^2 - 4 * A * C)) / 2 / A;

            u(ho,ve) = F1(ho,ve) / ro(ho,ve);
            v(ho,ve) = F3(ho,ve) / F1(ho,ve);
            p(ho,ve) = F2(ho,ve) - F1(ho,ve)*u(ho,ve);
            T(ho,ve) = p(ho,ve) / ro(ho,ve) / R;
        else % Boundary
            % Predicted dFde
            dF1dePre(ho + 1, ve) = dEtadX * (F1(ho+1, ve+1) - F1(ho + 1, ve)/ dEta) + 1 / h * ((G1(ho+1, ve+1) - G1(ho + 1, ve))/dEta);
            dF2dePre(ho + 1, ve) = dEtadX * (F2(ho+1, ve+1) - F2(ho + 1, ve)/ dEta) + 1 / h * ((G2(ho+1, ve+1) - G2(ho + 1, ve))/dEta);
            dF3dePre(ho + 1, ve) = dEtadX * (F3(ho+1, ve+1) - F3(ho + 1, ve)/ dEta) + 1 / h * ((G3(ho+1, ve+1) - G3(ho + 1, ve))/dEta);
            dF4dePre(ho + 1, ve) = dEtadX * (F4(ho+1, ve+1) - F4(ho + 1, ve)/ dEta) + 1 / h * ((G4(ho+1, ve+1) - G4(ho + 1, ve))/dEta);

            % Averaging
            dF1deAvg = (dF1de(ho, ve) + dF1dePre(ho + 1, ve)) / 2;
            dF2deAvg = (dF2de(ho, ve) + dF2dePre(ho + 1, ve)) / 2;
            dF3deAvg = (dF3de(ho, ve) + dF3dePre(ho + 1, ve)) / 2;
            dF4deAvg = (dF4de(ho, ve) + dF4dePre(ho + 1, ve)) / 2;

            % Predicted pressure
            pPre(ho, ve) = F2pre - F1pre^2 / predictedRo;

            % Predicted SF
            SF1pre(ho + 1,ve) = Cy * abs(pPre(ho, ve + 1) - 2*pPre(ho,ve) + pPre(ho, ve - 1)) / (pPre(ho, ve + 1) + 2*pPre(ho,ve) + pPre(ho, ve - 1)) * (F1pre(ho, ve + 1) - 2*F1pre(ho,ve) + F1pre(ho,ve-1));
            SF2pre(ho + 1,ve) = Cy * abs(pPre(ho, ve + 1) - 2*pPre(ho,ve) + pPre(ho, ve - 1)) / (pPre(ho, ve + 1) + 2*pPre(ho,ve) + pPre(ho, ve - 1)) * (F2pre(ho, ve + 1) - 2*F2pre(ho,ve) + F2pre(ho,ve-1));
            SF3pre(ho + 1,ve) = Cy * abs(pPre(ho, ve + 1) - 2*pPre(ho,ve) + pPre(ho, ve - 1)) / (pPre(ho, ve + 1) + 2*pPre(ho,ve) + pPre(ho, ve - 1)) * (F3pre(ho, ve + 1) - 2*F3pre(ho,ve) + F3pre(ho,ve-1));
            SF4pre(ho + 1,ve) = Cy * abs(pPre(ho, ve + 1) - 2*pPre(ho,ve) + pPre(ho, ve - 1)) / (pPre(ho, ve + 1) + 2*pPre(ho,ve) + pPre(ho, ve - 1)) * (F4pre(ho, ve + 1) - 2*F4pre(ho,ve) + F4pre(ho,ve-1));

            % New F
            F1(ho + 1, ve) = F1(ho, ve) + dF1deAvg * dE + SF1pre(ho+1,ve);
            F2(ho + 1, ve) = F2(ho, ve) + dF2deAvg * dE + SF2pre(ho+1,ve);
            F3(ho + 1, ve) = F3(ho, ve) + dF3deAvg * dE + SF3pre(ho+1,ve);
            F4(ho + 1, ve) = F4(ho, ve) + dF4deAvg * dE + SF4pre(ho+1,ve);

            % Magnitudes fisicas
            A = F3(ho, ve) ^2 / 2 / F1(ho, ve) - F4(ho, ve);
            B = gamma / (gamma - 1) * F1(ho, ve) * F2(ho, ve);
            C = - (gamma + 1) / 2 / (gamma - 1) * F1(ho, ve)^3;
            ro(ho,ve) = (-B + sqrt(B^2 - 4 * A * C)) / 2 / A;

            u(ho,ve) = F1(ho,ve) / ro(ho,ve);
            v(ho,ve) = F3(ho,ve) / F1(ho,ve);
            p(ho,ve) = F2(ho,ve) - F1(ho,ve)*u(ho,ve);
            T(ho,ve) = p(ho,ve) / ro(ho,ve) / R;
            M(ho,ve) = sqrt(v(ho,ve)^2 + u(ho,ve)^2) / sqrt(gamma * p(ho,ve) / ro(ho,ve));
            if x > E
                phi = theta - atan2(abs(v(ho,ve)) , u(ho,ve));
            else
                phi = atan2(v(ho,ve) , u(ho,ve));
            end

            f_cal = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (M(ho,ve)^2 - 1))) - atan(sqrt(M(ho,ve)^2 - 1));

            f_act = f_cal + phi;

            var = 10;
            tol = 1e-3;
            MAct(ho,ve) = 0;
            while (f_act - guessed > tol)
                guessed = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (MAct^2 - 1))) - atan(sqrt(MAct^2 - 1));
                MAct(ho,ve) = MAct(ho,ve) + 1e-3;
            end

            pAct(ho,ve) = p(ho,ve) * ((1 + ((gamma - 1) / 2) * M(ho,ve)^ 2) / (1 + ((gamma - 1) / 2) * MAct(ho,ve)^ 2))^(gamma / (gamma - 1));
            TAct(ho,ve) = T(ho,ve) * ((1 + ((gamma - 1) / 2) * M(ho,ve)^ 2) / (1 + ((gamma - 1) / 2) * MAct(ho,ve)^ 2));
            roAct(ho, ve) = pAct(ho,ve) / R / TAct(ho,ve);

            p(ho,ve) = pAct(ho,ve);
            T(ho,ve) = TAct(ho,ve);
            ro(ho,ve) = roAct(ho,ve);





        end
    end
end