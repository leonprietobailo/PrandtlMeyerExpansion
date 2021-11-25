clear;

j = 41;
gamma = 1.4;
xMax = 35;
E = 10;
H = 40;
theta = 5.352*pi/180;
R = 287;
Cy = 0.1;
Courant = 0.5;
close all;

base(:,1) = ones(j,1);

uVal = 0.678e3;
vVal = 0;
roVal = 0.123e1;
pVal = 0.101e6;
TVal = 0.286e3;
MVal = 0.200e1;

u = uVal * base;
v = vVal * base;
ro = roVal * base;
p = pVal * base;
T = TVal * base;
M = MVal * base;

F1Val = roVal * uVal;
F2Val = roVal * uVal ^ 2 + pVal;
F3Val = roVal * uVal * vVal;
F4Val = gamma / (gamma - 1) * pVal*uVal + roVal*uVal*(uVal^2 + vVal^2) / 2;

F1 = F1Val * base;
F2 = F2Val * base;
F3 = F3Val * base;
F4 = F4Val * base;

G1Val = roVal * F3Val / F1Val;
G2Val = F3Val;
G3Val = roVal*(F3Val / F1Val)^2+F2Val - F1Val ^2 / roVal;
G4Val = gamma / (gamma - 1) * (F2Val - F1Val^2/roVal)*F3Val/F1Val + roVal*F3Val/F1Val/2*((F1Val/roVal)^2+(F3Val/F1Val)^2);

G1 = G1Val * base;
G2 = G2Val * base;
G3 = G3Val * base;
G4 = G4Val * base;

% F1pre = zeros(j, i);
% F2pre = zeros(j, i);
% F3pre = zeros(j, i);
% F4pre = zeros(j, i);
%
% G1pre = zeros(j, i);
% G2pre = zeros(j, i);
% G3pre = zeros(j, i);
% G4pre = zeros(j, i);
%
% dF1de = zeros(j, i);
% dF2de = zeros(j, i);
% dF3de = zeros(j, i);
% dF4de = zeros(j, i);
%
% dF1dePre = zeros(j, i);
% dF2dePre = zeros(j, i);
% dF3dePre = zeros(j, i);
% dF4dePre = zeros(j, i);

test = 0;

pPre = p;
n = 0;
xi = 0;
%for ho=1:i-1
ho = 1;

while(xi < xMax)
    %% PREDICTOR STEP
    for ve=1:j
        % Init Vars
        [dEtadX, dXi, dEta, h] = initVars(j, xi, E, H, theta, ve, ho, M, Courant);
        
        if(ve ~= 1 && ve ~= j)
            % Fowrard differences.
            dF1de(ve, ho) = dEtadX * (F1(ve, ho) - F1(ve+1,ho)/ dEta) + 1 / h * ((G1(ve, ho) - G1(ve+1,ho))/dEta);
            dF2de(ve, ho) = dEtadX * (F2(ve, ho) - F2(ve+1,ho)/ dEta) + 1 / h * ((G2(ve, ho) - G2(ve+1,ho))/dEta);
            dF3de(ve, ho) = dEtadX * (F3(ve, ho) - F3(ve+1,ho)/ dEta) + 1 / h * ((G3(ve, ho) - G3(ve+1,ho))/dEta);
            dF4de(ve, ho) = dEtadX * (F4(ve, ho) - F4(ve+1,ho)/ dEta) + 1 / h * ((G4(ve, ho) - G4(ve+1,ho))/dEta);
        elseif ve == 1
            dF1de(ve, ho) = dEtadX * (F1(ve, ho) - F1(ve+1,ho)/ dEta) + 1 / h * ((G1(ve, ho) - G1(ve+1,ho))/dEta);
            dF2de(ve, ho) = dEtadX * (F2(ve, ho) - F2(ve+1,ho)/ dEta) + 1 / h * ((G2(ve, ho) - G2(ve+1,ho))/dEta);
            dF3de(ve, ho) = dEtadX * (F3(ve, ho) - F3(ve+1,ho)/ dEta) + 1 / h * ((G3(ve, ho) - G3(ve+1,ho))/dEta);
            dF4de(ve, ho) = dEtadX * (F4(ve, ho) - F4(ve+1,ho)/ dEta) + 1 / h * ((G4(ve, ho) - G4(ve+1,ho))/dEta);
        
        else
            % Rearwards differences.
            dF1de(ve, ho) = dEtadX * (F1(ve - 1, ho) - F1(ve,ho)/ dEta) + 1 / h * ((G1(ve - 1,ho) - G1(ve,ho))/dEta);
            dF2de(ve, ho) = dEtadX * (F2(ve - 1, ho) - F2(ve,ho)/ dEta) + 1 / h * ((G2(ve - 1,ho) - G2(ve,ho))/dEta);
            dF3de(ve, ho) = dEtadX * (F3(ve - 1, ho) - F3(ve,ho)/ dEta) + 1 / h * ((G3(ve - 1,ho) - G3(ve,ho))/dEta);
            dF4de(ve, ho) = dEtadX * (F4(ve - 1, ho) - F4(ve,ho)/ dEta) + 1 / h * ((G4(ve - 1,ho) - G4(ve,ho))/dEta);
        end
        
        % Predicted F
        F1pre(ve,ho+1) = F1(ve, ho) + dF1de(ve,ho) * dXi;
        F2pre(ve,ho+1) = F2(ve, ho) + dF2de(ve,ho) * dXi;
        F3pre(ve,ho+1) = F3(ve, ho) + dF3de(ve,ho) * dXi;
        F4pre(ve,ho+1) = F4(ve, ho) + dF4de(ve,ho) * dXi;
        
        % Predicted values:
        A = F3pre(ve,ho+1) ^2 / 2 / F1pre(ve,ho+1) - F4pre(ve,ho+1);
        B = gamma / (gamma - 1) * F1pre(ve,ho+1) * F2pre(ve,ho+1);
        C = - (gamma + 1) / 2 / (gamma - 1) * F1pre(ve,ho+1)^3;
        roPre(ve, ho + 1) = (-B + sqrt(B^2 - 4 * A * C)) / 2 / A;
        
        % Predicted G
        G1pre(ve,ho+1) = roPre(ve, ho + 1) * F3pre(ve,ho+1) / F1pre(ve,ho+1);
        G2pre(ve,ho+1) = F3pre(ve,ho+1);
        G3pre(ve,ho+1) = roPre(ve, ho + 1) * (F3pre(ve,ho+1) / F1pre(ve,ho+1))^2 + F2pre(ve,ho+1) - F1pre(ve,ho+1)^2 / roPre(ve, ho + 1);
        G4pre(ve,ho+1) = gamma / (gamma - 1) * (F2pre(ve,ho+1) - F1pre(ve,ho+1)^2 / roPre(ve, ho + 1)) * F3pre(ve,ho+1) / F1pre(ve,ho+1) + roPre(ve, ho + 1) / 2 * F3pre(ve,ho+1) / F1pre(ve,ho+1) * ((F1pre(ve,ho+1)/roPre(ve, ho + 1))^2 + (F3pre(ve,ho+1) / F1pre(ve,ho+1))^2);
    end
    
    for ve=1:j
        % Init Vars
        [dEtadX, dXi, dEta, h] = initVars(j, xi, E, H, theta, ve, ho, M, Courant);
        
        if(ve ~= 1)
            % Predicted dFde [TOP INFLOW]
            dF1dePre(ve,ho+1) = dEtadX * (F1pre(ve-1, ho+1) - F1pre(ve,ho+1)/ dEta) + 1 / h * ((G1pre(ve-1, ho+1) - G1pre(ve,ho+1))/dEta);
            dF2dePre(ve,ho+1) = dEtadX * (F2pre(ve-1, ho+1) - F2pre(ve,ho+1)/ dEta) + 1 / h * ((G2pre(ve-1, ho+1) - G2pre(ve,ho+1))/dEta);
            dF3dePre(ve,ho+1) = dEtadX * (F3pre(ve-1, ho+1) - F3pre(ve,ho+1)/ dEta) + 1 / h * ((G3pre(ve-1, ho+1) - G3pre(ve,ho+1))/dEta);
            dF4dePre(ve,ho+1) = dEtadX * (F4pre(ve-1, ho+1) - F4pre(ve,ho+1)/ dEta) + 1 / h * ((G4pre(ve-1, ho+1) - G4pre(ve,ho+1))/dEta);
            
        else
            % Predicted dFde [BOT]
            dF1dePre(ve,ho+1) = dEtadX * (F1pre(ve, ho+1) - F1pre(ve + 1, ho + 1)/ dEta) + 1 / h * ((G1pre(ve, ho+1) - G1pre(ve + 1, ho + 1))/dEta);
            dF2dePre(ve,ho+1) = dEtadX * (F2pre(ve, ho+1) - F2pre(ve + 1, ho + 1)/ dEta) + 1 / h * ((G2pre(ve, ho+1) - G2pre(ve + 1, ho + 1))/dEta);
            dF3dePre(ve,ho+1) = dEtadX * (F3pre(ve, ho+1) - F3pre(ve + 1, ho + 1)/ dEta) + 1 / h * ((G3pre(ve, ho+1) - G3pre(ve + 1, ho + 1))/dEta);
            dF4dePre(ve,ho+1) = dEtadX * (F4pre(ve, ho+1) - F4pre(ve + 1, ho + 1)/ dEta) + 1 / h * ((G4pre(ve, ho+1) - G4pre(ve + 1, ho + 1))/dEta);
        end
        % Averaging
        dF1deAvg = (dF1de(ve, ho) + dF1dePre(ve,ho+1)) / 2;
        dF2deAvg = (dF2de(ve, ho) + dF2dePre(ve,ho+1)) / 2;
        dF3deAvg = (dF3de(ve, ho) + dF3dePre(ve,ho+1)) / 2;
        dF4deAvg = (dF4de(ve, ho) + dF4dePre(ve,ho+1)) / 2;
        
        % New F
        F1(ve,ho+1) = F1(ve, ho) + dF1deAvg * dXi;
        F2(ve,ho+1) = F2(ve, ho) + dF2deAvg * dXi;
        F3(ve,ho+1) = F3(ve, ho) + dF3deAvg * dXi;
        F4(ve,ho+1) = F4(ve, ho) + dF4deAvg * dXi;
        
        % Magnitudes fisicas
        A = F3(ve, ho+1) ^2 / 2 / F1(ve, ho+1) - F4(ve, ho+1);
        B = gamma / (gamma - 1) * F1(ve, ho+1) * F2(ve, ho+1);
        C = - (gamma + 1) / 2 / (gamma - 1) * F1(ve, ho+1)^3;
        ro(ve,ho+1) = (-B + sqrt(B^2 - 4 * A * C)) / 2 / A;
        
        u(ve,ho+1) = F1(ve,ho+1) / ro(ve,ho+1);
        v(ve,ho+1) = F3(ve,ho+1) / F1(ve,ho+1);
        p(ve,ho+1) = F2(ve,ho+1) - F1(ve,ho+1)*u(ve,ho+1);
        T(ve,ho+1) = p(ve,ho+1) / ro(ve,ho+1) / R;
        M(ve,ho+1) = sqrt(v(ve,ho+1)^2 + u(ve,ho+1)^2) / sqrt(gamma * p(ve,ho+1) / ro(ve,ho+1));
    end
    
    for ve=1:j
        % Init Vars
        [dEtadX, dXi, dEta, h] = initVars(j, xi, E, H, theta, ve, ho, M, Courant);
        
        if(ve ~= 1 && ve ~= j)
            % Viscosity
            SF1(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / (p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho)) * (F1(ve - 1,ho) - 2 * F1(ve,ho) + F1(ve + 1, ho));
            SF2(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / (p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho)) * (F2(ve - 1,ho) - 2 * F2(ve,ho) + F2(ve + 1, ho));
            SF3(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / (p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho)) * (F3(ve - 1,ho) - 2 * F3(ve,ho) + F3(ve + 1, ho));
            SF4(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / (p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho)) * (F4(ve - 1,ho) - 2 * F4(ve,ho) + F4(ve + 1, ho));
%         elseif (ve == 1)
%             % Viscosity
%             SF1(ve,ho) = Cy * abs(p(ve + 1, ho) - p(ve, ho)) / (p(ve + 1, ho) + p(ve, ho)) * (F1(ve + 1, ho) - F1(ve,ho));
%             SF2(ve,ho) = Cy * abs(p(ve + 1, ho) - p(ve, ho)) / (p(ve + 1, ho) + p(ve, ho)) * (F2(ve + 1, ho) - F2(ve,ho));
%             SF3(ve,ho) = Cy * abs(p(ve + 1, ho) - p(ve, ho)) / (p(ve + 1, ho) + p(ve, ho)) * (F3(ve + 1, ho) - F3(ve,ho));
%             SF4(ve,ho) = Cy * abs(p(ve + 1, ho) - p(ve, ho)) / (p(ve + 1, ho) + p(ve, ho)) * (F4(ve + 1, ho) - F4(ve,ho));
%             
%         elseif (ve == j)
%             % Viscosity
%             SF1(ve,ho) = Cy * abs(p(ve - 1, ho) - p(ve, ho)) / (p(ve - 1, ho) + p(ve, ho)) * (F1(ve, ho) - F1(ve - 1,ho));
%             SF2(ve,ho) = Cy * abs(p(ve - 1, ho) - p(ve, ho)) / (p(ve - 1, ho) + p(ve, ho)) * (F2(ve, ho) - F2(ve - 1,ho));
%             SF3(ve,ho) = Cy * abs(p(ve - 1, ho) - p(ve, ho)) / (p(ve - 1, ho) + p(ve, ho)) * (F3(ve, ho) - F3(ve - 1,ho));
%             SF4(ve,ho) = Cy * abs(p(ve - 1, ho) - p(ve, ho)) / (p(ve - 1, ho) + p(ve, ho)) * (F4(ve, ho) - F4(ve - 1,ho));
%         end
        % Prdicted F
        F1pre(ve,ho+1) = F1(ve, ho) + dF1de(ve, ho) * dXi + SF1(ve, ho);
        F2pre(ve,ho+1) = F2(ve, ho) + dF2de(ve, ho) * dXi + SF2(ve, ho);
        F3pre(ve,ho+1) = F3(ve, ho) + dF3de(ve, ho) * dXi + SF3(ve, ho);
        F4pre(ve,ho+1) = F4(ve, ho) + dF4de(ve, ho) * dXi + SF4(ve, ho);
        
        A = F3pre(ve,ho+1) ^2 / 2 / F1pre(ve,ho+1) - F4pre(ve,ho+1);
        B = gamma / (gamma - 1) * F1pre(ve,ho+1) * F2pre(ve,ho+1);
        C = - (gamma + 1) / 2 / (gamma - 1) * F1pre(ve,ho+1)^3;
        roPre(ve, ho + 1) = (-B + sqrt(B^2 - 4 * A * C)) / 2 / A;
        
        pPre(ve,ho+1) = F2pre(ve,ho+1) - F1pre(ve,ho+1)^2 / roPre(ve, ho + 1);
        end
    end
    
    for ve =1:j
        % Init Vars
        [dEtadX, dXi, dEta, h] = initVars(j, xi, E, H, theta, ve, ho, M, Courant);
        
        if(ve ~= 1 && ve ~= j)
            % Viscosity
            SF1pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1,ho)) / (pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1,ho)) * (F1pre(ve - 1,ho) - 2 * F1pre(ve,ho) + F1pre(ve + 1, ho));
            SF2pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1,ho)) / (pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1,ho)) * (F2pre(ve - 1,ho) - 2 * F2pre(ve,ho) + F2pre(ve + 1, ho));
            SF3pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1,ho)) / (pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1,ho)) * (F3pre(ve - 1,ho) - 2 * F3pre(ve,ho) + F3pre(ve + 1, ho));
            SF4pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1,ho)) / (pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1,ho)) * (F4pre(ve - 1,ho) - 2 * F4pre(ve,ho) + F4pre(ve + 1, ho));
%         elseif (ve == 1)
%             % Viscosity
%             SF1pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho) - pPre(ve, ho)) / (pPre(ve + 1, ho) + pPre(ve, ho)) * (F1pre(ve + 1, ho) - F1pre(ve,ho));
%             SF2pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho) - pPre(ve, ho)) / (pPre(ve + 1, ho) + pPre(ve, ho)) * (F2pre(ve + 1, ho) - F2pre(ve,ho));
%             SF3pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho) - pPre(ve, ho)) / (pPre(ve + 1, ho) + pPre(ve, ho)) * (F3pre(ve + 1, ho) - F3pre(ve,ho));
%             SF4pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho) - pPre(ve, ho)) / (pPre(ve + 1, ho) + pPre(ve, ho)) * (F4pre(ve + 1, ho) - F4pre(ve,ho));
%             
%         elseif (ve == j)
%             % Viscosity
%             SF1pre(ve,ho+1) = Cy * abs(pPre(ve - 1, ho) - pPre(ve, ho)) / (pPre(ve - 1, ho) + pPre(ve, ho)) * (F1pre(ve, ho) - F1pre(ve - 1,ho));
%             SF2pre(ve,ho+1) = Cy * abs(pPre(ve - 1, ho) - pPre(ve, ho)) / (pPre(ve - 1, ho) + pPre(ve, ho)) * (F2pre(ve, ho) - F2pre(ve - 1,ho));
%             SF3pre(ve,ho+1) = Cy * abs(pPre(ve - 1, ho) - pPre(ve, ho)) / (pPre(ve - 1, ho) + pPre(ve, ho)) * (F3pre(ve, ho) - F3pre(ve - 1,ho));
%             SF4pre(ve,ho+1) = Cy * abs(pPre(ve - 1, ho) - pPre(ve, ho)) / (pPre(ve - 1, ho) + pPre(ve, ho)) * (F4pre(ve, ho) - F4pre(ve - 1,ho));
%         end
        
        % New F
        F1(ve,ho+1) = F1(ve, ho) + dF1deAvg * dXi + SF1pre(ve,ho+1);
        F2(ve,ho+1) = F2(ve, ho) + dF2deAvg * dXi + SF2pre(ve,ho+1);
        F3(ve,ho+1) = F3(ve, ho) + dF3deAvg * dXi + SF3pre(ve,ho+1);
        F4(ve,ho+1) = F4(ve, ho) + dF4deAvg * dXi + SF4pre(ve,ho+1);
        
        % Magnitudes fisicas
        A = F3(ve, ho+1) ^2 / 2 / F1(ve, ho+1) - F4(ve, ho+1);
        B = gamma / (gamma - 1) * F1(ve, ho+1) * F2(ve, ho+1);
        C = - (gamma + 1) / 2 / (gamma - 1) * F1(ve, ho+1)^3;
        ro(ve,ho+1) = (-B + sqrt(B^2 - 4 * A * C)) / 2 / A;
        
        u(ve,ho+1) = F1(ve,ho+1) / ro(ve,ho+1);
        v(ve,ho+1) = F3(ve,ho+1) / F1(ve,ho+1);
        p(ve,ho+1) = F2(ve,ho+1) - F1(ve,ho+1)*u(ve,ho+1);
        T(ve,ho+1) = p(ve,ho+1) / ro(ve,ho+1) / R;
        M(ve,ho+1) = sqrt(v(ve,ho+1)^2 + u(ve,ho+1)^2) / sqrt(gamma * p(ve,ho+1) / ro(ve,ho+1));
        
        end
        
        % MACH NUMBER BOUNDARY.
        if ve == 1
            
            if xi > E
                phi = theta - atan2(abs(v(ve,ho)) , u(ve,ho));
            else
                phi = atan2(v(ve,ho+1) , u(ve,ho+1));
            end
            
            f_cal = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (M(ve,ho+1)^2 - 1))) - atan(sqrt(M(ve,ho+1)^2 - 1));
            
            f_act = f_cal + phi;
            
            var = 10;
            tol = 1e-3;
            MAct(ve,ho+1) = 1;
            guessed = 10e10;
            while (abs(f_act - guessed) > tol)
                guessed = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (MAct(ve,ho+1)^2 - 1))) - atan(sqrt(MAct(ve,ho+1)^2 - 1));
                MAct(ve,ho+1) = MAct(ve,ho+1) + 1e-5;
            end
            
            pAct(ve,ho+1) = p(ve,ho+1) * ((1 + ((gamma - 1) / 2) * M(ve,ho+1)^ 2) / (1 + ((gamma - 1) / 2) * MAct(ve,ho+1)^ 2))^(gamma / (gamma - 1));
            TAct(ve,ho+1) = T(ve,ho+1) * ((1 + ((gamma - 1) / 2) * M(ve,ho+1)^ 2) / (1 + ((gamma - 1) / 2) * MAct(ve,ho+1)^ 2));
            roAct(ve, ho+1) = pAct(ve,ho+1) / R / TAct(ve,ho+1);
            
            p(ve,ho+1) = pAct(ve,ho+1);
            T(ve,ho+1) = TAct(ve,ho+1);
            ro(ve,ho+1) = roAct(ve,ho+1);
            M(ve, ho+1) = MAct(ve, ho + 1);
        end
        
        % New G
        G1(ve,ho+1) = ro(ve,ho+1) * F3(ve,ho+1) / F1(ve,ho+1);
        G2(ve,ho+1) = F3(ve,ho+1);
        G3(ve,ho+1) = ro(ve,ho+1) * (F3(ve,ho+1) / F1(ve,ho+1))^2 + F2(ve,ho+1) - F1(ve,ho+1)^2 / ro(ve,ho+1);
        G4(ve,ho+1) = gamma / (gamma - 1) * (F2(ve,ho+1) - F1(ve,ho+1)^2 / ro(ve,ho+1)) * F3(ve,ho+1) / F1(ve,ho+1) + ro(ve,ho+1) / 2 * F3(ve,ho+1) / F1(ve,ho+1) * ((F1(ve,ho+1)/ro(ve,ho+1))^2 + (F3(ve,ho+1) / F1(ve,ho+1))^2);
        
        
    end
    ho = ho + 1;
    xi = xi + dXi;
    test(end + 1) = dXi;
end

x = 1:length(M);
y = 1:j;
[X, Y] = meshgrid(x, y);
pcolor(X, Y, F4);