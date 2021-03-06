clear;
close all;

j = 41;
gamma = 1.4;
xMax = 65;
E = 10;
H = 40;
theta = 5.352*pi/180;
R = 287;
Cy = 0.5;
Courant = 0.5;
% close all;

base(:,1) = ones(j,1);

% 
% vVal = 0;
% roVal = 1.225;
% pVal = 101325;
% TVal = pVal/(roVal*R);
% MVal = 2;
% a_in = sqrt(gamma*R*TVal); % delete
% uVal = MVal * a_in;

uVal = .678e3;
vVal = 0;
roVal = 1.23;
pVal = .101e6;
TVal = .2861e3;
MVal = 2;

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

yP = (0:H/(j-1):H)';
x = zeros(j,1);

while(xi < xMax)
    %% Marching Step
    if(xi > E)
        tv = 1;
    end
        [dEtadX, dEta, h, dy, yP(:, ho + 1)] = initVars(j, xi, E, H, theta);
        [dXi] = computeStep(M, theta, ho, Courant, dy);
    %% PREDICTOR STEP
    for ve=1:j
        if(ve ~= 1 && ve ~= j)
            % Fowrard differences + Viscosity.
            dF1de(ve, ho) = dEtadX(ve) * ((F1(ve, ho) - F1(ve+1,ho))/ dEta) + 1 / h * ((G1(ve, ho) - G1(ve+1,ho))/dEta);
            dF2de(ve, ho) = dEtadX(ve) * ((F2(ve, ho) - F2(ve+1,ho))/ dEta) + 1 / h * ((G2(ve, ho) - G2(ve+1,ho))/dEta);
            dF3de(ve, ho) = dEtadX(ve) * ((F3(ve, ho) - F3(ve+1,ho))/ dEta) + 1 / h * ((G3(ve, ho) - G3(ve+1,ho))/dEta);
            dF4de(ve, ho) = dEtadX(ve) * ((F4(ve, ho) - F4(ve+1,ho))/ dEta) + 1 / h * ((G4(ve, ho) - G4(ve+1,ho))/dEta);
            
            SF1(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / (p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho)) * (F1(ve - 1,ho) - 2 * F1(ve,ho) + F1(ve + 1, ho));
            SF2(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / (p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho)) * (F2(ve - 1,ho) - 2 * F2(ve,ho) + F2(ve + 1, ho));
            SF3(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / (p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho)) * (F3(ve - 1,ho) - 2 * F3(ve,ho) + F3(ve + 1, ho));
            SF4(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / (p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho)) * (F4(ve - 1,ho) - 2 * F4(ve,ho) + F4(ve + 1, ho));
            
            F1pre(ve,ho+1) = F1(ve, ho) + dF1de(ve,ho) * dXi + SF1(ve,ho);
            F2pre(ve,ho+1) = F2(ve, ho) + dF2de(ve,ho) * dXi + SF2(ve,ho);
            F3pre(ve,ho+1) = F3(ve, ho) + dF3de(ve,ho) * dXi + SF3(ve,ho);
            F4pre(ve,ho+1) = F4(ve, ho) + dF4de(ve,ho) * dXi + SF4(ve,ho);
            
        elseif ve == 1
            % Forward differences.
            asd = ((F1(ve, ho) - F1(ve+1,ho)));
            dF1de(ve, ho) = dEtadX(ve) * ((F1(ve, ho) - F1(ve+1,ho))/ dEta) + 1 / h * ((G1(ve, ho) - G1(ve+1,ho))/dEta);
            dF2de(ve, ho) = dEtadX(ve) * ((F2(ve, ho) - F2(ve+1,ho))/ dEta) + 1 / h * ((G2(ve, ho) - G2(ve+1,ho))/dEta);
            dF3de(ve, ho) = dEtadX(ve) * ((F3(ve, ho) - F3(ve+1,ho))/ dEta) + 1 / h * ((G3(ve, ho) - G3(ve+1,ho))/dEta);
            dF4de(ve, ho) = dEtadX(ve) * ((F4(ve, ho) - F4(ve+1,ho))/ dEta) + 1 / h * ((G4(ve, ho) - G4(ve+1,ho))/dEta);
            
            F1pre(ve,ho+1) = F1(ve, ho) + dF1de(ve,ho) * dXi;
            F2pre(ve,ho+1) = F2(ve, ho) + dF2de(ve,ho) * dXi;
            F3pre(ve,ho+1) = F3(ve, ho) + dF3de(ve,ho) * dXi;
            F4pre(ve,ho+1) = F4(ve, ho) + dF4de(ve,ho) * dXi;
            
        else
            % Rearwards differences.
            dF1de(ve, ho) = dEtadX(ve) * ((F1(ve - 1, ho) - F1(ve,ho))/ dEta) + 1 / h * ((G1(ve - 1,ho) - G1(ve,ho))/dEta);
            dF2de(ve, ho) = dEtadX(ve) * ((F2(ve - 1, ho) - F2(ve,ho))/ dEta) + 1 / h * ((G2(ve - 1,ho) - G2(ve,ho))/dEta);
            dF3de(ve, ho) = dEtadX(ve) * ((F3(ve - 1, ho) - F3(ve,ho))/ dEta) + 1 / h * ((G3(ve - 1,ho) - G3(ve,ho))/dEta);
            dF4de(ve, ho) = dEtadX(ve) * ((F4(ve - 1, ho) - F4(ve,ho))/ dEta) + 1 / h * ((G4(ve - 1,ho) - G4(ve,ho))/dEta);
            
            F1pre(ve,ho+1) = F1(ve, ho) + dF1de(ve,ho) * dXi;
            F2pre(ve,ho+1) = F2(ve, ho) + dF2de(ve,ho) * dXi;
            F3pre(ve,ho+1) = F3(ve, ho) + dF3de(ve,ho) * dXi;
            F4pre(ve,ho+1) = F4(ve, ho) + dF4de(ve,ho) * dXi;
        end
        
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
        
        pPre(ve,ho+1) = F2pre(ve,ho+1) - F1pre(ve,ho+1)^2 / roPre(ve, ho + 1);
        
    end
    
    %% CORRECTOR STEP
    for ve=1:j
        if(ve ~= 1 && ve ~= j)
            % Predicted dFde [TOP INFLOW]
            dF1dePre(ve,ho+1) = dEtadX(ve) * ((F1pre(ve-1, ho+1) - F1pre(ve,ho+1))/ dEta) + 1 / h * ((G1pre(ve-1, ho+1) - G1pre(ve,ho+1))/dEta);
            dF2dePre(ve,ho+1) = dEtadX(ve) * ((F2pre(ve-1, ho+1) - F2pre(ve,ho+1))/ dEta) + 1 / h * ((G2pre(ve-1, ho+1) - G2pre(ve,ho+1))/dEta);
            dF3dePre(ve,ho+1) = dEtadX(ve) * ((F3pre(ve-1, ho+1) - F3pre(ve,ho+1))/ dEta) + 1 / h * ((G3pre(ve-1, ho+1) - G3pre(ve,ho+1))/dEta);
            dF4dePre(ve,ho+1) = dEtadX(ve) * ((F4pre(ve-1, ho+1) - F4pre(ve,ho+1))/ dEta) + 1 / h * ((G4pre(ve-1, ho+1) - G4pre(ve,ho+1))/dEta);
            
            dF1deAvg = (dF1de(ve, ho) + dF1dePre(ve,ho+1)) / 2;
            dF2deAvg = (dF2de(ve, ho) + dF2dePre(ve,ho+1)) / 2;
            dF3deAvg = (dF3de(ve, ho) + dF3dePre(ve,ho+1)) / 2;
            dF4deAvg = (dF4de(ve, ho) + dF4dePre(ve,ho+1)) / 2;
            
            SF1pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho+1) - 2 * pPre(ve, ho+1) + pPre(ve - 1,ho+1)) / (pPre(ve + 1, ho+1) + 2 * pPre(ve, ho+1) + pPre(ve - 1,ho+1)) * (F1pre(ve - 1,ho+1) - 2 * F1pre(ve,ho+1) + F1pre(ve + 1, ho+1));
            SF2pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho+1) - 2 * pPre(ve, ho+1) + pPre(ve - 1,ho+1)) / (pPre(ve + 1, ho+1) + 2 * pPre(ve, ho+1) + pPre(ve - 1,ho+1)) * (F2pre(ve - 1,ho+1) - 2 * F2pre(ve,ho+1) + F2pre(ve + 1, ho+1));
            SF3pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho+1) - 2 * pPre(ve, ho+1) + pPre(ve - 1,ho+1)) / (pPre(ve + 1, ho+1) + 2 * pPre(ve, ho+1) + pPre(ve - 1,ho+1)) * (F3pre(ve - 1,ho+1) - 2 * F3pre(ve,ho+1) + F3pre(ve + 1, ho+1));
            SF4pre(ve,ho+1) = Cy * abs(pPre(ve + 1, ho+1) - 2 * pPre(ve, ho+1) + pPre(ve - 1,ho+1)) / (pPre(ve + 1, ho+1) + 2 * pPre(ve, ho+1) + pPre(ve - 1,ho+1)) * (F4pre(ve - 1,ho+1) - 2 * F4pre(ve,ho+1) + F4pre(ve + 1, ho+1));
            
            F1(ve,ho+1) = F1(ve, ho) + dF1deAvg * dXi + SF1pre(ve,ho+1);
            F2(ve,ho+1) = F2(ve, ho) + dF2deAvg * dXi + SF2pre(ve,ho+1);
            F3(ve,ho+1) = F3(ve, ho) + dF3deAvg * dXi + SF3pre(ve,ho+1);
            F4(ve,ho+1) = F4(ve, ho) + dF4deAvg * dXi + SF4pre(ve,ho+1);
            
        elseif (ve == 1)
            % Predicted dFde [BOT]
            dF1dePre(ve,ho+1) = dEtadX(ve) * ((F1pre(ve, ho+1) - F1pre(ve + 1, ho + 1))/ dEta) + 1 / h * ((G1pre(ve, ho+1) - G1pre(ve + 1, ho + 1))/dEta);
            dF2dePre(ve,ho+1) = dEtadX(ve) * ((F2pre(ve, ho+1) - F2pre(ve + 1, ho + 1))/ dEta) + 1 / h * ((G2pre(ve, ho+1) - G2pre(ve + 1, ho + 1))/dEta);
            dF3dePre(ve,ho+1) = dEtadX(ve) * ((F3pre(ve, ho+1) - F3pre(ve + 1, ho + 1))/ dEta) + 1 / h * ((G3pre(ve, ho+1) - G3pre(ve + 1, ho + 1))/dEta);
            dF4dePre(ve,ho+1) = dEtadX(ve) * ((F4pre(ve, ho+1) - F4pre(ve + 1, ho + 1))/ dEta) + 1 / h * ((G4pre(ve, ho+1) - G4pre(ve + 1, ho + 1))/dEta);
            
            dF1deAvg = (dF1de(ve, ho) + dF1dePre(ve,ho+1)) / 2;
            dF2deAvg = (dF2de(ve, ho) + dF2dePre(ve,ho+1)) / 2;
            dF3deAvg = (dF3de(ve, ho) + dF3dePre(ve,ho+1)) / 2;
            dF4deAvg = (dF4de(ve, ho) + dF4dePre(ve,ho+1)) / 2;
            
            F1(ve,ho+1) = F1(ve, ho) + dF1deAvg * dXi;
            F2(ve,ho+1) = F2(ve, ho) + dF2deAvg * dXi;
            F3(ve,ho+1) = F3(ve, ho) + dF3deAvg * dXi;
            F4(ve,ho+1) = F4(ve, ho) + dF4deAvg * dXi;
        else
            dF1dePre(ve,ho+1) = dEtadX(ve) * ((F1pre(ve-1, ho+1) - F1pre(ve,ho+1))/ dEta) + 1 / h * ((G1pre(ve-1, ho+1) - G1pre(ve,ho+1))/dEta);
            dF2dePre(ve,ho+1) = dEtadX(ve) * ((F2pre(ve-1, ho+1) - F2pre(ve,ho+1))/ dEta) + 1 / h * ((G2pre(ve-1, ho+1) - G2pre(ve,ho+1))/dEta);
            dF3dePre(ve,ho+1) = dEtadX(ve) * ((F3pre(ve-1, ho+1) - F3pre(ve,ho+1))/ dEta) + 1 / h * ((G3pre(ve-1, ho+1) - G3pre(ve,ho+1))/dEta);
            dF4dePre(ve,ho+1) = dEtadX(ve) * ((F4pre(ve-1, ho+1) - F4pre(ve,ho+1))/ dEta) + 1 / h * ((G4pre(ve-1, ho+1) - G4pre(ve,ho+1))/dEta);
            
            dF1deAvg = (dF1de(ve, ho) + dF1dePre(ve,ho+1)) / 2;
            dF2deAvg = (dF2de(ve, ho) + dF2dePre(ve,ho+1)) / 2;
            dF3deAvg = (dF3de(ve, ho) + dF3dePre(ve,ho+1)) / 2;
            dF4deAvg = (dF4de(ve, ho) + dF4dePre(ve,ho+1)) / 2;
            
            F1(ve,ho+1) = F1(ve, ho) + dF1deAvg * dXi;
            F2(ve,ho+1) = F2(ve, ho) + dF2deAvg * dXi;
            F3(ve,ho+1) = F3(ve, ho) + dF3deAvg * dXi;
            F4(ve,ho+1) = F4(ve, ho) + dF4deAvg * dXi;
        end
        
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
        
        if ve == 1
            
            if xi > E
                phi = theta - atan2(abs(v(ve,ho+1)) , u(ve,ho+1));
            else
                phi = atan2(v(ve,ho+1) , u(ve,ho+1));
            end
            
            f_cal = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (M(ve,ho+1)^2 - 1))) - atan(sqrt(M(ve,ho+1)^2 - 1));
            
            f_act = f_cal + phi;
            
            var = 10;
            tol = 1e-3;
            guessed = 10e10;
            precision = 0.0000001;
            %PROVISIONA, CAMBIAR
            a_int = 1.1;
            b_int = 2.9;
            zero_f1 = sqrt((gamma + 1)/(gamma - 1))*(atan(sqrt(((gamma - 1)/(gamma + 1))*(a_int^2 - 1)))) - (atan(sqrt((a_int^2) - 1))) - f_act; % Function used to find its zero
            zero_f2 = sqrt((gamma + 1)/(gamma - 1))*(atan(sqrt(((gamma - 1)/(gamma + 1))*(((a_int + b_int)/2)^2 - 1)))) - (atan(sqrt((((a_int + b_int)/2)^2) - 1))) - f_act;
            while ((b_int-a_int)/2 > precision)
                if (zero_f1*zero_f2 <=0)
                    b_int = (a_int + b_int)/2;
                else
                    a_int = (a_int + b_int)/2;8
                end
                zero_f1 = sqrt((gamma + 1)/(gamma - 1))*(atan(sqrt(((gamma - 1)/(gamma + 1))*(a_int^2 - 1)))) - (atan(sqrt((a_int^2) - 1))) - f_act;
                zero_f2 = sqrt((gamma + 1)/(gamma - 1))*(atan(sqrt(((gamma - 1)/(gamma + 1))*(((a_int + b_int)/2)^2 - 1)))) - (atan(sqrt((((a_int + b_int)/2)^2) - 1))) - f_act;
            end
            M_act = (a_int + b_int)/2;
            
            pAct(ve,ho+1) = p(ve,ho+1) * ((1 + ((gamma - 1) / 2) * M_act^2) / (1 + ((gamma - 1) / 2) * M_act^2))^(gamma / (gamma - 1));
            TAct(ve,ho+1) = T(ve,ho+1) * ((1 + ((gamma - 1) / 2) * M_act^2) / (1 + ((gamma - 1) / 2) * M_act^2));
            roAct(ve, ho+1) = pAct(ve,ho+1) / R / TAct(ve,ho+1);
            
            p(ve,ho+1) = pAct(ve,ho+1);
            T(ve,ho+1) = TAct(ve,ho+1);
            ro(ve,ho+1) = roAct(ve,ho+1);
            M(ve, ho+1) = M_act;
            
            if (xi > E)
                v(ve,ho+1) = -(u(ve,ho+1)*tan(theta));
            else
                v(ve,ho+1) = 0;
            end
            
            F1c(ve,ho+1) = ro(ve,ho+1) * u(ve,ho+1);
            F2c(ve,ho+1) = ro(ve,ho+1) * u(ve,ho+1) ^ 2 + p(ve,ho+1);
            F3c(ve,ho+1) = ro(ve,ho+1) * u(ve,ho+1) * v(ve,ho+1);
            F4c(ve,ho+1) = gamma / (gamma - 1) * p(ve,ho+1)*u(ve,ho+1) + ro(ve,ho+1)*u(ve,ho+1)*(u(ve,ho+1)^2 + v(ve,ho+1)^2) / 2;
            
            F1(ve,ho+1) = F1c(ve,ho+1);
            F2(ve,ho+1) = F2c(ve,ho+1);
            F3(ve,ho+1) = F3c(ve,ho+1);
            F4(ve,ho+1) = F4c(ve,ho+1);
        
%             G1(ve,ho+1) = ro(ve,ho+1) * F3c(ve,ho+1) / F1c(ve,ho+1);
%             G2(ve,ho+1) = F3c(ve,ho+1);
%             G3(ve,ho+1) = ro(ve,ho+1) * (F3c(ve,ho+1) / F1c(ve,ho+1))^2 + F2c(ve,ho+1) - F1c(ve,ho+1)^2 / ro(ve,ho+1);
%             G4(ve,ho+1) = gamma / (gamma - 1) * (F2c(ve,ho+1) - F1c(ve,ho+1)^2 / ro(ve,ho+1)) * F3c(ve,ho+1) / F1c(ve,ho+1) + ro(ve,ho+1) / 2 * F3c(ve,ho+1) / F1c(ve,ho+1) * ((F1c(ve,ho+1)/ro(ve,ho+1))^2 + (F3c(ve,ho+1) / F1c(ve,ho+1))^2);
        
        else
        % New G
        G1(ve,ho+1) = ro(ve,ho+1) * F3(ve,ho+1) / F1(ve,ho+1);
        G2(ve,ho+1) = F3(ve,ho+1);
        G3(ve,ho+1) = ro(ve,ho+1) * (F3(ve,ho+1) / F1(ve,ho+1))^2 + F2(ve,ho+1) - F1(ve,ho+1)^2 / ro(ve,ho+1);
        G4(ve,ho+1) = gamma / (gamma - 1) * (F2(ve,ho+1) - F1(ve,ho+1)^2 / ro(ve,ho+1)) * F3(ve,ho+1) / F1(ve,ho+1) + ro(ve,ho+1) / 2 * F3(ve,ho+1) / F1(ve,ho+1) * ((F1(ve,ho+1)/ro(ve,ho+1))^2 + (F3(ve,ho+1) / F1(ve,ho+1))^2);
        end
        G1(ve,ho+1) = ro(ve,ho+1) * F3(ve,ho+1) / F1(ve,ho+1);
        G2(ve,ho+1) = F3(ve,ho+1);
        G3(ve,ho+1) = ro(ve,ho+1) * (F3(ve,ho+1) / F1(ve,ho+1))^2 + F2(ve,ho+1) - F1(ve,ho+1)^2 / ro(ve,ho+1);
        G4(ve,ho+1) = gamma / (gamma - 1) * (F2(ve,ho+1) - F1(ve,ho+1)^2 / ro(ve,ho+1)) * F3(ve,ho+1) / F1(ve,ho+1) + ro(ve,ho+1) / 2 * F3(ve,ho+1) / F1(ve,ho+1) * ((F1(ve,ho+1)/ro(ve,ho+1))^2 + (F3(ve,ho+1) / F1(ve,ho+1))^2);
    end
    ho = ho + 1;
    xi = xi + dXi;
    test(end + 1) = dXi;
    %x(:, end + 1) = xi * ones(j, 1);
end

% x = 1:size(u,2);
% y = 1:j;
% [X, Y] = meshgrid(x, y);
% p = pcolor(X, Y, M);
% p.EdgeAlpha = 0;


figure
subplot(2, 3, 1);
S = mesh(yP,M); % Mesh function plots in 3D, but this has no sense in this simulation. Move the axes to view it in 2d.
S.FaceColor = 'Flat';
xlabel('x steps');
zlabel('y [m]');
title('Mach Number [1]');
colorbar;
view([0 0])

subplot(2, 3, 2);
S = mesh(yP,T); % Mesh function plots in 3D, but this has no sense in this simulation. Move the axes to view it in 2d.
S.FaceColor = 'Flat';
xlabel('x steps');
zlabel('y [m]');
title('Temperature [K]');
colorbar;
view([0 0])

subplot(2, 3, 3);
S = mesh(yP,p); % Mesh function plots in 3D, but this has no sense in this simulation. Move the axes to view it in 2d.
S.FaceColor = 'Flat';
xlabel('x steps');
zlabel('y [m]');
title('Pressure [Pa]');
colorbar;
view([0 0])

subplot(2, 3, 4);
S = mesh(yP,u); % Mesh function plots in 3D, but this has no sense in this simulation. Move the axes to view it in 2d.
S.FaceColor = 'Flat';
xlabel('x steps');
zlabel('y [m]');
title('Horizontal Speed [m/s]');
colorbar;
view([0 0])

subplot(2, 3, 5);
S = mesh(yP,v); % Mesh function plots in 3D, but this has no sense in this simulation. Move the axes to view it in 2d.
S.FaceColor = 'Flat';
xlabel('x steps');
zlabel('y [m]');
title('Vertical Speed [m/s]');
colorbar;
view([0 0])

subplot(2, 3, 6);
S = mesh(yP,ro); % Mesh function plots in 3D, but this has no sense in this simulation. Move the axes to view it in 2d.
S.FaceColor = 'Flat';
xlabel('x steps');
zlabel('y [m]');
title('Density [kg/m^3]');
colorbar;
view([0 0])

sgtitle('Prandtl-Meyer Expansion Magnitudes') 