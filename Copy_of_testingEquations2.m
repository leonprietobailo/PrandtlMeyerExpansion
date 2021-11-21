clear;
i = 50;
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

pPre = p;
n = 0;
xi = 0;
%for ho=1:i-1
ho = 1;

while(xi < xMax)
    for ve=1:j
        % Init Vars
        [dEtadX, dXi, dEta] = initVars(j, xi, E, H, theta, ve);
        
        if(ve ~= j)
            % "Forward differences"
            dF1de(ve, ho) = dEtadX * (F1(ve, ho) - F1(ve+1,ho)/ dEta) + 1 / h * ((G1(ve, ho) - G1(ve+1,ho))/dEta);
            dF2de(ve, ho) = dEtadX * (F2(ve, ho) - F2(ve+1,ho)/ dEta) + 1 / h * ((G2(ve, ho) - G2(ve+1,ho))/dEta);
            dF3de(ve, ho) = dEtadX * (F3(ve, ho) - F3(ve+1,ho)/ dEta) + 1 / h * ((G3(ve, ho) - G3(ve+1,ho))/dEta);
            dF4de(ve, ho) = dEtadX * (F4(ve, ho) - F4(ve+1,ho)/ dEta) + 1 / h * ((G4(ve, ho) - G4(ve+1,ho))/dEta);
            
            % Predicted F
            F1pre(ve,ho+1) = F1(ve, ho) + dF1de(ve,ho) * dXi;
            F2pre(ve,ho+1) = F2(ve, ho) + dF2de(ve,ho) * dXi;
            F3pre(ve,ho+1) = F3(ve, ho) + dF3de(ve,ho) * dXi;
            F4pre(ve,ho+1) = F4(ve, ho) + dF4de(ve,ho) * dXi;
            
            % Predicted values:
            A = F3pre(ve,ho+1) ^2 / 2 / F1pre(ve,ho+1) - F4pre(ve,ho+1);
            B = gamma / (gamma - 1) * F1pre(ve,ho+1) * F2pre(ve,ho+1);
            C = - (gamma + 1) / 2 / (gamma - 1) * F1pre(ve,ho+1)^3;
            predictedRo = (-B + sqrt(B^2 - 4 * A * C)) / 2 / A;
            
            % Predicted G
            G1pre(ve,ho+1) = predictedRo * F3pre(ve,ho+1) / F1pre(ve,ho+1);
            G2pre(ve,ho+1) = F3pre(ve,ho+1);
            G3pre(ve,ho+1) = predictedRo * (F3pre(ve,ho+1) / F1pre(ve,ho+1))^2 + F2pre(ve,ho+1) - F1pre(ve,ho+1)^2 / predictedRo;
            G4pre(ve,ho+1) = gamma / (gamma - 1) * (F2pre(ve,ho+1) - F1pre(ve,ho+1)^2 / predictedRo) * F3pre(ve,ho+1) / F1pre(ve,ho+1) + predictedRo / 2 * F3pre(ve,ho+1) / F1pre(ve,ho+1) * ((F1pre(ve,ho+1)/predictedRo)^2 + (F3pre(ve,ho+1) / F1pre(ve,ho+1))^2);
            
            %         else % TOP LAYER
            %             dF1de(ve, ho) = dEtadX * (F1(ve - 1, ho) - F1(ve,ho)/ dEta) +
            %             1 / h * ((G1(ve - 1,ho) - G1(ve,ho))/dEta); dF2de(ve, ho) =
            %             dEtadX * (F2(ve - 1, ho) - F2(ve,ho)/ dEta) + 1 / h * ((G2(ve
            %             - 1,ho) - G2(ve,ho))/dEta); dF3de(ve, ho) = dEtadX * (F3(ve -
            %             1, ho) - F3(ve,ho)/ dEta) + 1 / h * ((G3(ve - 1,ho) -
            %             G3(ve,ho))/dEta); dF4de(ve, ho) = dEtadX * (F4(ve - 1, ho) -
            %             F4(ve,ho)/ dEta) + 1 / h * ((G4(ve - 1,ho) -
            %             G4(ve,ho))/dEta);
            %
            %             % Seguramente sobre estructura matriz. F1pre(ve,ho+1) =
            %             F1(ve, ho) + dF1de(ve, ho) * dXi; F2pre(ve,ho+1) = F2(ve, ho)
            %             + dF2de(ve, ho) * dXi; F3pre(ve,ho+1) = F3(ve, ho) +
            %             dF3de(ve, ho) * dXi; F4pre(ve,ho+1) = F4(ve, ho) + dF4de(ve,
            %             ho) * dXi;
            %
            %             % Predicted values: A = F3pre(ve,ho+1) ^2 / 2 /
            %             F1pre(ve,ho+1) - F4pre(ve,ho+1); B = gamma / (gamma - 1) *
            %             F1pre(ve,ho+1) * F2pre(ve,ho+1); C = - (gamma + 1) / 2 /
            %             (gamma - 1) * F1pre(ve,ho+1)^3; predictedRo = (-B + sqrt(B^2
            %             - 4 * A * C)) / 2 / A;
            %
            %             % Predicted G G1pre(ve,ho+1) = predictedRo * F3pre(ve,ho+1) /
            %             F1pre(ve,ho+1); G2pre(ve,ho+1) = F3pre(ve, ho+1);
            %             G3pre(ve,ho+1) = predictedRo * (F3pre(ve,ho+1) /
            %             F1pre(ve,ho+1))^2 + F2pre(ve,ho+1) - F1pre(ve,ho+1)^2 /
            %             predictedRo; G4pre(ve,ho+1) = gamma / (gamma - 1) *
            %             (F2pre(ve,ho+1) - F1pre(ve,ho+1)^2 / predictedRo) *
            %             F3pre(ve,ho+1) / F1pre(ve,ho+1) + predictedRo / 2 *
            %             F3pre(ve,ho+1) / F1pre(ve,ho+1) *
            %             ((F1pre(ve,ho+1)/predictedRo)^2 + (F3pre(ve,ho+1) /
            %             F1pre(ve,ho+1))^2);
        end
    end
    
    for ve=1:j
        % Init Vars
        [dEtadX, dXi, dEta] = initVars(j, xi, E, H, theta, ve);
        
        if(ve ~= j)
            % Predicted dFde
            dF1dePre(ve,ho+1) = dEtadX * (F1pre(ve-1, ho+1) - F1pre(ve,ho+1)/ dEta) + 1 / h * ((G1pre(ve-1, ho+1) - G1pre(ve,ho+1))/dEta);
            dF2dePre(ve,ho+1) = dEtadX * (F2pre(ve-1, ho+1) - F2pre(ve,ho+1)/ dEta) + 1 / h * ((G2pre(ve-1, ho+1) - G2pre(ve,ho+1))/dEta);
            dF3dePre(ve,ho+1) = dEtadX * (F3pre(ve-1, ho+1) - F3pre(ve,ho+1)/ dEta) + 1 / h * ((G3pre(ve-1, ho+1) - G3pre(ve,ho+1))/dEta);
            dF4dePre(ve,ho+1) = dEtadX * (F4pre(ve-1, ho+1) - F4pre(ve,ho+1)/ dEta) + 1 / h * ((G4pre(ve-1, ho+1) - G4pre(ve,ho+1))/dEta);
            
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
            
%         else
%             % Predicted dFde
%             dF1dePre(ve,ho+1) = dEtadX * (F1pre(ve - 1, ho+1) - F1pre(ve, ho + 1)/ dEta) + 1 / h * ((G1pre(ve - 1, ho+1) - G1pre(ve, ho + 1))/dEta);
%             dF2dePre(ve,ho+1) = dEtadX * (F2pre(ve - 1, ho+1) - F2pre(ve, ho + 1)/ dEta) + 1 / h * ((G2pre(ve - 1, ho+1) - G2pre(ve, ho + 1))/dEta);
%             dF3dePre(ve,ho+1) = dEtadX * (F3pre(ve - 1, ho+1) - F3pre(ve, ho + 1)/ dEta) + 1 / h * ((G3pre(ve - 1, ho+1) - G3pre(ve, ho + 1))/dEta);
%             dF4dePre(ve,ho+1) = dEtadX * (F4pre(ve - 1, ho+1) - F4pre(ve, ho + 1)/ dEta) + 1 / h * ((G4pre(ve - 1, ho+1) - G4pre(ve, ho + 1))/dEta);
%             
%             % Averaging
%             dF1deAvg = (dF1de(ve, ho) + dF1dePre(ve,ho+1)) / 2;
%             dF2deAvg = (dF2de(ve, ho) + dF2dePre(ve,ho+1)) / 2;
%             dF3deAvg = (dF3de(ve, ho) + dF3dePre(ve,ho+1)) / 2;
%             dF4deAvg = (dF4de(ve, ho) + dF4dePre(ve,ho+1)) / 2;
%             
%             % New F
%             F1(ve,ho+1) = F1(ve, ho) + dF1deAvg * dXi;
%             F2(ve,ho+1) = F2(ve, ho) + dF2deAvg * dXi;
%             F3(ve,ho+1) = F3(ve, ho) + dF3deAvg * dXi;
%             F4(ve,ho+1) = F4(ve, ho) + dF4deAvg * dXi;
%             
%             % Magnitudes fisicas
%             A = F3(ve, ho+1) ^2 / 2 / F1(ve, ho+1) - F4(ve, ho+1);
%             B = gamma / (gamma - 1) * F1(ve, ho+1) * F2(ve, ho+1);
%             C = - (gamma + 1) / 2 / (gamma - 1) * F1(ve, ho+1)^3;
%             ro(ve,ho+1) = (-B + sqrt(B^2 - 4 * A * C)) / 2 / A;
%             
%             u(ve,ho+1) = F1(ve,ho+1) / ro(ve,ho+1);
%             v(ve,ho+1) = F3(ve,ho+1) / F1(ve,ho+1);
%             p(ve,ho+1) = F2(ve,ho+1) - F1(ve,ho+1)*u(ve,ho+1);
%             T(ve,ho+1) = p(ve,ho+1) / ro(ve,ho+1) / R;
%             M(ve,ho+1) = sqrt(v(ve,ho+1)^2 + u(ve,ho+1)^2) / sqrt(gamma * p(ve,ho+1) / ro(ve,ho+1));
        end
    end
    
    for ve=1:j
        % Init Vars
        [dEtadX, dXi, dEta] = initVars(j, xi, E, H, theta, ve);
        
        if(ve ~= 1 && ve ~= j)
            % Viscosity
            SF1(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho) * (F1(ve - 1,ho) - 2 * F1(ve,ho) + F1(ve + 1, ho));
            SF2(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho) * (F2(ve - 1,ho) - 2 * F2(ve,ho) + F2(ve + 1, ho));
            SF3(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho) * (F3(ve - 1,ho) - 2 * F3(ve,ho) + F3(ve + 1, ho));
            SF4(ve,ho) = Cy * abs(p(ve + 1, ho) - 2 * p(ve, ho) + p(ve - 1,ho)) / p(ve + 1, ho) + 2 * p(ve, ho) + p(ve - 1,ho) * (F4(ve - 1,ho) - 2 * F4(ve,ho) + F4(ve + 1, ho));
            
            % Prdicted F
            F1pre(ve,ho+1) = F1(ve, ho) + dF1de(ve, ho) * dXi + SF1(ve, ho);
            F2pre(ve,ho+1) = F2(ve, ho) + dF2de(ve, ho) * dXi + SF2(ve, ho);
            F3pre(ve,ho+1) = F3(ve, ho) + dF3de(ve, ho) * dXi + SF3(ve, ho);
            F4pre(ve,ho+1) = F4(ve, ho) + dF4de(ve, ho) * dXi + SF4(ve, ho);
            pPre(ve,ho+1) = F2pre(ve,ho+1) - F1pre(ve,ho+1)^2 / predictedRo;
        end
    end
    
    for ve =1:j
        % Init Vars
        [dEtadX, dXi, dEta] = initVars(j, xi, E, H, theta, ve);
        
        % Predicted SF
        SF1pre(ve,ho) = Cy * abs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1,ho)) / pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1,ho) * (F1(ve - 1,ho) - 2 * F1(ve,ho) + F1(ve + 1, ho));
        SF2pre(ve,ho) = Cy * abs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1,ho)) / pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1,ho) * (F2(ve - 1,ho) - 2 * F2(ve,ho) + F2(ve + 1, ho));
        SF3pre(ve,ho) = Cy * abs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1,ho)) / pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1,ho) * (F3(ve - 1,ho) - 2 * F3(ve,ho) + F3(ve + 1, ho));
        SF4pre(ve,ho) = Cy * abs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1,ho)) / pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1,ho) * (F4(ve - 1,ho) - 2 * F4(ve,ho) + F4(ve + 1, ho));
        
        % New F
        F1(ve,ho+1) = F1(ve, ho) + dF1deAvg * dXi + SF1pre(ve,ho+1);
        F2(ve,ho+1) = F2(ve, ho) + dF2deAvg * dXi + SF2pre(ve,ho+1);
        F3(ve,ho+1) = F3(ve, ho) + dF3deAvg * dXi + SF3pre(ve,ho+1);
        F4(ve,ho+1) = F4(ve, ho) + dF4deAvg * dXi + SF4pre(ve,ho+1);
        
        % New G
        G1(ve,ho+1) = ro(ve,ho+1) * F3(ve,ho+1) / F1(ve,ho+1);
        G2(ve,ho+1) = F3(ve,ho+1);
        G3(ve,ho+1) = ro(ve,ho+1) * (F3(ve,ho+1) / F1(ve,ho+1))^2 + F2(ve,ho+1) - F1(ve,ho+1)^2 / ro(ve,ho+1);
        G4(ve,ho+1) = gamma / (gamma - 1) * (F2(ve,ho+1) - F1(ve,ho+1)^2 / ro(ve,ho+1)) * F3(ve,ho+1) / F1(ve,ho+1) + ro(ve,ho+1) / 2 * F3(ve,ho+1) / F1(ve,ho+1) * ((F1(ve,ho+1)/ro(ve,ho+1))^2 + (F3(ve,ho+1) / F1(ve,ho+1))^2);
    
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
    ho = ho + 1;
    xi = xi + dXi;
end

% x = 0:xMax/i:(xMax-xMax / i);
% % y = 0:xMax/i:(xMax-xMax / i);
% %x = 1:i;
% y = 1:j;
% [X, Y] = meshgrid(x, y);
% pcolor(xV, Y, M);