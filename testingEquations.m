i = 10;
j = 10;
gamma = 1.4;
xMax = 35;
E = 10;
hMin = 40;
theta = 5.352*pi/180;


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

for horizontal=1:i-1
   for vertical=1:j-1
       dx = xMax / i;
       x = dx * horizontal;
       if x < E
           ys = 0;
           h = hMin;
       else
           h = hMin + (x - E) * tan(theta);
           ys = -(x - E) * tan(theta);
       end
       dy = h / j;
       y = dy * vertical;
       eta = (y - ys) / h;
       dEta = 1 / j;
       
       if x < E
           dEtadX = 0;
       else
           dEtadX = (1 - eta) * tan(theta) / h;
       end
       dF1de = dEtadX * (F1(horizontal, vertical) - F1(horizontal, vertical + 1)/ dEta) + 1 / h * ((G1(horizontal, vertical) - G1(horizontal, vertical+1))/dEta);
       dF2de = dEtadX * (F2(horizontal, vertical) - F2(horizontal, vertical + 1)/ dEta) + 1 / h * ((G2(horizontal, vertical) - G2(horizontal, vertical+1))/dEta);
       dF3de = dEtadX * (F3(horizontal, vertical) - F3(horizontal, vertical + 1)/ dEta) + 1 / h * ((G3(horizontal, vertical) - G3(horizontal, vertical+1))/dEta);
       dF4de = dEtadX * (F4(horizontal, vertical) - F4(horizontal, vertical + 1)/ dEta) + 1 / h * ((G4(horizontal, vertical) - G4(horizontal, vertical+1))/dEta);

       F1(horizontal + 1, vertical) = F1(horizontal, vertical) + dF1de * 
   end
end