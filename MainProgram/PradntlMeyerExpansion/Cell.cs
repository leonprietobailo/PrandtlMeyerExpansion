using System;

namespace PradntlMeyerExpansion
{
    public class Cell
    {
//        double u, v, ro, p, T, M; // Magnitudes
//        double Cy;
//        int chara, ho, ve; // 0 -> Bottom, 1 -> Inflow, 2 -> Top
//        double F1, F2, F3, F4, G1, G2, G3, G4;
//        double F1pre, F2pre, F3pre, F4pre, G1pre, G2pre, G3pre, G4pre;
//        double roPre, pPre;
//        double dF1de, dF2de, dF3de, dF4de;
//        Rules r;


//        // Initial cell.
//        public Cell(Rules r)
//        {
//            u = r.getU();
//            v = r.getV();
//            ro = r.getRO();
//            p = r.getP();
//            T = r.getT();
//            M = r.getM();
//            Cy = r.getCy();
//        }

//        // Computed cell.
//        public Cell(Rules r)
//        {
//            u = 0;
//            v = 0;
//            ro = 0;
//            p = 0;
//            T = 0;
//            M = 0;
//            Cy = r.getCy();
//        }

//        public void PredictorStep(double dEtadX, double dEta, double h, double dXi, Grid g)
//        {
//            if (chara == 1)
//            {
//                double dF1de = dEtadX * ((g.GetCell(ve, ho - 1).getF(1) - g.GetCell(ve + 1, ho - 1).getF(1)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho-1).getG(1) - g.GetCell(ve + 1, ho -1).getG(1)) / dEta);
//                double dF2de = dEtadX * ((g.GetCell(ve, ho - 1).getF(2) - g.GetCell(ve + 1, ho - 1).getF(2)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(2) - g.GetCell(ve + 1, ho - 1).getG(2)) / dEta);
//                double dF3de = dEtadX * ((g.GetCell(ve, ho - 1).getF(3) - g.GetCell(ve + 1, ho - 1).getF(3)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(3) - g.GetCell(ve + 1, ho - 1).getG(3)) / dEta);
//                double dF4de = dEtadX * ((g.GetCell(ve, ho - 1).getF(4) - g.GetCell(ve + 1, ho - 1).getF(4)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(4) - g.GetCell(ve + 1, ho - 1).getG(4)) / dEta);

//                double SF1 = Cy * Math.Abs(g.GetCell(ve + 1, ho - 1).getP() - 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) / (g.GetCell(ve + 1, ho - 1).getP() + 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) * (g.GetCell(ve - 1, ho - 1).getF(1) - 2 * g.GetCell(ve, ho - 1).getF(1) + g.GetCell(ve + 1, ho - 1).getF(1));
//                double SF2 = Cy * Math.Abs(g.GetCell(ve + 1, ho - 1).getP() - 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) / (g.GetCell(ve + 1, ho - 1).getP() + 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) * (g.GetCell(ve - 1, ho - 1).getF(2) - 2 * g.GetCell(ve, ho - 1).getF(2) + g.GetCell(ve + 1, ho - 1).getF(2));
//                double SF3 = Cy * Math.Abs(g.GetCell(ve + 1, ho - 1).getP() - 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) / (g.GetCell(ve + 1, ho - 1).getP() + 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) * (g.GetCell(ve - 1, ho - 1).getF(3) - 2 * g.GetCell(ve, ho - 1).getF(3) + g.GetCell(ve + 1, ho - 1).getF(3));
//                double SF4 = Cy * Math.Abs(g.GetCell(ve + 1, ho - 1).getP() - 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) / (g.GetCell(ve + 1, ho - 1).getP() + 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) * (g.GetCell(ve - 1, ho - 1).getF(4) - 2 * g.GetCell(ve, ho - 1).getF(4) + g.GetCell(ve + 1, ho - 1).getF(4));

//                F1pre = F1 + dF1de * dXi + SF1;
//                F2pre = F2 + dF2de * dXi + SF2;
//                F3pre = F3 + dF3de * dXi + SF3;
//                F4pre = F4 + dF4de * dXi + SF4;
//            }
//            else if (chara == 0)
//            {
//                double dF1de = dEtadX * ((g.GetCell(ve, ho - 1).getF(1) - g.GetCell(ve + 1, ho - 1).getF(1)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(1) - g.GetCell(ve + 1, ho - 1).getG(1)) / dEta);
//                double dF2de = dEtadX * ((g.GetCell(ve, ho - 1).getF(2) - g.GetCell(ve + 1, ho - 1).getF(2)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(2) - g.GetCell(ve + 1, ho - 1).getG(2)) / dEta);
//                double dF3de = dEtadX * ((g.GetCell(ve, ho - 1).getF(3) - g.GetCell(ve + 1, ho - 1).getF(3)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(3) - g.GetCell(ve + 1, ho - 1).getG(3)) / dEta);
//                double dF4de = dEtadX * ((g.GetCell(ve, ho - 1).getF(4) - g.GetCell(ve + 1, ho - 1).getF(4)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(4) - g.GetCell(ve + 1, ho - 1).getG(4)) / dEta);

//                F1pre = F1 + dF1de * dXi;
//                F2pre = F2 + dF2de * dXi;
//                F3pre = F3 + dF3de * dXi;
//                F4pre = F4 + dF4de * dXi;
//            }
//            else
//            {
//                double dF1de = dEtadX * ((g.GetCell(ve - 1, ho - 1).getF(1) - g.GetCell(ve, ho - 1).getF(1)) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho - 1).getG(1) - g.GetCell(ve, ho - 1).getG(1)) / dEta);
//                double dF2de = dEtadX * ((g.GetCell(ve - 1, ho - 1).getF(2) - g.GetCell(ve, ho - 1).getF(2)) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho - 1).getG(2) - g.GetCell(ve, ho - 1).getG(2)) / dEta);
//                double dF3de = dEtadX * ((g.GetCell(ve - 1, ho - 1).getF(3) - g.GetCell(ve, ho - 1).getF(3)) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho - 1).getG(3) - g.GetCell(ve, ho - 1).getG(3)) / dEta);
//                double dF4de = dEtadX * ((g.GetCell(ve - 1, ho - 1).getF(4) - g.GetCell(ve, ho - 1).getF(4)) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho - 1).getG(4) - g.GetCell(ve, ho - 1).getG(4)) / dEta);

//                F1pre = F1 + dF1de * dXi;
//                F2pre = F2 + dF2de * dXi;
//                F3pre = F3 + dF3de * dXi;
//                F4pre = F4 + dF4de * dXi;
//            }

//            // Predicted values:
//            double A = Math.Pow(F3pre, 2) / 2 / F1pre - F4;
//            double B = r.getGamma() / (r.getGamma() - 1) * F1pre * F2pre;
//            double C = -(r.getGamma() + 1) / 2 / (r.getGamma() - 1) * Math.Pow(F1pre, 3);
//            roPre = (-B + Math.Sqrt(Math.Pow(B, 2) - 4 * A * C)) / 2 / A;


//            // Predicted G
//            G1pre = roPre * F3pre / F1pre;
//            G2pre = F3pre;
//            G3pre = roPre * Math.Pow(F3pre / F1pre, 2) + F2pre - Math.Pow(F1pre, 2) / roPre;
//            G4pre = r.getGamma() / (r.getGamma() - 1) * (F2pre - Math.Pow(F1pre, 2) / roPre) * F3pre / F1pre + roPre / 2 * F3pre / F1pre * (Math.Pow(F1pre / roPre, 2) + Math.Pow(F3pre / F1pre, 2));
//            pPre = F2pre - Math.Pow(F1pre, 2) / roPre;
//        }



//        public void CorrectorStep(double dEtadX, double dEta, double h, double dXi, Grid g)
//        {
//            if (chara == 1) {
//            // Predicted dFde[TOP INFLOW]
//                double dF1dePre = dEtadX * ((g.GetCell(ve - 1, ho + 1).getFpre(1) - g.GetCell(ve, ho + 1).getFpre(1)) / dEta) + 1 / h * ((g.GetCell(ve - 1, ho + 1).getGpre(1) - g.GetCell(ve, ho + 1).getGpre(1)) / dEta);
//                double dF2dePre = dEtadX * ((g.GetCell(ve - 1, ho + 1).getFpre(2) - g.GetCell(ve, ho + 1).getFpre(2)) / dEta) + 1 / h * ((g.GetCell(ve - 1, ho + 1).getGpre(2) - g.GetCell(ve, ho + 1).getGpre(2)) / dEta);
//                double dF3dePre = dEtadX * ((g.GetCell(ve - 1, ho + 1).getFpre(3) - g.GetCell(ve, ho + 1).getFpre(3)) / dEta) + 1 / h * ((g.GetCell(ve - 1, ho + 1).getGpre(3) - g.GetCell(ve, ho + 1).getGpre(3)) / dEta);
//                double dF4dePre = dEtadX * ((g.GetCell(ve - 1, ho + 1).getFpre(4) - g.GetCell(ve, ho + 1).getFpre(4)) / dEta) + 1 / h * ((g.GetCell(ve - 1, ho + 1).getGpre(4) - g.GetCell(ve, ho + 1).getGpre(4)) / dEta);

//                double dF1deAvg = (dF1de + dF1dePre) / 2;
//                double dF2deAvg = (dF2de + dF2dePre) / 2;
//                double dF3deAvg = (dF3de + dF3dePre) / 2;
//                double dF4deAvg = (dF4de + dF4dePre) / 2;

//                double SF1pre = Cy * Math.Sabs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1, ho)) / (pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1, ho)) * (F1pre(ve - 1, ho) - 2 * F1pre(ve, ho) + F1pre(ve + 1, ho));
//                SF2pre(ve, ho + 1) = Cy * abs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1, ho)) / (pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1, ho)) * (F2pre(ve - 1, ho) - 2 * F2pre(ve, ho) + F2pre(ve + 1, ho));
//                SF3pre(ve, ho + 1) = Cy * abs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1, ho)) / (pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1, ho)) * (F3pre(ve - 1, ho) - 2 * F3pre(ve, ho) + F3pre(ve + 1, ho));
//                SF4pre(ve, ho + 1) = Cy * abs(pPre(ve + 1, ho) - 2 * pPre(ve, ho) + pPre(ve - 1, ho)) / (pPre(ve + 1, ho) + 2 * pPre(ve, ho) + pPre(ve - 1, ho)) * (F4pre(ve - 1, ho) - 2 * F4pre(ve, ho) + F4pre(ve + 1, ho));

//                F1(ve, ho + 1) = F1(ve, ho) + dF1deAvg * dXi + SF1pre(ve, ho + 1);
//                F2(ve, ho + 1) = F2(ve, ho) + dF2deAvg * dXi + SF2pre(ve, ho + 1);
//                F3(ve, ho + 1) = F3(ve, ho) + dF3deAvg * dXi + SF3pre(ve, ho + 1);
//                F4(ve, ho + 1) = F4(ve, ho) + dF4deAvg * dXi + SF4pre(ve, ho + 1);
//            }
//            elseif(ve == 1)
//                % Predicted dFde[BOT]
//            dF1dePre(ve, ho + 1) = dEtadX(ve) * ((F1pre(ve, ho + 1) - F1pre(ve + 1, ho + 1)) / dEta) + 1 / h * ((G1pre(ve, ho + 1) - G1pre(ve + 1, ho + 1)) / dEta);
//            dF2dePre(ve, ho + 1) = dEtadX(ve) * ((F2pre(ve, ho + 1) - F2pre(ve + 1, ho + 1)) / dEta) + 1 / h * ((G2pre(ve, ho + 1) - G2pre(ve + 1, ho + 1)) / dEta);
//            dF3dePre(ve, ho + 1) = dEtadX(ve) * ((F3pre(ve, ho + 1) - F3pre(ve + 1, ho + 1)) / dEta) + 1 / h * ((G3pre(ve, ho + 1) - G3pre(ve + 1, ho + 1)) / dEta);
//            dF4dePre(ve, ho + 1) = dEtadX(ve) * ((F4pre(ve, ho + 1) - F4pre(ve + 1, ho + 1)) / dEta) + 1 / h * ((G4pre(ve, ho + 1) - G4pre(ve + 1, ho + 1)) / dEta);

//            dF1deAvg = (dF1de(ve, ho) + dF1dePre(ve, ho + 1)) / 2;
//            dF2deAvg = (dF2de(ve, ho) + dF2dePre(ve, ho + 1)) / 2;
//            dF3deAvg = (dF3de(ve, ho) + dF3dePre(ve, ho + 1)) / 2;
//            dF4deAvg = (dF4de(ve, ho) + dF4dePre(ve, ho + 1)) / 2;

//            F1(ve, ho + 1) = F1(ve, ho) + dF1deAvg * dXi;
//            F2(ve, ho + 1) = F2(ve, ho) + dF2deAvg * dXi;
//            F3(ve, ho + 1) = F3(ve, ho) + dF3deAvg * dXi;
//            F4(ve, ho + 1) = F4(ve, ho) + dF4deAvg * dXi;
//        else
//                dF1dePre(ve, ho + 1) = dEtadX(ve) * ((F1pre(ve - 1, ho + 1) - F1pre(ve, ho + 1)) / dEta) + 1 / h * ((G1pre(ve - 1, ho + 1) - G1pre(ve, ho + 1)) / dEta);
//            dF2dePre(ve, ho + 1) = dEtadX(ve) * ((F2pre(ve - 1, ho + 1) - F2pre(ve, ho + 1)) / dEta) + 1 / h * ((G2pre(ve - 1, ho + 1) - G2pre(ve, ho + 1)) / dEta);
//            dF3dePre(ve, ho + 1) = dEtadX(ve) * ((F3pre(ve - 1, ho + 1) - F3pre(ve, ho + 1)) / dEta) + 1 / h * ((G3pre(ve - 1, ho + 1) - G3pre(ve, ho + 1)) / dEta);
//            dF4dePre(ve, ho + 1) = dEtadX(ve) * ((F4pre(ve - 1, ho + 1) - F4pre(ve, ho + 1)) / dEta) + 1 / h * ((G4pre(ve - 1, ho + 1) - G4pre(ve, ho + 1)) / dEta);

//            dF1deAvg = (dF1de(ve, ho) + dF1dePre(ve, ho + 1)) / 2;
//            dF2deAvg = (dF2de(ve, ho) + dF2dePre(ve, ho + 1)) / 2;
//            dF3deAvg = (dF3de(ve, ho) + dF3dePre(ve, ho + 1)) / 2;
//            dF4deAvg = (dF4de(ve, ho) + dF4dePre(ve, ho + 1)) / 2;

//            F1(ve, ho + 1) = F1(ve, ho) + dF1deAvg * dXi;
//            F2(ve, ho + 1) = F2(ve, ho) + dF2deAvg * dXi;
//            F3(ve, ho + 1) = F3(ve, ho) + dF3deAvg * dXi;
//            F4(ve, ho + 1) = F4(ve, ho) + dF4deAvg * dXi;
//            end

//            % Magnitudes fisicas
//            A = F3(ve, ho + 1) ^ 2 / 2 / F1(ve, ho + 1) - F4(ve, ho + 1);
//            B = gamma / (gamma - 1) * F1(ve, ho + 1) * F2(ve, ho + 1);
//            C = -(gamma + 1) / 2 / (gamma - 1) * F1(ve, ho + 1) ^ 3;
//            ro(ve, ho + 1) = (-B + sqrt(B ^ 2 - 4 * A * C)) / 2 / A;

//            u(ve, ho + 1) = F1(ve, ho + 1) / ro(ve, ho + 1);
//            v(ve, ho + 1) = F3(ve, ho + 1) / F1(ve, ho + 1);
//            p(ve, ho + 1) = F2(ve, ho + 1) - F1(ve, ho + 1) * u(ve, ho + 1);
//            T(ve, ho + 1) = p(ve, ho + 1) / ro(ve, ho + 1) / R;
//            M(ve, ho + 1) = sqrt(v(ve, ho + 1) ^ 2 + u(ve, ho + 1) ^ 2) / sqrt(gamma * p(ve, ho + 1) / ro(ve, ho + 1));

//            if ve == 1
    

//            if xi > E
//                phi = theta - atan2(abs(v(ve, ho)), u(ve, ho));
//                else
//                    phi = atan2(v(ve, ho + 1), u(ve, ho + 1));
//            end

//            f_cal = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (M(ve, ho + 1) ^ 2 - 1))) - atan(sqrt(M(ve, ho + 1) ^ 2 - 1));

//            f_act = f_cal + phi;

//            var = 10;
//            tol = 1e-3;
//            guessed = 10e10;
//            precision = 0.0000001;
//            % PROVISIONA, CAMBIAR
//             a_int = 1.1;
//            b_int = 2.9;
//            zero_f1 = sqrt((gamma + 1) / (gamma - 1)) * (atan(sqrt(((gamma - 1) / (gamma + 1)) * (a_int ^ 2 - 1)))) - (atan(sqrt((a_int ^ 2) - 1))) - f_act; % Function used to find its zero
//                        zero_f2 = sqrt((gamma + 1) / (gamma - 1)) * (atan(sqrt(((gamma - 1) / (gamma + 1)) * (((a_int + b_int) / 2) ^ 2 - 1)))) - (atan(sqrt((((a_int + b_int) / 2) ^ 2) - 1))) - f_act;
//            while ((b_int - a_int) / 2 > precision)
//                if (zero_f1 * zero_f2 <= 0)
//                    b_int = (a_int + b_int) / 2;
//                else
//                    a_int = (a_int + b_int) / 2;
//            end
//            zero_f1 = sqrt((gamma + 1) / (gamma - 1)) * (atan(sqrt(((gamma - 1) / (gamma + 1)) * (a_int ^ 2 - 1)))) - (atan(sqrt((a_int ^ 2) - 1))) - f_act;
//            zero_f2 = sqrt((gamma + 1) / (gamma - 1)) * (atan(sqrt(((gamma - 1) / (gamma + 1)) * (((a_int + b_int) / 2) ^ 2 - 1)))) - (atan(sqrt((((a_int + b_int) / 2) ^ 2) - 1))) - f_act;
//            end
//            M_act = (a_int + b_int) / 2;

//            pAct(ve, ho + 1) = p(ve, ho + 1) * ((1 + ((gamma - 1) / 2) * M(ve, ho + 1) ^ 2) / (1 + ((gamma - 1) / 2) * M_act ^ 2)) ^ (gamma / (gamma - 1));
//            TAct(ve, ho + 1) = T(ve, ho + 1) * ((1 + ((gamma - 1) / 2) * M(ve, ho + 1) ^ 2) / (1 + ((gamma - 1) / 2) * M_act ^ 2));
//            roAct(ve, ho + 1) = pAct(ve, ho + 1) / R / TAct(ve, ho + 1);

//            p(ve, ho + 1) = pAct(ve, ho + 1);
//            T(ve, ho + 1) = TAct(ve, ho + 1);
//            ro(ve, ho + 1) = roAct(ve, ho + 1);
//            M(ve, ho + 1) = M_act;

//            if (xi > E)
//                v(ve, ho + 1) = -(u(ve, ho + 1) * tan(theta));
//            else
//                v(ve, ho + 1) = 0;
//            end


//            F1c(ve, ho + 1) = ro(ve, ho + 1) * u(ve, ho + 1);
//            F2c(ve, ho + 1) = ro(ve, ho + 1) * u(ve, ho + 1) ^ 2 + p(ve, ho + 1);
//            F3c(ve, ho + 1) = ro(ve, ho + 1) * u(ve, ho + 1) * v(ve, ho + 1);
//            F4c(ve, ho + 1) = gamma / (gamma - 1) * p(ve, ho + 1) * u(ve, ho + 1) + ro(ve, ho + 1) * u(ve, ho + 1) * (u(ve, ho + 1) ^ 2 + v(ve, ho + 1) ^ 2) / 2;


//% F1(ve, ho + 1) = F1c(ve, ho + 1);
//% F2(ve, ho + 1) = F2c(ve, ho + 1);
//% F3(ve, ho + 1) = F3c(ve, ho + 1);
//% F4(ve, ho + 1) = F4c(ve, ho + 1);

//            G1(ve, ho + 1) = ro(ve, ho + 1) * F3c(ve, ho + 1) / F1c(ve, ho + 1);
//            G2(ve, ho + 1) = F3c(ve, ho + 1);
//            G3(ve, ho + 1) = ro(ve, ho + 1) * (F3c(ve, ho + 1) / F1c(ve, ho + 1)) ^ 2 + F2c(ve, ho + 1) - F1c(ve, ho + 1) ^ 2 / ro(ve, ho + 1);
//            G4(ve, ho + 1) = gamma / (gamma - 1) * (F2c(ve, ho + 1) - F1c(ve, ho + 1) ^ 2 / ro(ve, ho + 1)) * F3c(ve, ho + 1) / F1c(ve, ho + 1) + ro(ve, ho + 1) / 2 * F3c(ve, ho + 1) / F1c(ve, ho + 1) * ((F1c(ve, ho + 1) / ro(ve, ho + 1)) ^ 2 + (F3c(ve, ho + 1) / F1c(ve, ho + 1)) ^ 2);


//        else
//        % New G
//        G1(ve, ho + 1) = ro(ve, ho + 1) * F3(ve, ho + 1) / F1(ve, ho + 1);
//            G2(ve, ho + 1) = F3(ve, ho + 1);
//            G3(ve, ho + 1) = ro(ve, ho + 1) * (F3(ve, ho + 1) / F1(ve, ho + 1)) ^ 2 + F2(ve, ho + 1) - F1(ve, ho + 1) ^ 2 / ro(ve, ho + 1);
//            G4(ve, ho + 1) = gamma / (gamma - 1) * (F2(ve, ho + 1) - F1(ve, ho + 1) ^ 2 / ro(ve, ho + 1)) * F3(ve, ho + 1) / F1(ve, ho + 1) + ro(ve, ho + 1) / 2 * F3(ve, ho + 1) / F1(ve, ho + 1) * ((F1(ve, ho + 1) / ro(ve, ho + 1)) ^ 2 + (F3(ve, ho + 1) / F1(ve, ho + 1)) ^ 2);
//            end
        //end
        }

        //public double getU()
        //{
        //    return u;
        //}
        //public double getV()
        //{
        //    return v;
        //}
        //public double getRO()
        //{
        //    return ro;
        //}
        //public double getP()
        //{
        //    return p;
        //}
        //public double getT()
        //{
        //    return T;
        //}
        //public double getM()
        //{
        //    return M;
        //}

        //public double getF(int i)
        //{
        //    if (i == 1) { return F1; }
        //    else if (i == 2) { return F2; }
        //    else if (i == 3) { return F3; }
        //    else { return F4; }
        //}

        //public double getG(int i)
        //{
        //    if (i == 1) { return G1; }
        //    else if (i == 2) { return G2; }
        //    else if (i == 3) { return G3; }
        //    else { return G4; }
        //}

        //public double getFpre(int i)
        //{
        //    if (i == 1) { return F1pre; }
        //    else if (i == 2) { return F2pre; }
        //    else if (i == 3) { return F3pre; }
        //    else { return F4pre; }
        //}

        //public double getGpre(int i)
        //{
        //    if (i == 1) { return G1pre; }
        //    else if (i == 2) { return G2pre; }
        //    else if (i == 3) { return G3pre; }
        //    else { return G4pre; }
        //}


    //}
}
