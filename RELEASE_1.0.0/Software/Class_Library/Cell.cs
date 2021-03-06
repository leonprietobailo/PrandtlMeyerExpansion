using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Class_Library
{
    public class Cell
    {
        double u, v, ro, p, T, M;                                           // Magnitudes fisicas del fluido.
        int chara;                                                          // 0 -> Inferior (superficie), 1 -> Inflow, 2 -> Superior.
        int ho, ve;                                                         // Posicion de la celda en la malla.
        double F1, F2, F3, F4, G1, G2, G3, G4;                              // Valores de F y G.
        double F1pre, F2pre, F3pre, F4pre, G1pre, G2pre, G3pre, G4pre;      // Valores predichos de F y G.
        double roPre, pPre;                                                 // Valores predichos de densidad y presion.
        double dF1de, dF2de, dF3de, dF4de;                                  // Diferenciales de F.
        Rules r;                                                            // Reglas.


        // Inicialización de la celda.
        public Cell(Rules rIn, int veIn, int hoIn)
        {
            r = rIn;
            // Definir las condiciones iniciales de las celdas en la primera columna del grid.
            if (ho == 0)
            {
                u = r.getU();
                v = r.getV();
                ro = r.getRO();
                p = r.getP();
                T = r.getT();
                M = r.getM();

                F1 = ro * u;
                F2 = ro * Math.Pow(u, 2) + p;
                F3 = ro * u * v;
                F4 = r.getGamma() / (r.getGamma() - 1) * p * u + ro * u * (Math.Pow(u, 2) + Math.Pow(v, 2)) / 2.0;

                G1 = ro * F3 / F1;
                G2 = F3;
                G3 = ro * Math.Pow((F3 / F1), 2) + F2 - Math.Pow(F1, 2) / ro;
                G4 = r.getGamma() / (r.getGamma() - 1) * (F2 - Math.Pow(F1, 2) / ro) * F3 / F1 + ro / 2.0 * F3 / F1 * Math.Pow(F1 / ro, 2) + Math.Pow((F3 / F1), 2);

            }
            // Definir valores por defecto para el resto de celdas.
            else
            {
                u = 0;
                v = 0;
                ro = 0;
                p = 0;
                T = 0;
                M = 0;
            }
            // Asignar "Rol" a las celdas para distinguir si estan en la fila superior, inferior o en el inflow.
            if (veIn == 0)
            {
                chara = 0;
            }
            else if (veIn == (r.getJ() - 1))
            {
                chara = 2;
            }
            else
            {
                chara = 1;
            }
            ho = hoIn;
            ve = veIn;
        }

        // Cálculo del Predictor Step del método de MacCormack.
        public void PredictorStep(double dEtadX, double dEta, double h, double dXi, Grid g)
        {
            // Inflow.
            if (chara == 1)
            {
                // Cálculo de los diferenciales de F respecto xi.
                dF1de = dEtadX * ((g.GetCell(ve, ho - 1).getF(1) - g.GetCell(ve + 1, ho - 1).getF(1)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(1) - g.GetCell(ve + 1, ho - 1).getG(1)) / dEta);
                dF2de = dEtadX * ((g.GetCell(ve, ho - 1).getF(2) - g.GetCell(ve + 1, ho - 1).getF(2)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(2) - g.GetCell(ve + 1, ho - 1).getG(2)) / dEta);
                dF3de = dEtadX * ((g.GetCell(ve, ho - 1).getF(3) - g.GetCell(ve + 1, ho - 1).getF(3)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(3) - g.GetCell(ve + 1, ho - 1).getG(3)) / dEta);
                dF4de = dEtadX * ((g.GetCell(ve, ho - 1).getF(4) - g.GetCell(ve + 1, ho - 1).getF(4)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(4) - g.GetCell(ve + 1, ho - 1).getG(4)) / dEta);

                // Añadimos coeficiente de viscosidad para evitar oscilaciones.
                double SF1 = r.getCy() * Math.Abs(g.GetCell(ve + 1, ho - 1).getP() - 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) / (g.GetCell(ve + 1, ho - 1).getP() + 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) * (g.GetCell(ve - 1, ho - 1).getF(1) - 2 * g.GetCell(ve, ho - 1).getF(1) + g.GetCell(ve + 1, ho - 1).getF(1));
                double SF2 = r.getCy() * Math.Abs(g.GetCell(ve + 1, ho - 1).getP() - 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) / (g.GetCell(ve + 1, ho - 1).getP() + 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) * (g.GetCell(ve - 1, ho - 1).getF(2) - 2 * g.GetCell(ve, ho - 1).getF(2) + g.GetCell(ve + 1, ho - 1).getF(2));
                double SF3 = r.getCy() * Math.Abs(g.GetCell(ve + 1, ho - 1).getP() - 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) / (g.GetCell(ve + 1, ho - 1).getP() + 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) * (g.GetCell(ve - 1, ho - 1).getF(3) - 2 * g.GetCell(ve, ho - 1).getF(3) + g.GetCell(ve + 1, ho - 1).getF(3));
                double SF4 = r.getCy() * Math.Abs(g.GetCell(ve + 1, ho - 1).getP() - 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) / (g.GetCell(ve + 1, ho - 1).getP() + 2 * g.GetCell(ve, ho - 1).getP() + g.GetCell(ve - 1, ho - 1).getP()) * (g.GetCell(ve - 1, ho - 1).getF(4) - 2 * g.GetCell(ve, ho - 1).getF(4) + g.GetCell(ve + 1, ho - 1).getF(4));

                // Calculamos valores predichos de F.
                F1pre = g.GetCell(ve, ho - 1).getF(1) + dF1de * dXi + SF1;
                F2pre = g.GetCell(ve, ho - 1).getF(2) + dF2de * dXi + SF2;
                F3pre = g.GetCell(ve, ho - 1).getF(3) + dF3de * dXi + SF3;
                F4pre = g.GetCell(ve, ho - 1).getF(4) + dF4de * dXi + SF4;
            }
            // Inferior.
            else if (chara == 0)
            {
                // Cálculo de los diferenciales de F respecto xi.
                dF1de = dEtadX * ((g.GetCell(ve, ho - 1).getF(1) - g.GetCell(ve + 1, ho - 1).getF(1)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(1) - g.GetCell(ve + 1, ho - 1).getG(1)) / dEta);
                dF2de = dEtadX * ((g.GetCell(ve, ho - 1).getF(2) - g.GetCell(ve + 1, ho - 1).getF(2)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(2) - g.GetCell(ve + 1, ho - 1).getG(2)) / dEta);
                dF3de = dEtadX * ((g.GetCell(ve, ho - 1).getF(3) - g.GetCell(ve + 1, ho - 1).getF(3)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(3) - g.GetCell(ve + 1, ho - 1).getG(3)) / dEta);
                dF4de = dEtadX * ((g.GetCell(ve, ho - 1).getF(4) - g.GetCell(ve + 1, ho - 1).getF(4)) / dEta) + 1.0 / h * ((g.GetCell(ve, ho - 1).getG(4) - g.GetCell(ve + 1, ho - 1).getG(4)) / dEta);

                // Calculamos valores predichos de F.
                F1pre = g.GetCell(ve, ho - 1).getF(1) + dF1de * dXi;
                F2pre = g.GetCell(ve, ho - 1).getF(2) + dF2de * dXi;
                F3pre = g.GetCell(ve, ho - 1).getF(3) + dF3de * dXi;
                F4pre = g.GetCell(ve, ho - 1).getF(4) + dF4de * dXi;
            }
            // Superior.
            else
            {
                // Cálculo de los diferenciales de F respecto xi.
                dF1de = dEtadX * ((g.GetCell(ve - 1, ho - 1).getF(1) - g.GetCell(ve, ho - 1).getF(1)) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho - 1).getG(1) - g.GetCell(ve, ho - 1).getG(1)) / dEta);
                dF2de = dEtadX * ((g.GetCell(ve - 1, ho - 1).getF(2) - g.GetCell(ve, ho - 1).getF(2)) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho - 1).getG(2) - g.GetCell(ve, ho - 1).getG(2)) / dEta);
                dF3de = dEtadX * ((g.GetCell(ve - 1, ho - 1).getF(3) - g.GetCell(ve, ho - 1).getF(3)) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho - 1).getG(3) - g.GetCell(ve, ho - 1).getG(3)) / dEta);
                dF4de = dEtadX * ((g.GetCell(ve - 1, ho - 1).getF(4) - g.GetCell(ve, ho - 1).getF(4)) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho - 1).getG(4) - g.GetCell(ve, ho - 1).getG(4)) / dEta);

                // Calculamos valores predichos de F.
                F1pre = g.GetCell(ve, ho - 1).getF(1) + dF1de * dXi;
                F2pre = g.GetCell(ve, ho - 1).getF(2) + dF2de * dXi;
                F3pre = g.GetCell(ve, ho - 1).getF(3) + dF3de * dXi;
                F4pre = g.GetCell(ve, ho - 1).getF(4) + dF4de * dXi;
            }

            // Cálculo de las magnitudes predichas, para posteriormente utilizarlas en el corrector step.
            double A = Math.Pow(F3pre, 2) / 2.0 / F1pre - F4pre;
            double B = r.getGamma() / (r.getGamma() - 1) * F1pre * F2pre;
            double C = -(r.getGamma() + 1) / 2.0 / (r.getGamma() - 1) * Math.Pow(F1pre, 3);
            roPre = (-B + Math.Sqrt(Math.Pow(B, 2) - 4 * A * C)) / 2.0 / A;
            G1pre = roPre * F3pre / F1pre;
            G2pre = F3pre;
            G3pre = roPre * Math.Pow(F3pre / F1pre, 2) + F2pre - Math.Pow(F1pre, 2) / roPre;
            G4pre = r.getGamma() / (r.getGamma() - 1) * (F2pre - Math.Pow(F1pre, 2) / roPre) * F3pre / F1pre + roPre / 2 * F3pre / F1pre * (Math.Pow(F1pre / roPre, 2) + Math.Pow(F3pre / F1pre, 2));
            pPre = F2pre - Math.Pow(F1pre, 2) / roPre;
        }


        // Cálculo del Predictor Step del método de MacCormack.
        public void CorrectorStep(double dEtadX, double dEta, double h, double dXi, double Xi, Grid g)
        {
            // Inflow.
            if (chara == 1)
            {
                // Cálculo de los diferenciales de F respecto xi medio predicho.
                double dF1dePre = dEtadX * ((g.GetCell(ve - 1, ho).getFpre(1) - F1pre) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho).getGpre(1) - G1pre) / dEta);
                double dF2dePre = dEtadX * ((g.GetCell(ve - 1, ho).getFpre(2) - F2pre) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho).getGpre(2) - G2pre) / dEta);
                double dF3dePre = dEtadX * ((g.GetCell(ve - 1, ho).getFpre(3) - F3pre) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho).getGpre(3) - G3pre) / dEta);
                double dF4dePre = dEtadX * ((g.GetCell(ve - 1, ho).getFpre(4) - F4pre) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho).getGpre(4) - G4pre) / dEta);

                // Promedio de los diferenciales.
                double dF1deAvg = (dF1de + dF1dePre) / 2.0;
                double dF2deAvg = (dF2de + dF2dePre) / 2.0;
                double dF3deAvg = (dF3de + dF3dePre) / 2.0;
                double dF4deAvg = (dF4de + dF4dePre) / 2.0;

                // Añadimos coeficiente de viscosidad para evitar oscilaciones.
                double SF1pre = r.getCy() * Math.Abs(g.GetCell(ve + 1, ho).getPpre() - 2.0 * pPre + g.GetCell(ve - 1, ho - 1).getPpre()) / (g.GetCell(ve + 1, ho).getPpre() + 2 * pPre + g.GetCell(ve - 1, ho).getPpre()) * (g.GetCell(ve - 1, ho).getFpre(1) - 2 * F1pre + g.GetCell(ve + 1, ho).getFpre(1));
                double SF2pre = r.getCy() * Math.Abs(g.GetCell(ve + 1, ho).getPpre() - 2.0 * pPre + g.GetCell(ve - 1, ho - 1).getPpre()) / (g.GetCell(ve + 1, ho).getPpre() + 2 * pPre + g.GetCell(ve - 1, ho).getPpre()) * (g.GetCell(ve - 1, ho).getFpre(2) - 2 * F2pre + g.GetCell(ve + 1, ho).getFpre(2));
                double SF3pre = r.getCy() * Math.Abs(g.GetCell(ve + 1, ho).getPpre() - 2.0 * pPre + g.GetCell(ve - 1, ho - 1).getPpre()) / (g.GetCell(ve + 1, ho).getPpre() + 2 * pPre + g.GetCell(ve - 1, ho).getPpre()) * (g.GetCell(ve - 1, ho).getFpre(3) - 2 * F3pre + g.GetCell(ve + 1, ho).getFpre(3));
                double SF4pre = r.getCy() * Math.Abs(g.GetCell(ve + 1, ho).getPpre() - 2.0 * pPre + g.GetCell(ve - 1, ho - 1).getPpre()) / (g.GetCell(ve + 1, ho).getPpre() + 2 * pPre + g.GetCell(ve - 1, ho).getPpre()) * (g.GetCell(ve - 1, ho).getFpre(4) - 2 * F4pre + g.GetCell(ve + 1, ho).getFpre(4));

                // Cálculo de los valores de F para las celdas.
                F1 = g.GetCell(ve, ho - 1).getF(1) + dF1deAvg * dXi + SF1pre;
                F2 = g.GetCell(ve, ho - 1).getF(2) + dF2deAvg * dXi + SF2pre;
                F3 = g.GetCell(ve, ho - 1).getF(3) + dF3deAvg * dXi + SF3pre;
                F4 = g.GetCell(ve, ho - 1).getF(4) + dF4deAvg * dXi + SF4pre;
            }
            // Inferior.
            else if (chara == 0)
            {
                double dF1dePre = dEtadX * ((F1pre - g.GetCell(ve + 1, ho).getFpre(1)) / dEta) + 1.0 / h * ((G1pre - g.GetCell(ve + 1, ho).getGpre(1)) / dEta);
                double dF2dePre = dEtadX * ((F2pre - g.GetCell(ve + 1, ho).getFpre(2)) / dEta) + 1.0 / h * ((G2pre - g.GetCell(ve + 1, ho).getGpre(2)) / dEta);
                double dF3dePre = dEtadX * ((F3pre - g.GetCell(ve + 1, ho).getFpre(3)) / dEta) + 1.0 / h * ((G3pre - g.GetCell(ve + 1, ho).getGpre(3)) / dEta);
                double dF4dePre = dEtadX * ((F4pre - g.GetCell(ve + 1, ho).getFpre(4)) / dEta) + 1.0 / h * ((G4pre - g.GetCell(ve + 1, ho).getGpre(4)) / dEta);

                // Promedio de los diferenciales.
                double dF1deAvg = (dF1de + dF1dePre) / 2.0;
                double dF2deAvg = (dF2de + dF2dePre) / 2.0;
                double dF3deAvg = (dF3de + dF3dePre) / 2.0;
                double dF4deAvg = (dF4de + dF4dePre) / 2.0;

                // Cálculo de los valores de F para las celdas.
                F1 = g.GetCell(ve, ho - 1).getF(1) + dF1deAvg * dXi;
                F2 = g.GetCell(ve, ho - 1).getF(2) + dF2deAvg * dXi;
                F3 = g.GetCell(ve, ho - 1).getF(3) + dF3deAvg * dXi;
                F4 = g.GetCell(ve, ho - 1).getF(4) + dF4deAvg * dXi;
            }
            // Superior.
            else
            {
                // Cálculo de los diferenciales de F respecto xi medio predicho.
                double dF1dePre = dEtadX * ((g.GetCell(ve - 1, ho).getFpre(1) - F1pre) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho).getGpre(1) - G1pre) / dEta);
                double dF2dePre = dEtadX * ((g.GetCell(ve - 1, ho).getFpre(2) - F2pre) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho).getGpre(2) - G2pre) / dEta);
                double dF3dePre = dEtadX * ((g.GetCell(ve - 1, ho).getFpre(3) - F3pre) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho).getGpre(3) - G3pre) / dEta);
                double dF4dePre = dEtadX * ((g.GetCell(ve - 1, ho).getFpre(4) - F4pre) / dEta) + 1.0 / h * ((g.GetCell(ve - 1, ho).getGpre(4) - G4pre) / dEta);

                // Promedio de los diferenciales.
                double dF1deAvg = (dF1de + dF1dePre) / 2.0;
                double dF2deAvg = (dF2de + dF2dePre) / 2.0;
                double dF3deAvg = (dF3de + dF3dePre) / 2.0;
                double dF4deAvg = (dF4de + dF4dePre) / 2.0;

                // Cálculo de los valores de F para las celdas.
                F1 = g.GetCell(ve, ho - 1).getF(1) + dF1deAvg * dXi;
                F2 = g.GetCell(ve, ho - 1).getF(2) + dF2deAvg * dXi;
                F3 = g.GetCell(ve, ho - 1).getF(3) + dF3deAvg * dXi;
                F4 = g.GetCell(ve, ho - 1).getF(4) + dF4deAvg * dXi;
            }

            // Revertimos el cambio de F para hallar las magnitudes fisicas.
            double A = Math.Pow(F3, 2) / 2.0 / F1 - F4;
            double B = r.getGamma() / (r.getGamma() - 1) * F1 * F2;
            double C = -(r.getGamma() + 1) / 2 / (r.getGamma() - 1) * Math.Pow(F1, 3);
            ro = (-B + Math.Sqrt(Math.Pow(B, 2) - 4 * A * C)) / 2.0 / A;
            u = F1 / ro;
            v = F3 / F1;
            p = F2 - F1 * u;
            T = p / ro / r.getR();
            M = Math.Sqrt(Math.Pow(v, 2) + Math.Pow(u, 2)) / Math.Sqrt(r.getGamma() * p / ro);

            // Cálculo del número mach para la frontera inferior. 
            if (chara == 0)
            {
                // Calculamos el ángulo de los vectores de velocidad.
                double phi;
                if (Xi > r.getE())
                {
                    phi = r.getTheta() - Math.Atan2(Math.Abs(v), u);
                }
                else
                {
                    phi = Math.Atan2(v, u);
                }

                // Calculamos el valor de la funcion de Prandtl-Meyer.
                double f_cal = Math.Sqrt((r.getGamma() + 1) / (r.getGamma() - 1)) * Math.Atan(Math.Sqrt((r.getGamma() - 1) / (r.getGamma() + 1) * (Math.Pow(M, 2) - 1))) - Math.Atan(Math.Sqrt(Math.Pow(M, 2) - 1));
                // Hallamos el valor real añadiendole el ángulo.
                double f_act = f_cal + phi;
                // Buscamos el valor de M correspondiente al valor real que anula la funcion de Prandtl-Meyer.
                double precision = 0.0000001;
                double a_int = 1;
                double b_int = 10;
                double zero_f1 = Math.Sqrt((r.getGamma() + 1) / (r.getGamma() - 1)) * (Math.Atan(Math.Sqrt(((r.getGamma() - 1) / (r.getGamma() + 1)) * (Math.Pow(a_int, 2) - 1)))) - (Math.Atan(Math.Sqrt((Math.Pow(a_int, 2) - 1)))) - f_act;
                double zero_f2 = Math.Sqrt((r.getGamma() + 1) / (r.getGamma() - 1)) * (Math.Atan(Math.Sqrt(((r.getGamma() - 1) / (r.getGamma() + 1)) * (Math.Pow((a_int + b_int) / 2, 2) - 1)))) - (Math.Atan(Math.Sqrt((((Math.Pow((a_int + b_int) / 2, 2) - 1)))))) - f_act;
                // Guessing del numero mach.
                while ((b_int - a_int) / 2 > precision)
                {
                    if (zero_f1 * zero_f2 <= 0)
                    {
                        b_int = (a_int + b_int) / 2.0;
                    }
                    else
                    {
                        a_int = (a_int + b_int) / 2.0;
                    }

                    zero_f1 = Math.Sqrt((r.getGamma() + 1) / (r.getGamma() - 1)) * (Math.Atan(Math.Sqrt(((r.getGamma() - 1) / (r.getGamma() + 1)) * (Math.Pow(a_int, 2) - 1)))) - (Math.Atan(Math.Sqrt((Math.Pow(a_int, 2) - 1)))) - f_act;
                    zero_f2 = Math.Sqrt((r.getGamma() + 1) / (r.getGamma() - 1)) * (Math.Atan(Math.Sqrt(((r.getGamma() - 1) / (r.getGamma() + 1)) * (Math.Pow((a_int + b_int) / 2, 2) - 1)))) - (Math.Atan(Math.Sqrt((((Math.Pow((a_int + b_int) / 2, 2) - 1)))))) - f_act;
                }

                // Cálculo de las magnitudes reales.
                double M_act = (a_int + b_int) / 2.0;
                double pAct = p * Math.Pow((1 + ((r.getGamma() - 1) / 2.0) * Math.Pow(M_act, 2)) / (1 + ((r.getGamma() - 1) / 2.0) * Math.Pow(M_act, 2)), (r.getGamma() / (r.getGamma() - 1)));
                double TAct = T * ((1 + ((r.getGamma() - 1) / 2.0) * Math.Pow(M_act, 2)) / (1 + ((r.getGamma() - 1) / 2) * Math.Pow(M_act, 2)));
                double roAct = pAct / r.getR() / TAct;
                p = pAct;
                T = TAct;
                ro = roAct;
                M = M_act;

                // Cálculo de la velocidad vertical.
                if (Xi > r.getE())
                {
                    v = -(u * Math.Tan(r.getTheta()));
                }
                else
                {
                    v = 0;
                }

                // Cálculo de los valores de F para recuperar las magnitudes fiscias.
                F1 = ro * u;
                F2 = ro * Math.Pow(u, 2) + p;
                F3 = ro * u * v;
                F4 = r.getGamma() / (r.getGamma() - 1) * p * u + ro * u * (Math.Pow(u, 2) + Math.Pow(v, 2)) / 2.0;
            }

            // Cálculo de G.
            G1 = ro * F3 / F1;
            G2 = F3;
            G3 = ro * Math.Pow((F3 / F1), 2) + F2 - Math.Pow(F1, 2) / ro;
            G4 = r.getGamma() / (r.getGamma() - 1) * (F2 - Math.Pow(F1, 2) / ro) * F3 / F1 + ro / 2.0 * F3 / F1 * Math.Pow(F1 / ro, 2) + Math.Pow((F3 / F1), 2);
        }

        // Recuperar velocidad horizontal.
        public double getU() { return u; }
        // Recuperar velocidad vertical.
        public double getV() { return v; }
        // Recuperar densidad.
        public double getRO() { return ro; }
        // Recuperar presion.
        public double getP() { return p; }
        // Recuprerar temperatura.
        public double getT() { return T; }
        // Recuperar numero mach.
        public double getM() { return M; }
        // Recuperar valores de F.
        public double getF(int i)
        {
            if (i == 1) { return F1; }
            else if (i == 2) { return F2; }
            else if (i == 3) { return F3; }
            else { return F4; }
        }
        // Recuperar valores de G.
        public double getG(int i)
        {
            if (i == 1) { return G1; }
            else if (i == 2) { return G2; }
            else if (i == 3) { return G3; }
            else { return G4; }
        }
        // Recuperar valores de F predicho.
        public double getFpre(int i)
        {
            if (i == 1) { return F1pre; }
            else if (i == 2) { return F2pre; }
            else if (i == 3) { return F3pre; }
            else { return F4pre; }
        }
        // Recuperar valores de G predicho.
        public double getGpre(int i)
        {
            if (i == 1) { return G1pre; }
            else if (i == 2) { return G2pre; }
            else if (i == 3) { return G3pre; }
            else { return G4pre; }
        }
        // Recuperar presión predicha.
        public double getPpre() { return pPre; }
    }
}
