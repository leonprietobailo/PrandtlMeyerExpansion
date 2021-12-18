using Microsoft.Win32;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace PradntlMeyerExpansion
{
    public class Grid
    {
        List<Cell[]> Mesh = new List<Cell[]>();                             // Matriz almacenadora de los objetos "Cell".
        List<double> xP = new List<double>(0);                              // Vector con los valores de la posición horizontal.
        Rules r;                                                            // Reglas asociadas a la malla.
        double dEta, h, dy, xi, dXi;                                        // Valores caractreristicos de cada columna de la malla.
        int ho;                                                             // Posición horizontal de la malla.
        List<double[]> yP = new List<double[]>();                           // Matriz con los  valores de la posción vertical. 
        double[] dEtadX;                                                    // Diferencial de la posición horizontal.
        
        // Inicialización de la malla.
        public Grid(Rules rIn)
        {
            r = rIn;
            dEtadX = new double[r.getJ()];
            ho = 1;
            xi = 0;
            Cell[] row = new Cell[r.getJ()];
            double[] rowY = new double[r.getJ()];
            for (int i = 0; i < r.getJ(); i++) { row[i] = new Cell(r, i, ho); }
            for (int i = 0; i < r.getJ(); i++) { rowY[i] = i * r.getH() / (r.getJ() - 1); }
            Mesh.Add(row);
            xP.Add(0);
            yP.Add(rowY);
        }

        // Cálculo de la expansión de Prandtl-Meyer.
        public void PrandtlMeyerExpansion()
        {
            while (xi < r.getxMax())
            {
                Cell[] row = new Cell[r.getJ()];
                for (int i = 0; i < r.getJ(); i++) { row[i] = new Cell(r, i, ho); }
                Mesh.Add(row);
                InitVars();
                ComputeStepSize();

                // Cálculo del Predictor y Corrector Step.
                for (int i = 0; i < r.getJ(); i++) { Mesh[ho][i].PredictorStep(dEtadX[i], dEta, h, dXi, this); }
                for (int i = 0; i < r.getJ(); i++) { Mesh[ho][i].CorrectorStep(dEtadX[i], dEta, h, dXi, xi, this); }
                ho++;
                xi += dXi;
                xP.Add(xi);
            }
        }

        // Inicialización de las variables de la columna.
        private void InitVars()
        {
            double[] y = new double[r.getJ()];
            for (int i = 0; i < r.getJ(); i++)
            {
                double ys;
                if (xi < r.getE())
                {
                    ys = 0;
                    h = r.getH();
                }
                else
                {
                    ys = -(xi - r.getE()) * Math.Tan(r.getTheta());
                    h = r.getH() + (xi - r.getE()) * Math.Tan(r.getTheta());
                }
                dy = h / (r.getJ() - 1);
                y[i] = ys + dy * i;
                double eta = (y[i] - ys) / h;
                dEta = 1.0 / (r.getJ() - 1);

                if (xi < r.getE())
                {
                    dEtadX[i] = 0.0;
                }
                else
                {
                    dEtadX[i] = (1 - eta) * Math.Tan(r.getTheta()) / h;
                }
            }
            yP.Add(y);
        }

        // Cálculo del tamaño del step.
        private void ComputeStepSize()
        {
            double[] anglePlus = new double[r.getJ()];
            double[] angleMin = new double[r.getJ()];
            double mu;
            for (int i = 0; i < r.getJ(); i++)
            {
                mu = Math.Asin(1 / Mesh[ho - 1][i].getM());
                anglePlus[i] = Math.Abs(Math.Tan(r.getTheta() + mu));
                angleMin[i] = Math.Abs(Math.Tan(r.getTheta() - mu));
            }
            double[] angleMax = new double[2];
            angleMax[0] = anglePlus.Max();
            angleMax[1] = angleMin.Max();
            double MaximumValue = angleMax.Max();
            dXi = r.getCourant() * dy / MaximumValue;
        }

        // Obtener el valor medio de las magnitudes fisicas en el downstream.
        public (double u, double v, double ro, double p, double T, double M) getDownstream()
        {
            Cell ret = Mesh[Mesh.Count - 1][0];
            double u = 0;
            double v = 0;
            double RO = 0;
            double P = 0;
            double T = 0;
            double M = 0;
            int counter = 0;
            for(int i=1; i < Math.Floor(23.0/41.0 * r.getJ()); i++)
            {
                u += Mesh[Mesh.Count - 1][i].getU();
                v += Mesh[Mesh.Count - 1][i].getV();
                RO += Mesh[Mesh.Count - 1][i].getRO();
                P += Mesh[Mesh.Count - 1][i].getP();
                T += Mesh[Mesh.Count - 1][i].getT();
                M += Mesh[Mesh.Count - 1][i].getM();
                counter++;
            }
            return (u/counter, v / counter, RO / counter, P / counter, T / counter, M / counter);
        }

        // Guardar CSV de una magnitud en concreto.
        public void saveCSV(int magnitude)
        {
            SaveFileDialog diag = new SaveFileDialog();
            diag.Filter = "(*.csv)|*.*";
            diag.DefaultExt = "csv";
            if (diag.ShowDialog() == true)
            {
                for (int i = r.getJ() - 1; i >= 0; i--)
                {
                    for (int n = 0; n < Mesh.Count; n++)
                    {
                        if (magnitude == 0) { File.AppendAllText(diag.FileName, Mesh[n][i].getU().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }
                        if (magnitude == 1) { File.AppendAllText(diag.FileName, Mesh[n][i].getV().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }
                        if (magnitude == 2) { File.AppendAllText(diag.FileName, Mesh[n][i].getRO().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }
                        if (magnitude == 3) { File.AppendAllText(diag.FileName, Mesh[n][i].getP().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }
                        if (magnitude == 4) { File.AppendAllText(diag.FileName, Mesh[n][i].getT().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }
                        if (magnitude == 5) { File.AppendAllText(diag.FileName, Mesh[n][i].getM().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }

                        if (n < Mesh.Count - 1)
                        {
                            File.AppendAllText(diag.FileName, ",");
                        }
                        else
                        {
                            File.AppendAllText(diag.FileName, "\n");
                        }
                    }
                }
            }
        }

        // Getters de la clase.
        public Cell GetCell(int ve, int ho) { return Mesh[ho][ve]; }
        public List<double[]> GetYP() { return yP; }
        public List<double> GetXP() { return xP; }

    }
}

