using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace PradntlMeyerExpansion
{
    public class Grid
    {
        List<Cell[]> Mesh = new List<Cell[]>();
        List<double> xP = new List<double>(0);
        Rules r;
        double dEta, h, dy, xi, dXi;
        int ho;
        List<double[]> yP = new List<double[]>();
        double[] dEtadX;

        public Grid(Rules rIn)
        {
            r = rIn;
            dEtadX = new double[r.getJ()];
            ho = 1;
            xi = 0;
            Cell[] row = new Cell[r.getJ()];
            double[] rowY = new double[r.getJ()];
            for (int i = 0; i < r.getJ(); i++) { row[i] = new Cell(r, i, ho); }
            for (int i = 0; i < r.getJ(); i++) { rowY[i] = i * r.getH()/(r.getJ()-1); }
            Mesh.Add(row);
            xP.Add(0);
            yP.Add(rowY);
        }
        

        public void PrandtlMeyerExpansion()
        {
            while (xi < r.getxMax())
            {
                Cell[] row = new Cell[r.getJ()];
                for (int i = 0; i < r.getJ(); i++) { row[i] = new Cell(r,i,ho); }
                Mesh.Add(row);
                InitVars();
                ComputeStepSize();
                for (int i = 0; i < r.getJ(); i++) { Mesh[ho][i].PredictorStep(dEtadX[i], dEta, h, dXi, this); } //PREDICTOR STEP
                for (int i = 0; i < r.getJ(); i++) { Mesh[ho][i].CorrectorStep(dEtadX[i], dEta, h, dXi, xi, this); } //CORRECTOR STEP
                ho++;
                xi += dXi;
                xP.Add(xi);
            }
        }

        private void InitVars()
        {
            double[] y = new double[r.getJ()];
            for(int i = 0; i < r.getJ(); i++)
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

        private void ComputeStepSize()
        {
            double[] anglePlus = new double[r.getJ()];
            double[] angleMin = new double[r.getJ()];
            double mu;
            for(int i = 0; i < r.getJ(); i++)
            {
                mu = Math.Asin(1 / Mesh[ho-1][i].getM());
                anglePlus[i] = Math.Abs(Math.Tan(r.getTheta() + mu));
                angleMin[i] = Math.Abs(Math.Tan(r.getTheta() - mu));
            }
            double[] angleMax = new double[2];
            angleMax[0] = anglePlus.Max();
            angleMax[1] = angleMin.Max();
            double MaximumValue = angleMax.Max();
            dXi = r.getCourant() * dy / MaximumValue;
        }

        public (double u, double v, double ro, double p , double T, double M) getDownstream()
        {
            Cell ret = Mesh[Mesh.Count - 1][0];
            return (ret.getU(), ret.getV(), ret.getRO(), ret.getP(), ret.getT(), ret.getM());
        }


        public Cell GetCell(int ve, int ho) { return Mesh[ho][ve]; }

        //public int GetHorizontalPoints() { return ho; }
        public List<double[]> GetYP() { return yP; }
        public List<double> GetXP() { return xP; }

        public void data()
        {
            FileStream emptyFile = File.Create("test.txt");
            emptyFile.Close();
            for (int i=0;i<Mesh[0].Length;i++)
            {
                File.AppendAllText("text.txt", Convert.ToString(Mesh[18][i].getM()));
                File.AppendAllText("text.txt", "\n");
            }
            File.AppendAllText("text.txt", "\n");
            for (int i = 0; i < Mesh[0].Length; i++)
            {
                File.AppendAllText("text.txt", Convert.ToString(Mesh[Mesh.Count-1][i].getM()));
                File.AppendAllText("text.txt", "\n");
            }
        }

    }
}

