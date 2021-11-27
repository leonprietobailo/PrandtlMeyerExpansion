﻿using System;
using System.Collections.Generic;
using System.Linq;

namespace PradntlMeyerExpansion
{
    public class Grid
    {
        List<Cell[]> Mesh = new List<Cell[]>();
        Rules r;
        double dEta, h, dy, xi, dXi;
        int ho;
        List<double[]> yP = new List<double[]>();
        double[] dEtadX;

        public Grid(Rules rIn)
        {
            r = rIn;
            dEtadX = new double[r.getJ()];
            ho = 0;
            xi = 0;
        }

        public void PrandtlMeyerExpansion()
        {
            while (xi < r.getxMax())
            {
                InitVars();
                ComputeStepSize();
                for (int i = 1; i < r.getJ(); i++) { Mesh[ho][i].PredictorStep(dEtadX[i], dEta, h, dXi, this); } //PREDICTOR STEP
                for (int i = 1; i < r.getJ(); i++) { Mesh[ho][i].CorrectorStep(dEtadX[i], dEta, h, dXi, xi, this); } //CORRECTOR STEP
                ho++;
                xi += dXi;
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
                y[i] = ys + dy * (i - 1);
                double eta = (y[i] - ys) / h;
                dEta = 1 / (r.getJ() - 1);

                if (xi < r.getE())
                {
                    dEtadX[i] = 0;
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
                mu = Math.Asin(1 / Mesh[ho][i].getM());
                anglePlus[i] = Math.Abs(Math.Tan(r.getTheta() + mu));
                angleMin[i] = Math.Abs(Math.Tan(r.getTheta() - mu));
            }
            double[] angleMax = new double[2];
            angleMax[0] = anglePlus.Max();
            angleMax[1] = angleMin.Max();
            double MaximumValue = angleMax.Max();
            dXi = r.getCourant() * dy / MaximumValue;
        }


        public Cell GetCell(int ve, int ho)
        {
            return Mesh[ho][ve];
        }


    }
}

