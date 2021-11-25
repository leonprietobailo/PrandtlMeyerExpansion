using System;

namespace PradntlMeyerExpansion
{
    public class Rules
    {
        double u0, v0, ro0, p0, T0, M0; // Initial Magnitudes
        double Cy0, gamma0;

        public Rules(double u, double v, double ro, double p, double T, double M, double Cy, double gamma)
        {
            u0 = u;
            v0 = v;
            ro0 = ro;
            p0 = p;
            T0 = T;
            M0 = M;
            Cy0 = Cy;
            gamma0 = gamma;
        }

        public double getU()
        {
            return u0;
        }
        public double getV()
        {
            return v0;
        }
        public double getRO()
        {
            return ro0;
        }
        public double getP()
        {
            return p0;
        }
        public double getT()
        {
            return T0;
        }
        public double getM()
        {
            return M0;
        }
        public double getCy()
        {
            return Cy0;
        }
        public double getGamma()
        {
            return gamma0;
        }
    }
}
