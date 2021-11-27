using System;

namespace PradntlMeyerExpansion
{
    public class Rules
    {
        double u0, v0, ro0, p0, T0, M0; // Initial Magnitudes
        double Cy0, gamma0, R0, E0, theta0, xMax0, H0, Courant0;
        int j0;

        public Rules(double u, double v, double ro, double p, double T, double M, double Cy, double gamma, double R, double E, double theta, int j, double xMax, double H, double Courant)
        {
            u0 = u;
            v0 = v;
            ro0 = ro;
            p0 = p;
            T0 = T;
            M0 = M;
            Cy0 = Cy;
            gamma0 = gamma;
            R0 = R;
            E0 = E;
            theta0 = theta;
            j0 = j;
            xMax0 = xMax;
            H0 = H;
            Courant0 = Courant;

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
        public double getR()
        {
            return R0;
        }

        public double getH() { return H0; }
        public double getxMax() { return xMax0; }
        public int getJ() { return j0; }
        public double getCourant() { return Courant0; }
        public double getE() { return E0; }
        public double getTheta() { return theta0; }


    }
}
