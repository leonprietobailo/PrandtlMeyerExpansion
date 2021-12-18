using Microsoft.Win32;
using System;
using System.IO;

namespace PradntlMeyerExpansion
{
    public class Rules
    {
        double u0, v0, ro0, p0, T0, M0;                                     // Magnitudes iniciales para el flujo.
        double Cy0, gamma0, R0, E0, theta0, xMax0, H0, Courant0;            // Condiciones de la malla y de estabilidad.
        int j0;                                                             // Numero de celdas verticales.

        // Inicialización de reglas, con atributos null.
        public Rules() { }
        // Inicializacion de reglas, con valores heredados de otro objeto reglas y con un angulo Theta diferente.
        public Rules(Rules r, double thetaIn)
        {
            u0 = r.u0;
            v0 = r.v0;
            ro0 = r.ro0;
            p0 = r.p0;
            T0 = r.T0;
            Cy0 = r.Cy0;
            gamma0 = r.gamma0;
            R0 = r.R0;
            E0 = r.E0;
            theta0 = thetaIn;
            xMax0 = r.xMax0;
            H0 = r.H0;
            Courant0 = r.Courant0;
            j0 = r.j0;
            M0 = Math.Sqrt(u0*u0 + v0*v0) / Math.Sqrt(gamma0 * R0 * T0);
        }
        // Inicializacion de reglas, con atributos inicializados a partir del constructor.
        public Rules(double u, double v, double ro, double p, double T, double Cy, double gamma, double R, double E, double theta, int j, double xMax, double H, double Courant)
        {
            u0 = u;
            v0 = v;
            ro0 = ro;
            p0 = p;
            T0 = T;
            Cy0 = Cy;
            gamma0 = gamma;
            R0 = R;
            E0 = E;
            theta0 = theta;
            j0 = j;
            xMax0 = xMax;
            H0 = H;
            Courant0 = Courant;
            M0 = Math.Sqrt(u0 * u0 + v0 * v0) / Math.Sqrt(gamma0 * R0 * T0);
        }

        // Almacenar los valores del objeto en un archivo en formato de texto.
        public void saveRules()
        {
            SaveFileDialog diag = new SaveFileDialog();
            diag.Filter = "(*.sim)|*.*";
            diag.DefaultExt = "sim";

            if (diag.ShowDialog() == true)
            {
                FileStream emptyFile = File.Create(diag.FileName);
                emptyFile.Close();
                File.AppendAllText(diag.FileName, Convert.ToString(u0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(v0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(ro0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(p0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(T0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(M0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(Cy0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(gamma0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(R0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(E0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(theta0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(j0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(xMax0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(H0));
                File.AppendAllText(diag.FileName, "-");

                File.AppendAllText(diag.FileName, Convert.ToString(Courant0));
            }
        }

        // Inicializar la clase reglas con un archivo en formato de texto.
        public int loadRules()
        {
            try
            {
                OpenFileDialog dig = new OpenFileDialog();
                dig.Multiselect = false;
                dig.Filter = "(*.txt*)|*.*";
                dig.DefaultExt = ".txt";
                if (dig.ShowDialog() == true)
                {
                    StreamReader readFile = new StreamReader(dig.FileName);
                    string strReadline = readFile.ReadLine();
                    string[] readed = strReadline.Split('-');

                    u0 = Convert.ToDouble(readed[0]);
                    v0 = Convert.ToDouble(readed[1]);
                    ro0 = Convert.ToDouble(readed[2]);
                    p0 = Convert.ToDouble(readed[3]);
                    T0 = Convert.ToDouble(readed[4]);
                    M0 = Convert.ToDouble(readed[5]);
                    Cy0 = Convert.ToDouble(readed[6]);
                    gamma0 = Convert.ToDouble(readed[7]);
                    R0 = Convert.ToDouble(readed[8]);
                    E0 = Convert.ToDouble(readed[9]);
                    theta0 = Convert.ToDouble(readed[10]);
                    j0 = Convert.ToInt32(readed[11]);
                    xMax0 = Convert.ToDouble(readed[12]);
                    H0 = Convert.ToDouble(readed[13]);
                    Courant0 = Convert.ToDouble(readed[14]);
                    return 0;
                }
                else
                {
                    return -1;
                }
            }
            //catch (FileFormatException)
            //{
            //    return -1;
            //}
            //catch (FormatException)
            //{
            //    return -1;
            //}
            catch
            {
                return -1;
            }
        }

        // Getters de la clase reglas.
        public double getU() { return u0; }
        public double getV() { return v0; }
        public double getRO() { return ro0; }
        public double getP() { return p0; }
        public double getT() { return T0; }
        public double getM() { return M0; }
        public double getCy() { return Cy0; }
        public double getGamma() { return gamma0; }
        public double getR() { return R0; }
        public double getH() { return H0; }
        public double getxMax() { return xMax0; }
        public int getJ() { return j0; }
        public double getCourant() { return Courant0; }
        public double getE() { return E0; }
        public double getTheta() { return theta0; }
    }
}
