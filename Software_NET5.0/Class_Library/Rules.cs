using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Class_Library
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
            M0 = Math.Sqrt(u0 * u0 + v0 * v0) / Math.Sqrt(gamma0 * R0 * T0);
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
            // Cálculo del numero Mach.
            M0 = Math.Sqrt(u0 * u0 + v0 * v0) / Math.Sqrt(gamma0 * R0 * T0);
        }

        // Almacenar los valores del objeto en un archivo en formato de texto.
        public void saveRules(string filename)
        {

            // Creamos y cerramos un nuervo archivo con el titulo espeficicado por el usuario.
            FileStream emptyFile = File.Create(filename);
            emptyFile.Close();
            // Añadimos velocidad horizontal y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(u0));
            File.AppendAllText(filename, "-");
            // Añadimos velocidad vertical y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(v0));
            File.AppendAllText(filename, "-");
            // Añadimos denisodad y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(ro0));
            File.AppendAllText(filename, "-");
            // Añadimos presion y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(p0));
            File.AppendAllText(filename, "-");
            // Añadimos temperatura y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(T0));
            File.AppendAllText(filename, "-");
            // Añadimos numero Mach y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(M0));
            File.AppendAllText(filename, "-");
            // Añadimos condicion de courant para la viscosidad y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(Cy0));
            File.AppendAllText(filename, "-");
            // Añadimos coeficiente adiabatico y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(gamma0));
            File.AppendAllText(filename, "-");
            // Añadimos constante de los gases ideales y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(R0));
            File.AppendAllText(filename, "-");
            // Añadimos ubicación de la esquina de expansión y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(E0));
            File.AppendAllText(filename, "-");
            // Añadimos ángulo y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(theta0));
            File.AppendAllText(filename, "-");
            // Añadimos divisiones verticales y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(j0));
            File.AppendAllText(filename, "-");
            // Añadimos longitud màxima y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(xMax0));
            File.AppendAllText(filename, "-");
            // Añadimos altura inicial y separamos con guión.
            File.AppendAllText(filename, Convert.ToString(H0));
            File.AppendAllText(filename, "-");
            // Añadimos condición de courant.
            File.AppendAllText(filename, Convert.ToString(Courant0));
        }

        // Inicializar la clase reglas con un archivo en formato de texto.
        public int loadRules(string filename)
        {
            // Leemos el archivo.
            StreamReader readFile = new StreamReader(filename);
            // Seleccionamos la primera linea.
            string strReadline = readFile.ReadLine();
            // Separamos por guiones.
            string[] readed = strReadline.Split('-');
            // Almacenamos velocidad horizontal en el objeto rules.
            u0 = Convert.ToDouble(readed[0]);
            // Almacenamos velocidad vertical en el objeto rules.
            v0 = Convert.ToDouble(readed[1]);
            // Almacenamos densidad en el objeto rules.
            ro0 = Convert.ToDouble(readed[2]);
            // Almacenamos presión en el objeto rules.
            p0 = Convert.ToDouble(readed[3]);
            // Almacenamos temperatura en el objeto rules.
            T0 = Convert.ToDouble(readed[4]);
            // Almacenamos numero Mach en el objeto rules.
            M0 = Convert.ToDouble(readed[5]);
            // Almacenamos condicion de courant viscosa en el objeto rules.
            Cy0 = Convert.ToDouble(readed[6]);
            // Almacenamos coeficiente adiabatico en el objeto rules.
            gamma0 = Convert.ToDouble(readed[7]);
            // Almacenamos coeficiente de los gases ideales en el objeto rules.
            R0 = Convert.ToDouble(readed[8]);
            // Almacenamos ubicacion de la esquina de expansión en el objeto rules.
            E0 = Convert.ToDouble(readed[9]);
            // Almacenamos ángulo en el objeto rules.
            theta0 = Convert.ToDouble(readed[10]);
            // Almacenamos numero de divisiones verticales en el objeto rules.
            j0 = Convert.ToInt32(readed[11]);
            // Almacenamos distancia máxima en el objeto rules.
            xMax0 = Convert.ToDouble(readed[12]);
            // Almacenamos altura inicial en el objeto rules.
            H0 = Convert.ToDouble(readed[13]);
            // Almacenamos condición de courant en el objeto rules.
            Courant0 = Convert.ToDouble(readed[14]);
            // Retornamos 0 si la carga ha sido exitosa.
            return 0;
        }

        // Retornamos velocidad horizontal.
        public double getU() { return u0; }
        // Retornamos velocidad vertical.
        public double getV() { return v0; }
        // Retornamos densidad.
        public double getRO() { return ro0; }
        // Retornamos presión.
        public double getP() { return p0; }
        // Retornamos temperatura.
        public double getT() { return T0; }
        // Retornamos numero Mach.
        public double getM() { return M0; }
        // Retornamos condición de courant viscosa.
        public double getCy() { return Cy0; }
        // Retornamos coeficiente adiabatico.
        public double getGamma() { return gamma0; }
        // Retornamos constante de los gases ideales.
        public double getR() { return R0; }
        // Retornamos altura inicial.
        public double getH() { return H0; }
        // Retornamos longitud máxima.
        public double getxMax() { return xMax0; }
        // Retornamos numero de divisiones verticales.
        public int getJ() { return j0; }
        // Retornamos constante de courant.
        public double getCourant() { return Courant0; }
        // Retornamos ubucacion de la esquina de expansión.
        public double getE() { return E0; }
        // Retornamos ángulo de la esquina.
        public double getTheta() { return theta0; }
    }
}
