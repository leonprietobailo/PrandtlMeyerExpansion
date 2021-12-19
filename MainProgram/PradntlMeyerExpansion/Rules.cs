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
            // Cálculo del numero Mach.
            M0 = Math.Sqrt(u0 * u0 + v0 * v0) / Math.Sqrt(gamma0 * R0 * T0);
        }

        // Almacenar los valores del objeto en un archivo en formato de texto.
        public void saveRules()
        {
            // Crear dialogo para guardar el archivo.
            SaveFileDialog diag = new SaveFileDialog();
            // Mostrar formato customizado (.sim).
            diag.Filter = "(*.sim)|*.*";
            // Formato por defecto del archivo de texto .sim.
            diag.DefaultExt = "sim";
            // Si se guarda el archivo.
            if (diag.ShowDialog() == true)
            {
                // Creamos y cerramos un nuervo archivo con el titulo espeficicado por el usuario.
                FileStream emptyFile = File.Create(diag.FileName);
                emptyFile.Close();
                // Añadimos velocidad horizontal y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(u0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos velocidad vertical y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(v0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos denisodad y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(ro0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos presion y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(p0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos temperatura y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(T0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos numero Mach y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(M0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos condicion de courant para la viscosidad y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(Cy0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos coeficiente adiabatico y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(gamma0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos constante de los gases ideales y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(R0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos ubicación de la esquina de expansión y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(E0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos ángulo y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(theta0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos divisiones verticales y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(j0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos longitud màxima y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(xMax0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos altura inicial y separamos con guión.
                File.AppendAllText(diag.FileName, Convert.ToString(H0));
                File.AppendAllText(diag.FileName, "-");
                // Añadimos condición de courant.
                File.AppendAllText(diag.FileName, Convert.ToString(Courant0));
            }
        }
        // Inicializar la clase reglas con un archivo en formato de texto.
        public int loadRules()
        {

            // Abrimos un formulario para cargar el archivo.
            OpenFileDialog dig = new OpenFileDialog();
            // Desactivamos la opcion de cargar mas de un fichero.
            dig.Multiselect = false;
            // Filtramos los archivos para ver el formato customizado .sim.
            dig.Filter = "(*.sim*)|*.*";
            // Añadimos extensión.
            dig.DefaultExt = ".sim";
            // Si hay un archvo cargado.
            if (dig.ShowDialog() == true)
            {
                // Leemos el archivo.
                StreamReader readFile = new StreamReader(dig.FileName);
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
            else
            {
                // Retornamos -1 si el usuario no ha seleccionado un archivo.
                return -1;
            }
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
