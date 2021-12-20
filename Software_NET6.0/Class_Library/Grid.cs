using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Microsoft.Win32;

namespace Class_Library
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
            // Añadimos las reglas introducidas 
            r = rIn;
            // Generamos un vector de diferencial de xi.
            dEtadX = new double[r.getJ()];
            // Añadimos a partir de que posición hay que empezar a rellenar las celdas (posicion 0 se llenarán a continuación).
            ho = 1;
            xi = 0;
            // Creamos el vector de celdas donde almacenaremos la primera columna de la malla.
            Cell[] row = new Cell[r.getJ()];
            // Creamos el vector correspondiente a la posicion vertical de las celdas de la primera columna de la malla.
            double[] rowY = new double[r.getJ()];
            // Añadimos celdas al vector.
            for (int i = 0; i < r.getJ(); i++) { row[i] = new Cell(r, i, ho); }
            // Añadimos posiciones verticales al vector.
            for (int i = 0; i < r.getJ(); i++) { rowY[i] = i * r.getH() / (r.getJ() - 1); }
            // Añadimos vectores a las listas.
            Mesh.Add(row);
            xP.Add(0);
            yP.Add(rowY);
        }

        // Cálculo de la expansión de Prandtl-Meyer.
        public void PrandtlMeyerExpansion()
        {
            // Calcularemos las magnitudes del flujo hasta que lleguemos a xMax.
            while (xi < r.getxMax())
            {
                // Creamos nuestra nueva columna.
                Cell[] row = new Cell[r.getJ()];
                // Rellenamos columna con las celdas por defecto.
                for (int i = 0; i < r.getJ(); i++) { row[i] = new Cell(r, i, ho); }
                // Añadimos la columna a la malla.
                Mesh.Add(row);
                // Inicializamos las variables necesarias para el calculo del Predictor y Corrector step.
                InitVars();
                // Calculamos el tamaño del Step (dXi) que estamos haciendo.
                ComputeStepSize();
                // Cálculo del Predictor y Corrector Step.
                for (int i = 0; i < r.getJ(); i++) { Mesh[ho][i].PredictorStep(dEtadX[i], dEta, h, dXi, this); }
                for (int i = 0; i < r.getJ(); i++) { Mesh[ho][i].CorrectorStep(dEtadX[i], dEta, h, dXi, xi, this); }
                // Sumamos uno a nuestra posicion horizontal de la malla.
                ho++;
                // Añadimos a xi el valor calculado en ComputeStepSize (dXi).
                xi += dXi;
                // Añadimos el valor actual de xi al vector de puntos horizontales.
                xP.Add(xi);
            }
        }

        // Inicialización de las variables de la columna.
        private void InitVars()
        {
            // Generamos vector donde almacenaremos las posiciones verticales de las celdas.
            double[] y = new double[r.getJ()];
            // Recorremos la columna para calcular el estado de las variables en cada una de las posiciones verticales.
            for (int i = 0; i < r.getJ(); i++)
            {
                double ys;
                // Calculamos el valor de la posicion vertical en la capa inferior y la altura de la columna.
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
                // Calculamos el diferencial de posicion vertical.
                dy = h / (r.getJ() - 1);
                // Añadimos la posicion vertical en el vecctor generado al principio de la funcion.
                y[i] = ys + dy * i;
                // Calculamos el valor de eta correspondiente a la posicion (normalización).
                double eta = (y[i] - ys) / h;
                // Calculamos diferencial de eta.
                dEta = 1.0 / (r.getJ() - 1);
                // Calculamos el valor del diferencial de Eta respecto x (o xi) y lo añadimos al vector dEtadX.
                if (xi < r.getE())
                {
                    dEtadX[i] = 0.0;
                }
                else
                {
                    dEtadX[i] = (1 - eta) * Math.Tan(r.getTheta()) / h;
                }
            }
            // Añadimos el vector de posiciones verticales a la matriz de posiciones verticales.
            yP.Add(y);
        }

        // Cálculo del tamaño del step.
        private void ComputeStepSize()
        {
            // Generamos un vector de angulos añadiendo y sustraendo mu, respectivamente.
            double[] anglePlus = new double[r.getJ()];
            double[] angleMin = new double[r.getJ()];
            double mu;
            // Recorremos columna para hallar el valor maximo de los angulos.
            for (int i = 0; i < r.getJ(); i++)
            {
                // Cálculo de mu.
                mu = Math.Asin(1 / Mesh[ho - 1][i].getM());
                // Almacenamos todos los angulos añadiendo mu.
                anglePlus[i] = Math.Abs(Math.Tan(r.getTheta() + mu));
                // Almacenamos todos los angulos sustraendo mu.
                angleMin[i] = Math.Abs(Math.Tan(r.getTheta() - mu));
            }
            // Generamos un vector unicamente para hallar los valores maximos de ambos vectores.
            double[] angleMax = new double[2];
            // Añadimos el maximo de ambos vectores.
            angleMax[0] = anglePlus.Max();
            angleMax[1] = angleMin.Max();
            // Hallamos el maximo del maximo de ambos vectores (maximo total).
            double MaximumValue = angleMax.Max();
            // Calculamos dXi dividiendo por el maximo totoal.
            dXi = r.getCourant() * dy / MaximumValue;
        }

        // Obtener el valor medio de las magnitudes fisicas en el downstream.
        public (double u, double v, double ro, double p, double T, double M) getDownstream()
        {
            // Inicialización de variables.
            double u = 0;
            double v = 0;
            double RO = 0;
            double P = 0;
            double T = 0;
            double M = 0;
            int counter = 0;
            // Recorremos los puntos del downstream segun el criterio de anderson (empezando en la segunda casilla hasta completar el 23/41 de las filas).
            for (int i = 1; i < Math.Floor(23.0 / 41.0 * r.getJ()); i++)
            {
                // Sumamos a cada variable el valor que hay en cada celda que estamos mirando.
                u += Mesh[Mesh.Count - 1][i].getU();
                v += Mesh[Mesh.Count - 1][i].getV();
                RO += Mesh[Mesh.Count - 1][i].getRO();
                P += Mesh[Mesh.Count - 1][i].getP();
                T += Mesh[Mesh.Count - 1][i].getT();
                M += Mesh[Mesh.Count - 1][i].getM();
                counter++;
            }
            // Retornamos la media de los valores.
            return (u / counter, v / counter, RO / counter, P / counter, T / counter, M / counter);
        }

        // Guardar CSV de una magnitud en concreto.
        public void saveCSV(int magnitude, string filename)
        {

            // Recorremos todas las posiciones verticales.
            for (int i = r.getJ() - 1; i >= 0; i--)
            {
                // Recorremos todas las posiciones horizontales.
                for (int n = 0; n < Mesh.Count; n++)
                {
                    // Guardamos el valor en el csv en funcion de la magnitud que queramos almacenar. (0 -> Velocidad Horizontal, 1 -> Velocidad Vertical, 2 -> Densidad, 3 -> Presion, 4 -> Temperatura y 5 -> Mach).
                    if (magnitude == 0) { File.AppendAllText(filename, Mesh[n][i].getU().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }
                    if (magnitude == 1) { File.AppendAllText(filename, Mesh[n][i].getV().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }
                    if (magnitude == 2) { File.AppendAllText(filename, Mesh[n][i].getRO().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }
                    if (magnitude == 3) { File.AppendAllText(filename, Mesh[n][i].getP().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }
                    if (magnitude == 4) { File.AppendAllText(filename, Mesh[n][i].getT().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }
                    if (magnitude == 5) { File.AppendAllText(filename, Mesh[n][i].getM().ToString(CultureInfo.CreateSpecificCulture("en-GB"))); }

                    // Si no estamos en la ultima posicion añadimos una coma para separar los valores.
                    if (n < Mesh.Count - 1) { File.AppendAllText(filename, ","); }
                    // Si estamos en la ultima posicion, hacemos un salto de linea.
                    else { File.AppendAllText(filename, "\n"); }
                }
            }
        }
        

        // Retornamos una celda en la posicion solicitada.
        public Cell GetCell(int ve, int ho) { return Mesh[ho][ve]; }
        // Retornamos la matriz de posiciones verticales.
        public List<double[]> GetYP() { return yP; }
        // Retornamos el vector de posiciones horizontales.
        public List<double> GetXP() { return xP; }
    }
}
