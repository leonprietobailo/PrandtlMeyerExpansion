using System;
using System.Collections.Generic;
using System.Data;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;
using LiveCharts;
using LiveCharts.Wpf;
using LiveCharts.Defaults;
using System.Diagnostics;

namespace PradntlMeyerExpansion
{
    public partial class MainWindow : Window
    {
        Grid mesh,meshCte;                                            // Creación de objeto Grid (mesh: caso general y modificable. meshCte: caso específio Anderson)
        Rules r,rCte;                                                 // Creación de objeto Rules (r: caso general y modificable. rCte: caso específio Anderson)
        Canvas plano = new Canvas();                                  // Creación de objeto Canvas
        double angulo, x1, y1, MaxY;                                  // Valores característicos del plano físico
        int divisionesy;                                              // Valor característico del plano físico
        double maxvalue, minvalue;                                    // Valores máximo y mínimo de cada magnitud del fluido
        Polygon[,] Polygons;                                          // Matriz que contiene todos los polígonos del plano físico
        SolidColorBrush mySolidColorBrush = new SolidColorBrush();    // Creación de un nuevo color
        SolidColorBrush color = new SolidColorBrush();                // Creación de un nuevo color
        SolidColorBrush whiteBrush = new SolidColorBrush();           // Creación de un nuevo color

        // Creación de las colecciones que contendrán los valores de las gráficas del estudio avanzado
        ChartValues<ObservablePoint> uEVOList = new ChartValues<ObservablePoint>();     // Velocidad horizontal 
        ChartValues<ObservablePoint> vEVOList = new ChartValues<ObservablePoint>();     // Velocidad vertical
        ChartValues<ObservablePoint> roEVOList = new ChartValues<ObservablePoint>();    // Densidad
        ChartValues<ObservablePoint> pEVOList = new ChartValues<ObservablePoint>();     // Presión
        ChartValues<ObservablePoint> TEVOList = new ChartValues<ObservablePoint>();     // Temperatura
        ChartValues<ObservablePoint> MEVOList = new ChartValues<ObservablePoint>();     // Número Mach

        //Creación de las series para plotear
        public SeriesCollection SeriesCollection, uCollection, vCollection, roCOllection, pCollection, TCollection, MCollection; 

        // Creación de las tablas que contendrán los resultados obtenidos por el simulador
        DataTable UTable = new DataTable();                           // Velocidad horizontal
        DataTable VTable = new DataTable();                           // Velocidad vertical
        DataTable rhoTable = new DataTable();                         // Densidad
        DataTable pTable = new DataTable();                           // Presión
        DataTable TTable = new DataTable();                           // Temperatura
        DataTable MTable = new DataTable();                           // Número Mach

        // Creación de las tablas que contendrán los resultados obtenidos por el Anderson
        DataTable AndersonUTable = new DataTable();                   // Velocidad horizontal
        DataTable AndersonVTable = new DataTable();                   // Velocidad vertical
        DataTable AndersonRhoTable = new DataTable();                 // Densidad
        DataTable AndersonpTable = new DataTable();                   // Presión
        DataTable AndersonTTable = new DataTable();                   // Temperatura
        DataTable AndersonMTable = new DataTable();                   // Número Mach

        LinearGradientBrush myHorizontalGradient = new LinearGradientBrush();

        public MainWindow()
        {
            InitializeComponent();

            mySolidColorBrush.Color = Color.FromRgb(75, 70, 70);      // Color gris
            color.Color = Color.FromRgb(82, 222, 197);                // Color negro
            whiteBrush.Color = Colors.White;                          // Color blanco

            // Cambiamos la visibilidad d elas diferentes pestañas
            Introduction.Visibility = Visibility.Visible;
            Introduction.Background = mySolidColorBrush;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
            EstudioAvanzado.Visibility = Visibility.Hidden;
            EstudioAvanzada.Visibility = Visibility.Hidden;
            RunSim.Visibility = Visibility.Hidden;
            Grid.Visibility = Visibility.Hidden;
            Info.Visibility = Visibility.Hidden;

            // Cambiamos el color del fondo de los diferentes botones del simulador (Home, Simulation, Validation, Video Tutorial, About Us)
            topic.Background = mySolidColorBrush;
            simulation.Background = Brushes.Transparent;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = Brushes.Transparent;

            // Añadimos el canvas creado anteriormente a la pestaña de simulación
            Simulation.Children.Add(plano);

            // Asociamops al gradiente los diferentes colores
            myHorizontalGradient.StartPoint = new Point(0, 0.3);
            myHorizontalGradient.EndPoint = new Point(1, 0.3);
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(240, 255, 0), 1));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(76, 175, 80), 0.5));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(117, 95, 160), 0));
            Gradient.Fill = myHorizontalGradient;

            // Creación tablas del Anderson y del simulador
            // Se introducen los valores utilizados por el Anderson en las rules rCte
            rCte = new Rules(678, 0, 1.23, 0.101e6, 286.1, 0.5, 1.4, 287, 10, 5.352 * Math.PI / 180, 41, 65, 40, 0.5);
            // Se crea un grid con las rules rCte anteriormente creadas 
            meshCte = new Grid(rCte);
            // Se llama al evento que realiza el cálculo de la expanción de Prandtl-Meyer
            meshCte.PrandtlMeyerExpansion();
            // Se llama al método que crea las tablas
            CreateTables();

            // Errores
            double uReal = 710;
            double vReal = -66.5;
            double roReal = 0.984;
            double pReal = 0.739e5;
            double TReal = 262;
            double MReal = 2.2;

            double uEVO, vEVO, roEVO, pEVO, TEVO, MEVO;
            (uEVO, vEVO, roEVO, pEVO, TEVO, MEVO) = meshCte.getDownstream();

            // Cálculo del error relativo
            double relU = Math.Round((uEVO - uReal) / uReal * 100, 4);
            double relV = Math.Round((vEVO - vReal) / vReal * 100,4);
            double relRO = Math.Round((roEVO - roReal) / roReal * 100,4);
            double relP = Math.Round((pEVO - pReal) / pReal * 100,4);
            double relT = Math.Round((TEVO - TReal) / TReal * 100,4);
            double relM = Math.Round((MEVO - MReal) / MReal * 100,4);

            // Error relativo obtendio por el Anderson
            double relUAnd = 0.45;
            double relVAnd = 3.76;
            double relROAnd = 0.813;
            double relPAnd = 1.08;
            double relTAnd = 0.038;
            double relMAnd = 0.45;

            // Poner errores del Anderson en el formulario.
            ErrorAU.Content = relUAnd.ToString();
            ErrorAV.Content = relVAnd.ToString();
            ErrorARO.Content = relROAnd.ToString();
            ErrorAP.Content = relPAnd.ToString();
            ErrorAT.Content = relTAnd.ToString();
            ErrorAM.Content = relMAnd.ToString();

            // Poner errores del Simulador en el formulario.
            ErrorSU.Content = relU.ToString();
            ErrorSV.Content = relV.ToString();
            ErrorSRO.Content = relRO.ToString();
            ErrorSP.Content = relP.ToString();
            ErrorST.Content = relT.ToString();
            ErrorSM.Content = relM.ToString();
        }

        // Evento que muestra la pestaña de Introdución y oculta el resto
        private void lntro_Click(object sender, RoutedEventArgs e)
        {
            // Se cambia la visibilidad de las diferentes pestañas
            Introduction.Visibility = Visibility.Visible;
            Introduction.Background = mySolidColorBrush;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;

            // Se cambia el color de fondo de los botones principales del simulador
            topic.Background = mySolidColorBrush;
            simulation.Background = Brushes.Transparent;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = Brushes.Transparent;
        }

        // Evento que muestra la pestaña de Simulación y oculta el resto
        private void Simualtion_Click(object sender, RoutedEventArgs e)
        {
            // Se cambia la visibilidad de las diferentes pestañas
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Visible;
            Simulation.Background = mySolidColorBrush;
            EstudioAvanzado.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
            Grids.Visibility = Visibility.Hidden;

            // Se cambia el color de fondo de los botones principales del simulador
            topic.Background = Brushes.Transparent;
            simulation.Background = mySolidColorBrush;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = Brushes.Transparent;

            if (r != null)
            {
                EstudioAvanzada.Visibility = Visibility.Visible;
                Grids.Visibility = Visibility.Visible;
                plano.Visibility = Visibility.Visible;
            }
        }

        // Evento que muestra la pestaña de Validación y oculta el resto
        private void Validation_Click(object sender, RoutedEventArgs e)
        {
            // Se cambia la visibilidad de las diferentes pestañas
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Visible;
            Validation.Background = mySolidColorBrush;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;

            // Se cambia el color de fondo de los botones principales del simulador
            topic.Background = Brushes.Transparent;
            simulation.Background = Brushes.Transparent;
            ValidationButton.Background = mySolidColorBrush;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = Brushes.Transparent;

            // Se actiba el radio button "u" para poder ver las tablas de la velocidad horizontal de manera predeterminada cuando se presiona el botón de Validación
            Utable.IsChecked = true;
        }

        // Evento que muestra la pestaña del Video Tutorial y oculta el resto
        private void vTutorial_Click(object sender, RoutedEventArgs e)
        {
            // Se cambia la visibilidad de las diferentes pestañas
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Visible;
            VideoTutorial.Background = mySolidColorBrush;
            AboutUs.Visibility = Visibility.Hidden;

            // Se cambia el color de fondo de los botones principales del simulador
            topic.Background = Brushes.Transparent;
            simulation.Background = Brushes.Transparent;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = mySolidColorBrush;
            aboutUs.Background = Brushes.Transparent;
        }

        // Evento que muestra la pestaña de About Us y oculta el resto
        private void aboutUs_Click(object sender, RoutedEventArgs e)
        {
            // Se cambia la visibilidad de las diferentes pestañas
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Visible;
            AboutUs.Background = mySolidColorBrush;

            // Se cambia el color de fondo de los botones principales del simulador
            topic.Background = Brushes.Transparent;
            simulation.Background = Brushes.Transparent;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = mySolidColorBrush;
        }

        // Evento que cambia los parámetros a los de default utilizados en el Anderson
        private void Default_Click_1(object sender, RoutedEventArgs e)
        {
            E1.Text = "10";
            x11.Text = "65";
            H1.Text = "40";
            theta1.Text = "5,352";
            dEta1.Text = "0,025";
            Cy1.Text = "0,5";
            gamma1.Text = "1,4";
            R1.Text = "287";
            C1.Text = "0,5";
            p1.Text = "10100";
            rho1.Text = "1,23";
            T1.Text = "286,1";
            u1.Text = "678";
            v1.Text = "0";
        }

        // Evento que muestra la pestaña con el plano físico de la simulación y oculta la del estudio avanzado
        private void Grid_Click(object sender, RoutedEventArgs e)
        {

            Grids.Visibility = Visibility.Visible;
            EstudioAvanzado.Visibility = Visibility.Hidden;
            EstudioAvanzada.Visibility = Visibility.Visible;
            Grid.Visibility = Visibility.Hidden;
            plano.Visibility = Visibility.Visible;
        }

        // Evento que oculta la pestaña con el plano físico de la simulación y muestra la del estudio avanzado
        private void estudio_Click(object sender, RoutedEventArgs e)
        {
            // Se cambia la visibilidad de las diferentes pestañas
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Visible;
            Simulation.Background = mySolidColorBrush;
            EstudioAvanzado.Visibility = Visibility.Visible;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
            Grids.Visibility = Visibility.Hidden;
            plano.Visibility = Visibility.Hidden;
            EstudioAvanzada.Visibility = Visibility.Hidden;
            Grid.Visibility = Visibility.Visible;
            RunSim.Visibility = Visibility.Hidden;

            // Se cambia el color de fondo de los botones principales del simulador
            topic.Background = Brushes.Transparent;
            simulation.Background = mySolidColorBrush;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = Brushes.Transparent;
        }

        // Evento que carga la frontera y las divisiones verticales del plano físico
        private void LoadGrid(object sender, RoutedEventArgs e)
        {
            try
            {
                // Se limpia el plano 
                plano.Children.Clear();

                // Se recogen todos los parámetros de simulación y se crean las reglas y el grid con ellos
                double u = Convert.ToDouble(u1.Text);
                double v = Convert.ToDouble(v1.Text);
                double rho = Convert.ToDouble(rho1.Text);
                double p = Convert.ToDouble(p1.Text);
                double T = Convert.ToDouble(T1.Text);
                double Cy = Convert.ToDouble(Cy1.Text);
                double gamma = Convert.ToDouble(gamma1.Text);
                double R = Convert.ToDouble(R1.Text);
                double E = Convert.ToDouble(E1.Text);
                double theta = Convert.ToDouble(theta1.Text);
                double dEta = Convert.ToDouble(dEta1.Text);
                int j1 = Convert.ToInt32(1 / dEta + 1);
                double xmax = Convert.ToDouble(x11.Text);
                double H = Convert.ToDouble(H1.Text);
                double C = Convert.ToDouble(C1.Text);
                r = new Rules(u, v, rho, p, T, Cy, gamma, R, E, theta * Math.PI / 180, j1, xmax, H, C);
                mesh = new Grid(r);

                // Se recogen ciertos parámetros necesarios par dibujar el plano físico
                E = r.getE();
                x1 = r.getxMax();
                y1 = r.getH();
                angulo = r.getTheta();
                divisionesy = r.getJ() - 1;

                // Se calcula el primer incremento de la coordenada y
                double ay1;
                ay1 = y1 / divisionesy;
                // Se calcula el segundo incremento de la coordenada y
                double ay2;
                ay2 = (y1 + Math.Tan(angulo) * (x1 - E)) / divisionesy;
                // Se recoge el valor máximo de la coordenada y
                MaxY = y1 + Math.Tan(angulo) * (x1 - E);


                for (int j = 0; j < divisionesy; j++)
                {
                    // Se crea un nuevo polígono
                    Polygon polygon = new Polygon();
                    // Se crean puntos con las coordenadas de cada vértice de cada polígono
                    System.Windows.Point Point11 = new System.Windows.Point(0, j * ay1 * 11 * 45.1028 / MaxY);
                    System.Windows.Point Point21 = new System.Windows.Point(E * 715 / x1, j * ay1 * 11 * 45.1028 / MaxY);
                    System.Windows.Point Point31 = new System.Windows.Point(x1 * 715 / x1, j * ay2 * 11 * 45.1028 / MaxY);
                    System.Windows.Point Point41 = new System.Windows.Point(x1 * 715 / x1, (j + 1) * ay2 * 11 * 45.1028 / MaxY);
                    System.Windows.Point Point51 = new System.Windows.Point(E * 715 / x1, (j + 1) * ay1 * 11 * 45.1028 / MaxY);
                    System.Windows.Point Point61 = new System.Windows.Point(0, (j + 1) * ay1 * 11 * 45.1028 / MaxY);
                    // Se crean las colleciones de puntos
                    PointCollection polygonPoints1 = new PointCollection();
                    // Se añaden los puntos anteriormente creados a la coleccion de puntos
                    polygonPoints1.Add(Point11);
                    polygonPoints1.Add(Point21);
                    polygonPoints1.Add(Point31);
                    polygonPoints1.Add(Point41);
                    polygonPoints1.Add(Point51);
                    polygonPoints1.Add(Point61);
                    // Se asocia la colección de puntos anteriormente creada a los puntos del polígono
                    polygon.Points = polygonPoints1;
                    // Se añade el color del borde del polígono y se le asocia un grosor
                    polygon.Stroke = whiteBrush;
                    polygon.StrokeThickness = 0.1;
                    // Se crea un margen para el poano físico
                    System.Windows.Thickness planoPM = new System.Windows.Thickness(520, 100, 92, 292);
                    plano.Margin = planoPM;
                    // Se añaden los polígonos al plano físico
                    plano.Children.Add(polygon);
                    // Se coloca el canvas
                    Canvas.SetTop(polygon, 0);
                    Canvas.SetLeft(polygon, 0);
                }
                // Se cambia la visibilidad de las diferentes pestañas
                RunSim.Visibility = Visibility.Visible;
                Grid.Visibility = Visibility.Hidden;
                EstudioAvanzada.Visibility = Visibility.Hidden;
                Grids.Visibility = Visibility.Hidden;
                plano.Visibility = Visibility.Visible;
                EstudioAvanzado.Visibility = Visibility.Hidden;
                Error_Label.Visibility = Visibility.Hidden;
            }
            catch
            {
                Error_Label.Visibility = Visibility.Visible;
                RunSim.Visibility = Visibility.Hidden;
            }
        }

        // Evento que carga el plano físico
        private void Run_Click(object sender, RoutedEventArgs e)
        {
            // Se vacian las listas del estudio avanzado para volver a realizar los calculos.
            uEVOList.Clear();
            vEVOList.Clear();
            roEVOList.Clear();
            pEVOList.Clear();
            TEVOList.Clear();
            MEVOList.Clear();
            // Se cambia la visibilidad de las diferentes pestañas
            Grid.Visibility = Visibility.Hidden;
            EstudioAvanzada.Visibility = Visibility.Visible;
            Grids.Visibility = Visibility.Visible;

            // Se calcula la expansión de Prandtl-Meyer
            plano.Children.Clear();
            mesh.PrandtlMeyerExpansion();

            // Se desactiva el radio button "u"
            uRB.IsChecked = false;

            // Se crea un nuevo vector de polígonos
            Polygons = new Polygon[mesh.GetXP().Count - 1, divisionesy];


            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    // Se crea un nuevo polígono
                    Polygon polygon = new Polygon();
                    // Se crean puntos con las coordenadas de cada vértice de cada polígono
                    System.Windows.Point Point11 = new System.Windows.Point(mesh.GetXP()[i] * 715 / x1, (mesh.GetYP()[0][40] - mesh.GetYP()[i][j]) * 11 * 45.1028 / MaxY);
                    System.Windows.Point Point21 = new System.Windows.Point(mesh.GetXP()[i + 1] * 715 / x1, (mesh.GetYP()[0][40] - mesh.GetYP()[i + 1][j]) * 11 * 45.1028 / MaxY);
                    System.Windows.Point Point31 = new System.Windows.Point(mesh.GetXP()[i + 1] * 715 / x1, (mesh.GetYP()[0][40] - mesh.GetYP()[i + 1][j + 1]) * 11 * 45.1028 / MaxY);
                    System.Windows.Point Point41 = new System.Windows.Point(mesh.GetXP()[i] * 715 / x1, (mesh.GetYP()[0][40] - mesh.GetYP()[i][j + 1]) * 11 * 45.1028 / MaxY);
                    // Se crean las colleciones de puntos
                    PointCollection polygonPoints1 = new PointCollection();
                    // Se añaden los puntos anteriormente creados a la coleccion de puntos
                    polygonPoints1.Add(Point11);
                    polygonPoints1.Add(Point21);
                    polygonPoints1.Add(Point31);
                    polygonPoints1.Add(Point41);
                    // Se asocia la colección de puntos anteriormente creada a los puntos del polígono
                    polygon.Points = polygonPoints1;
                    // Se añade el color del borde del polígono y se le asocia un grosor
                    polygon.Stroke = whiteBrush;
                    polygon.StrokeThickness = 0.1;
                    // Se crea un margen para el poano físico
                    System.Windows.Thickness planoPM = new System.Windows.Thickness(520, 100, 92, 292);
                    plano.Margin = planoPM;
                    // Se añaden los polígonos al plano físico
                    plano.Children.Add(polygon);
                    // Se añade un tag  a cada polígono creado
                    polygon.Tag = new Point(i, j);
                    // Se asocia a cada polígono un evento cuando el usuario pone el mouse encima del plano físico
                    polygon.MouseEnter += new MouseEventHandler(polygon_MouseEnter);
                    // Se asocia a cada polígono un evento cuando el usuario pone el mouse fuera del plano físico
                    polygon.MouseLeave += new MouseEventHandler(polygon_MouseLeave);
                    // Se añade cada polígono al vector de polígonos
                    Polygons[i, j] = polygon;
                }
            }
            // Se activa el radio button "u" para poder ver el plano físico con los valores de la velocidad horizontal por defecto
            uRB.IsChecked = true;
            computeEvolutionChange();
            setChartNumbers();
            
        }

        // Evento que colorea los polígonos del plano físico de acuerdo a los valores de la velocidad horizontal
        private void u_Checked(object sender, RoutedEventArgs e)
        {
            // Asociamos un valor máximo y mínimo teniendo un gran margen para poder calcularlos posteriormente
            maxvalue = -100000000;
            minvalue = 100000000;
            // Calculamos los valores máximo y mínimo de la velocidad horizontal
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = 0; j < divisionesy; j++)
                {
                    if (mesh.GetCell(j, i).getU() > maxvalue)
                    {
                        maxvalue = mesh.GetCell(j, i).getU();
                    }

                    if (mesh.GetCell(j, i).getU() < minvalue)
                    {
                        minvalue = mesh.GetCell(j, i).getU();
                    }
                }
            }
            // Coloreamos cada polígono 
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    // Se obtiene una media con los valores de las cuatro esquinas del polígono
                    double mediaM = (mesh.GetCell(j, i).getU() + mesh.GetCell(j, i + 1).getU() + mesh.GetCell(j + 1, i + 1).getU() + mesh.GetCell(j + 1, i).getU()) / 4.0;
                    
                    // Se obtiene el color asocado a R
                    byte first = 0;
                    if (mediaM < ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (41 * (mediaM - (maxvalue + minvalue) / 2)) / (minvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM > ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (164 * (mediaM - (maxvalue + minvalue) / 2)) / (maxvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM == ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(76);
                    }
                    // Se obtiene el color asocado a G
                    byte second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    // Se obtiene el color asocado a B
                    byte third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));
                    // Se crea el color RGB con los valores anteriormente calculados
                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    // Se colorea el interior de cada polígono
                    Polygons[i, j].Fill = mySolidColorBrush;
                }
            }
            // Se muestra el valor maxímo en la leyenda del gradiente
            minValueGra.Content = Convert.ToString(Math.Round(maxvalue) + " m/s");
            // Se muestra el valor mínimo en la leyenda del gradiente
            maxValueGra.Content = Convert.ToString(Math.Round(minvalue) + " m/s");
        }

        // Evento que colorea los polígonos del plano físico de acuerdo a los valores de la velocidad vertical
        private void v_Checked(object sender, RoutedEventArgs e)
        {
            // Asociamos un valor máximo y mínimo teniendo un gran margen para poder calcularlos posteriormente
            maxvalue = -100000000;
            minvalue = 100000000;
            // Calculamos los valores máximo y mínimo de la velocidad vertical
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = 0; j < divisionesy; j++)
                {
                    if (mesh.GetCell(j, i).getV() > maxvalue)
                    {
                        maxvalue = mesh.GetCell(j, i).getV();
                    }

                    if (mesh.GetCell(j, i).getV() < minvalue)
                    {
                        minvalue = mesh.GetCell(j, i).getV();
                    }
                }
            }
            // Coloreamos cada polígono 
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    // Se obtiene una media con los valores de las cuatro esquinas del polígono
                    double mediaM = (mesh.GetCell(j, i).getV() + mesh.GetCell(j, i + 1).getV() + mesh.GetCell(j + 1, i + 1).getV() + mesh.GetCell(j + 1, i).getV()) / 4.0;

                    // Se obtiene el color asocado a R
                    byte first = 0;
                    if (mediaM < ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (41 * (mediaM - (maxvalue + minvalue) / 2)) / (minvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM > ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (164 * (mediaM - (maxvalue + minvalue) / 2)) / (maxvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM == ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(76);
                    }
                    // Se obtiene el color asocado a G
                    byte second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    // Se obtiene el color asocado a B
                    byte third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));
                    // Se crea el color RGB con los valores anteriormente calculados
                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    // Se colorea el interior de cada polígono
                    Polygons[i, j].Fill = mySolidColorBrush;
                }
            }
            // Se muestra el valor maxímo en la leyenda del gradiente
            minValueGra.Content = Convert.ToString(Math.Round(maxvalue) + " m/s");
            // Se muestra el valor mínimo en la leyenda del gradiente
            maxValueGra.Content = Convert.ToString(Math.Round(minvalue) + " m/s");
        }

        // Evento que colorea los polígonos del plano físico de acuerdo a los valores de la densidad
        private void rho_Checked(object sender, RoutedEventArgs e)
        {
            // Asociamos un valor máximo y mínimo teniendo un gran margen para poder calcularlos posteriormente
            maxvalue = -100000000;
            minvalue = 100000000;

            // Calculamos los valores máximo y mínimo de la densidad
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = 0; j < divisionesy; j++)
                {
                    if (mesh.GetCell(j, i).getRO() > maxvalue)
                    {
                        maxvalue = mesh.GetCell(j, i).getRO();
                    }

                    if (mesh.GetCell(j, i).getRO() < minvalue)
                    {
                        minvalue = mesh.GetCell(j, i).getRO();
                    }
                }
            }

            // Coloreamos cada polígono 
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    // Se obtiene una media con los valores de las cuatro esquinas del polígono
                    double mediaM = (mesh.GetCell(j, i).getRO() + mesh.GetCell(j, i + 1).getRO() + mesh.GetCell(j + 1, i + 1).getRO() + mesh.GetCell(j + 1, i).getRO()) / 4.0;

                    // Se obtiene el color asocado a R
                    byte first = 0;
                    if (mediaM < ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (41 * (mediaM - (maxvalue + minvalue) / 2)) / (minvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM > ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (164 * (mediaM - (maxvalue + minvalue) / 2)) / (maxvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM == ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(76);
                    }
                    // Se obtiene el color asocado a G
                    byte second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    // Se obtiene el color asocado a B
                    byte third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));
                    // Se crea el color RGB con los valores anteriormente calculados
                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    // Se colorea el interior de cada polígono
                    Polygons[i, j].Fill = mySolidColorBrush;
                }
            }
            // Se muestra el valor maxímo en la leyenda del gradiente
            minValueGra.Content = Convert.ToString(Math.Round(maxvalue, 4) + " kg/m3");
            // Se muestra el valor mínimo en la leyenda del gradiente
            maxValueGra.Content = Convert.ToString(Math.Round(minvalue, 4) + " kg/m3");
        }

        // Evento que colorea los polígonos del plano físico de acuerdo a los valores de la presión
        private void p_Checked(object sender, RoutedEventArgs e)
        {
            // Asociamos un valor máximo y mínimo teniendo un gran margen para poder calcularlos posteriormente
            double maxvalue = -100000000;
            double minvalue = 100000000;

            // Calculamos los valores máximo y mínimo de la presión
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = 0; j < divisionesy; j++)
                {
                    if (mesh.GetCell(j, i).getP() > maxvalue)
                    {
                        maxvalue = mesh.GetCell(j, i).getP();
                    }

                    if (mesh.GetCell(j, i).getP() < minvalue)
                    {
                        minvalue = mesh.GetCell(j, i).getP();
                    }
                }
            }

            // Coloreamos cada polígono 
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    // Se obtiene una media con los valores de las cuatro esquinas del polígono
                    double mediaM = (mesh.GetCell(j, i).getP() + mesh.GetCell(j, i + 1).getP() + mesh.GetCell(j + 1, i + 1).getP() + mesh.GetCell(j + 1, i).getP()) / 4.0;
                    
                    // Se obtiene el color asocado a R
                    byte first = 0;
                    if (mediaM < ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (41 * (mediaM - (maxvalue + minvalue) / 2)) / (minvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM > ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (164 * (mediaM - (maxvalue + minvalue) / 2)) / (maxvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM == ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(76);
                    }
                    // Se obtiene el color asocado a G
                    byte second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    // Se obtiene el color asocado a B
                    byte third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));
                    // Se crea el color RGB con los valores anteriormente calculados
                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    // Se colorea el interior de cada polígono
                    Polygons[i, j].Fill = mySolidColorBrush;
                }
            }
            // Se muestra el valor maxímo en la leyenda del gradiente
            minValueGra.Content = Convert.ToString(Math.Round(maxvalue) + " Pa");
            // Se muestra el valor mínimo en la leyenda del gradiente
            maxValueGra.Content = Convert.ToString(Math.Round(minvalue) + " Pa");
        }

        // Evento que colorea los polígonos del plano físico de acuerdo a los valores de la temperatura
        private void T_Checked(object sender, RoutedEventArgs e)
        {
            // Asociamos un valor máximo y mínimo teniendo un gran margen para poder calcularlos posteriormente
            double maxvalue = -100000000;
            double minvalue = 100000000;

            // Calculamos los valores máximo y mínimo de la temperatura
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = 0; j < divisionesy; j++)
                {
                    if (mesh.GetCell(j, i).getT() > maxvalue)
                    {
                        maxvalue = mesh.GetCell(j, i).getT();
                    }

                    if (mesh.GetCell(j, i).getT() < minvalue)
                    {
                        minvalue = mesh.GetCell(j, i).getT();
                    }
                }
            }

            // Coloreamos cada polígono 
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    // Se obtiene una media con los valores de las cuatro esquinas del polígono
                    double mediaM = (mesh.GetCell(j, i).getT() + mesh.GetCell(j, i + 1).getT() + mesh.GetCell(j + 1, i + 1).getT() + mesh.GetCell(j + 1, i).getT()) / 4.0;

                    // Se obtiene el color asocado a R
                    byte first = 0;
                    if (mediaM < ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (41 * (mediaM - (maxvalue + minvalue) / 2)) / (minvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM > ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (164 * (mediaM - (maxvalue + minvalue) / 2)) / (maxvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM == ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(76);
                    }
                    // Se obtiene el color asocado a G
                    byte second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    // Se obtiene el color asocado a B
                    byte third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));
                    // Se crea el color RGB con los valores anteriormente calculados
                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    // Se colorea el interior de cada polígono
                    Polygons[i, j].Fill = mySolidColorBrush;
                }
            }
            // Se muestra el valor mínimo en la leyenda del gradiente
            maxValueGra.Content = Convert.ToString(Math.Round(minvalue) + " K");
            // Se muestra el valor maxímo en la leyenda del gradiente
            minValueGra.Content = Convert.ToString(Math.Round(maxvalue) + " K");
        }

        // Evento que colorea los polígonos del plano físico de acuerdo a los valores del número de Mach
        private void M_Checked(object sender, RoutedEventArgs e)
        {
            // Asociamos un valor máximo y mínimo teniendo un gran margen para poder calcularlos posteriormente
            double maxvalue = -100000000;
            double minvalue = 100000000;

            // Calculamos los valores máximo y mínimo del número de Mach
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = 0; j < divisionesy; j++)
                {
                    if (mesh.GetCell(j, i).getM() > maxvalue)
                    {
                        maxvalue = mesh.GetCell(j, i).getM();
                    }

                    if (mesh.GetCell(j, i).getM() < minvalue)
                    {
                        minvalue = mesh.GetCell(j, i).getM();
                    }
                }
            }

            // Coloreamos cada polígono 
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    // Se obtiene una media con los valores de las cuatro esquinas del polígono
                    double mediaM = (mesh.GetCell(j, i).getM() + mesh.GetCell(j, i + 1).getM() + mesh.GetCell(j + 1, i + 1).getM() + mesh.GetCell(j + 1, i).getM()) / 4.0;

                    // Se obtiene el color asocado a R
                    byte first = 0;
                    if (mediaM < ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (41 * (mediaM - (maxvalue + minvalue) / 2)) / (minvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM > ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(Math.Round(76 + (164 * (mediaM - (maxvalue + minvalue) / 2)) / (maxvalue - (maxvalue + minvalue) / 2)));
                    }
                    else if (mediaM == ((maxvalue + minvalue) / 2))
                    {
                        first = Convert.ToByte(76);
                    }
                    // Se obtiene el color asocado a G
                    byte second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    // Se obtiene el color asocado a B
                    byte third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));
                    // Se crea el color RGB con los valores anteriormente calculados
                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    // Se colorea el interior de cada polígono
                    Polygons[i, j].Fill = mySolidColorBrush;
                }
            }
            // Se muestra el valor mínimo en la leyenda del gradiente
            maxValueGra.Content = Convert.ToString(Math.Round(minvalue, 4));
            // Se muestra el valor maxímo en la leyenda del gradiente
            minValueGra.Content = Convert.ToString(Math.Round(maxvalue, 4));
        }

        // Evento que muestra los valores de las diferentes magnitudes físicas cuando el usuario pone el mouse sobre alguna celda del plano físico
        private void polygon_MouseEnter(object sender, MouseEventArgs e)
        {
            // Obtenemos la ubicación del puntero
            Polygon pol = (Polygon)sender;
            Point p = (Point)pol.Tag;
            // Cambiamos la visibilidad de las labels
            Info.Visibility = Visibility.Visible;
            // Mostramos los valores de las diferentes magnitudes físicas del fluido en las labels
            uData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getU(), 2) + " m/s";
            vData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getV(), 2) + " m/s";
            rhoData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getRO(), 2) + " kg/m3";
            pData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getP(), 2) + " Pa";
            TData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getT(), 2) + " K";
            MData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getM(), 2);
            xData.Content = p.X + 1;
            yData.Content = p.Y + 1;
        }

        // Evento que oculta los valores de las diferentes magnitudes físicas cuando el usuario retira el mouse del plano físico
        private void polygon_MouseLeave(object sender, MouseEventArgs e)
        {
            // Cambiamos la visibilidad de las labels
            Info.Visibility = Visibility.Hidden;
        }

        // Creación de las tablas para la validación
        private void CreateTables()
        {
            int divisionesycte = rCte.getJ();       // Otenemos el número de filas que tendrá la tabla

            // u Table
            UTable.Columns.Add("y-x");              // Añadimos la  primera columna a la tabla de resultados del simulador
            // u Anderson Table
            AndersonUTable.Columns.Add("y-x");      // Añadimos la  primera columna a la tabla del Anderson de u
            AndersonUTable.Columns.Add("18");       // Añadimos la  segunda columna a la tabla del Anderson de u
            AndersonUTable.Columns.Add("89");       // Añadimos la  tercera columna a la tabla del Anderson de u

            // Introducimos en dos listas los valores del Anderson de u
            List<double> uAcolumna18 = new List<double> { 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 679, 683, 691, 701, 707 };
            List<double> uAcolumna89 = new List<double> { 678, 678, 678, 679, 679, 680, 681, 683, 685, 688, 690, 693, 696, 699, 702, 705, 707, 709, 711, 713, 713, 713, 712, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 710, 705 };

            // v Table
            VTable.Columns.Add("y-x");
            // v Anderson Table
            AndersonVTable.Columns.Add("y-x");
            AndersonVTable.Columns.Add("18");
            AndersonVTable.Columns.Add("89");
            List<double> vAcolumna18 = new List<double> { 0, -0.636e-4, -0.107e-3, -0.342e-4, -0.128e-3, -0.848e-4, 0.401e-4, 0.161e-3, 0.160e-3, .242e-03, -.607e-04, -.193e-04, .125e-03, .354e-05, .120e-03, .118e-03, .217e-10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -.325e-13, -.642e-04, -.598e-04, .180e-04, -.195e-04, -.702e-04, .472e-04, -.167e-03, .326e-05, -.148e-1, -.131e+1, -.869e+1, -.266e+02, -.494e+02, -.662e+02 };
            List<double> vAcolumna89 = new List<double> { 0, -0.229, -0.499, -0.105e1, -0.203e1, -0.361e01, -0.591e+1, -0.901e+1, -0.129e+2, -0.175e+2, -0.227e+2, -0.283e+2, -0.343e+2, -0.405e+2, -0.468e+2, -0.531e+2, -0.591e+2, -0.647e+2, -0.693e+2, -0.726e+2, -0.740e+2, -0.732e+2, -0.708e+2, -0.683e+2, -0.672e+2, -0.678e+2, -0.690e+2, -0.696e+2, -0.694e+2, -0.688e+2, -0.686e+2, -0.688e+2, -0.690e+2, -0.690e+2, -0.689e+2, -0.688e+2, -0.689e+2, -0.688e+2, -0.690e+2, -0.682e+2, -0.661e+2 };

            // rho Table
            rhoTable.Columns.Add("y-x");
            // rho Anderson Table
            AndersonRhoTable.Columns.Add("y-x");
            AndersonRhoTable.Columns.Add("18");
            AndersonRhoTable.Columns.Add("89");
            List<double> rhoAcolumna18 = new List<double> { 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.123e+1, 0.122e+1, 0.119e+1, 0.112e+1, 0.104e+1, 0.992 };
            List<double> rhoAcolumna89 = new List<double> { 0.123e1, 0.123e1, 0.123e1, 0.123e1, 0.122e1, 0.121e1, 0.121e1, 0.129e1, 0.118e1, 0.116e1, 0.114e1, 0.112e1, 0.110e1, 0.107e1, 0.105e1, 0.103e1, 0.101e1, 0.990, 0.975, 0.974, 0.960, 0.963, 0.970, 0.978, 0.982, 0.980, 0.986, 0.974, 0.975, 0.977, 0.977, 0.977, .976, .976, .976, .976, .976, .977, .979, 0.107e+1, 0.109e+1 };

            // p Table
            pTable.Columns.Add("y-x");
            // p Anderson Table 
            AndersonpTable.Columns.Add("y-x");
            AndersonpTable.Columns.Add("18");
            AndersonpTable.Columns.Add("89");
            List<double> pAcolumna18 = new List<double> { 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.787e5, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.101e6, 0.100e6, 0.969e5, 0.891e5, 0.795e5, 0.734e5 };
            List<double> pAcolumna89 = new List<double> { 0.101e6, 0.101e6, 0.101e6, 0.100e6, 0.100e6, 0.993e5, 0.982e5, 0.968e5, 0.950e5, 0.930e5, 0.907e5, 0.883e5, 0.859e5, 0.834e5, 0.810e5, 0.787e5, 0.765e5, 0.746e5, 0.730e5, 0.719e5, 0.714e5, 0.717e5, 0.725e5, 0.733e5, 0.737e5, 0.735e5, 0.731e5, 0.729e5, 0.729e5, 0.731e5, 0.732e5, 0.731e5, 0.731e5, 0.731e5, 0.731e5, 0.731e5, 0.731e5, 0.731e5, 0.732e5, 0.730e5, 0.731e5 };

            // T Table
            TTable.Columns.Add("y-x");
            // T Anderson Table 
            AndersonTTable.Columns.Add("y-x");
            AndersonTTable.Columns.Add("18");
            AndersonTTable.Columns.Add("89");
            List<double> TAcolumna18 = new List<double> { 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.283e3, 0.277e3, 0.267e3, 0.258e3 };
            List<double> TAcolumna89 = new List<double> { 0.286e3, 0.286e3, 0.286e3, 0.286e3, 0.285e3, 0.285e3, 0.284e3, 0.283e3, 0.281e3, 0.286e3, 0.279e3, 0.277e3, 0.275e3, 0.273e3, 0.271e3, 0.269e3, 0.266e3, 0.264e3, 0.262e3, 0.261e3, 0.260e3, 0.259e3, 0.259e3, 0.260e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.261e3, 0.263e3, 0.237e3, 0.233e3 };

            // M Table
            MTable.Columns.Add("y-x");
            // M Anderson Table
            AndersonMTable.Columns.Add("y-x");
            AndersonMTable.Columns.Add("18");
            AndersonMTable.Columns.Add("89");
            List<double> MAcolumna18 = new List<double> { 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.203e1, 0.208e1, 0.215e1, 0.220e1 };
            List<double> MAcolumna89 = new List<double> { 0.200e1, 0.200e1, 0.200e1, 0.200e1, 0.201e1, 0.201e1, 0.202e1, 0.203e1, 0.204e1, 0.205e1, 0.207e1, 0.209e1, 0.210e1, 0.212e1, 0.214e1, 0.216e1, 0.218e1, 0.219e1, 0.221e1, 0.222e1, 0.222e1, 0.222e1, 0.221e1, 0.221e1, 0.220e1, 0.220e1, 0.221e1, 0.221e1, 0.221e1, 0.221e1, 0.221e1, 0.221e1, 0.221e1, 0.221e1, 0.221e1, 0.221e1, 0.221e1, 0.221e1, 0.220e1, 0.231e1, 0.231e1 };

            // Añadimos todas las columnas a cada una de las tablas de las magnitudes del fluido
            for (int i = 0; i < meshCte.GetXP().Count; i++)
            {
                // u Table
                UTable.Columns.Add(Convert.ToString(i + 1));
                // v Table
                VTable.Columns.Add(Convert.ToString(i + 1));
                // rho Table
                rhoTable.Columns.Add(Convert.ToString(i + 1));
                // p Table
                pTable.Columns.Add(Convert.ToString(i + 1));
                // T Table
                TTable.Columns.Add(Convert.ToString(i + 1));
                // M Table
                MTable.Columns.Add(Convert.ToString(i + 1));
            }
            int s = 0;
            for (int j = divisionesycte - 1; j >= 0; j--)
            {
                // u Table
                // Añadimos el número de las filas a la tabla de resultados del simulador
                DataRow uRow = UTable.NewRow();
                uRow["y-x"] = Convert.ToString(j + 1);
                UTable.Rows.Add(uRow);
                // u Anderson Table
                // Añadimos los valores de las listas creadas anteriormente a las tablas del Anderson
                DataRow AuRow = AndersonUTable.NewRow();
                AuRow["y-x"] = Convert.ToString(j + 1);
                AuRow["18"] = Convert.ToString(uAcolumna18[s]);
                AuRow["89"] = Convert.ToString(uAcolumna89[s]);
                AndersonUTable.Rows.Add(AuRow);

                // v Table
                DataRow vRow = VTable.NewRow();
                vRow["y-x"] = Convert.ToString(j + 1);
                VTable.Rows.Add(vRow);
                // v Anderson Table
                DataRow AvRow = AndersonVTable.NewRow();
                AvRow["y-x"] = Convert.ToString(j + 1);
                AvRow["18"] = Convert.ToString(vAcolumna18[s]);
                AvRow["89"] = Convert.ToString(vAcolumna89[s]);
                AndersonVTable.Rows.Add(AvRow);

                // rho Table
                DataRow rhoRow = rhoTable.NewRow();
                rhoRow["y-x"] = Convert.ToString(j + 1);
                rhoTable.Rows.Add(rhoRow);
                // rho Anderson Table
                DataRow ArhoRow = AndersonRhoTable.NewRow();
                ArhoRow["y-x"] = Convert.ToString(j + 1);
                ArhoRow["18"] = Convert.ToString(rhoAcolumna18[s]);
                ArhoRow["89"] = Convert.ToString(rhoAcolumna89[s]);
                AndersonRhoTable.Rows.Add(ArhoRow);

                // p Table
                DataRow pRow = pTable.NewRow();
                pRow["y-x"] = Convert.ToString(j + 1);
                pTable.Rows.Add(pRow);
                // p Anderson Table
                DataRow ApRow = AndersonpTable.NewRow();
                ApRow["y-x"] = Convert.ToString(j + 1);
                ApRow["18"] = Convert.ToString(pAcolumna18[s]);
                ApRow["89"] = Convert.ToString(pAcolumna89[s]);
                AndersonpTable.Rows.Add(ApRow);

                // T Table
                DataRow TRow = TTable.NewRow();
                TRow["y-x"] = Convert.ToString(j + 1);
                TTable.Rows.Add(TRow);
                // T Anderson Table
                DataRow ATRow = AndersonTTable.NewRow();
                ATRow["y-x"] = Convert.ToString(j + 1);
                ATRow["18"] = Convert.ToString(TAcolumna18[s]);
                ATRow["89"] = Convert.ToString(TAcolumna89[s]);
                AndersonTTable.Rows.Add(ATRow);

                // M Table
                DataRow MRow = MTable.NewRow();
                MRow["y-x"] = Convert.ToString(j + 1);
                MTable.Rows.Add(MRow);
                // M Anderson Table
                DataRow AMRow = AndersonMTable.NewRow();
                AMRow["y-x"] = Convert.ToString(j + 1);
                AMRow["18"] = Convert.ToString(MAcolumna18[s]);
                AMRow["89"] = Convert.ToString(MAcolumna89[s]);
                AndersonMTable.Rows.Add(AMRow);
                s++;
                for (int i = 0; i < (meshCte.GetXP().Count); i++)
                {
                    // u Table
                    // Rellenamos la tabla con los valores obtenidos por el simualador
                    uRow[Convert.ToString(i + 1)] = Convert.ToString(Math.Round(meshCte.GetCell(j, i).getU(), 4));
                    // v Table
                    vRow[Convert.ToString(i + 1)] = Convert.ToString(Math.Round(meshCte.GetCell(j, i).getV(), 4));
                    // rho Table
                    rhoRow[Convert.ToString(i + 1)] = Convert.ToString(Math.Round(meshCte.GetCell(j, i).getRO(), 4));
                    // p Table
                    pRow[Convert.ToString(i + 1)] = Convert.ToString(Math.Round(meshCte.GetCell(j, i).getP(), 4));
                    // T Table
                    TRow[Convert.ToString(i + 1)] = Convert.ToString(Math.Round(meshCte.GetCell(j, i).getT(), 4));
                    // M Table
                    MRow[Convert.ToString(i + 1)] = Convert.ToString(Math.Round(meshCte.GetCell(j, i).getM(), 4));
                }
            }
        }

        // Evento que muestra las tablas de la velocidad horizontal cuando se presiona el radio button "u"
        private void uTable_Checked(object sender, RoutedEventArgs e)
        {
            // Se muestra la tabla
            gridData.Visibility = Visibility.Visible;
            // Se pasa la información del data table al DataGrid
            gridData.DataContext = UTable.DefaultView;
            gridAndersonData.DataContext = AndersonUTable.DefaultView;
        }

        // Evento que muestra las tablas de la velocidad vertical cuando se presiona el radio button "v"
        private void vTable_Checked(object sender, RoutedEventArgs e)
        {
            // Se pasa la información del data table al DataGrid
            gridData.DataContext = VTable.DefaultView;
            gridAndersonData.DataContext = AndersonVTable.DefaultView;
        }

        // Evento que muestra las tablas de la densidad cuando se presiona el radio button "ρ"
        private void rhoTable_Checked(object sender, RoutedEventArgs e)
        {
            // Se pasa la información del data table al DataGrid
            gridData.DataContext = rhoTable.DefaultView;
            gridAndersonData.DataContext = AndersonRhoTable.DefaultView;
        }

        // Evento que muestra las tablas de la presión cuando se presiona el radio button "p"
        private void pTable_Checked(object sender, RoutedEventArgs e)
        {
            // Se pasa la información del data table al DataGrid
            gridData.DataContext = pTable.DefaultView;
            gridAndersonData.DataContext = AndersonpTable.DefaultView;
        }

        // Evento que muestra las tablas de la temperatura cuando se presiona el radio button "T"
        private void TTable_Checked(object sender, RoutedEventArgs e)
        {
            // Se pasa la información del data table al DataGrid
            gridData.DataContext = TTable.DefaultView;
            gridAndersonData.DataContext = AndersonTTable.DefaultView;
        }

        // Evento que muestra las tablas del Mach Number cuando se presiona el radio button "M"
        private void MTable_Checked(object sender, RoutedEventArgs e)
        {
            // Se pasa la información del data table al DataGrid
            gridData.DataContext = MTable.DefaultView;
            gridAndersonData.DataContext = AndersonMTable.DefaultView;
        }


        private void selected_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            gridData.SelectedIndex = gridAndersonData.SelectedIndex;
        }

        private void gridData_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            gridAndersonData.SelectedIndex = gridData.SelectedIndex;
        }


        

        private void ShowOnGM_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                r = new Rules();
                int result = r.loadRules();
                if (result == 0)
                {
                    Run_Click(null, null);
                    Error_Label.Visibility = Visibility.Hidden;
                    RunSim.Visibility = Visibility.Visible;
                }
                //else
                //{
                //    RunSim.Visibility = Visibility.Hidden;
                //}
            }
            catch
            {
                Error_Label.Visibility = Visibility.Visible;
                RunSim.Visibility = Visibility.Hidden;
            }
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            r.saveRules();
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            int mode = 0;
            if (Convert.ToBoolean(uRB.IsChecked)) { mode = 0; }
            else if (Convert.ToBoolean(vRB.IsChecked)) { mode = 1; }
            else if (Convert.ToBoolean(roRB.IsChecked)) { mode = 2; }
            else if (Convert.ToBoolean(pRB.IsChecked)) { mode = 3; }
            else if (Convert.ToBoolean(TRB.IsChecked)) { mode = 4; }
            else if (Convert.ToBoolean(MRB.IsChecked)) { mode = 5; }
            mesh.saveCSV(mode);
        }

        private void computeEvolutionChange()
        {
            for(double theta = 10; theta >= 0; theta -= 0.1) // De - 10 grados a 0 grados.
            {
                double uEVO, vEVO, roEVO, pEVO, TEVO, MEVO;
                Rules rEVO = new Rules(r, theta * Math.PI / 180);
                Grid gEVO = new Grid(rEVO);
                gEVO.PrandtlMeyerExpansion();
                (uEVO, vEVO, roEVO, pEVO, TEVO, MEVO) = gEVO.getDownstream();
                uEVOList.Add(new ObservablePoint(-Convert.ToDouble(Math.Round(theta, 4)), Convert.ToDouble(Math.Round(uEVO, 4))));
                vEVOList.Add(new ObservablePoint(-Convert.ToDouble(Math.Round(theta, 4)), Convert.ToDouble(Math.Round(vEVO, 4))));
                roEVOList.Add(new ObservablePoint(-Convert.ToDouble(Math.Round(theta, 4)), Convert.ToDouble(Math.Round(roEVO, 4))));
                pEVOList.Add(new ObservablePoint(-Convert.ToDouble(Math.Round(theta, 4)), Convert.ToDouble(Math.Round(pEVO, 4))));
                TEVOList.Add(new ObservablePoint(-Convert.ToDouble(Math.Round(theta, 4)), Convert.ToDouble(Math.Round(TEVO, 4))));
                MEVOList.Add(new ObservablePoint(-Convert.ToDouble(Math.Round(theta, 4)), Convert.ToDouble(Math.Round(MEVO, 4))));
            }
            uChartAxis.Value = r.getU();
            vChartAxis.Value = r.getV();
            RoChartAxis.Value = r.getRO();
            pChartAxis.Value = r.getP();
            TchartAxis.Value = r.getT();
            MChartAxis.Value = r.getM();
        }

        

        
        

        private void setChartNumbers()
        {
            uCollection = new SeriesCollection
            {
                new LineSeries
                {
                    Values = uEVOList,
                    ScalesYAt = 0,
                    Stroke = new SolidColorBrush(Color.FromRgb(75, 75, 255)),
                    Fill = Brushes.Transparent,
                    

                }
            };

            vCollection = new SeriesCollection
            {
                new LineSeries
                {
                    Values = vEVOList,
                    ScalesYAt = 0,
                    Stroke = new SolidColorBrush(Color.FromRgb(75, 75, 255)),
                    Fill = Brushes.Transparent,
                }
            };

            roCOllection = new SeriesCollection
            {
                new LineSeries
                {
                    Values = roEVOList,
                    ScalesYAt = 0,
                    Stroke = new SolidColorBrush(Color.FromRgb(75, 75, 255)),
                    Fill = Brushes.Transparent,
                }
            };

            pCollection = new SeriesCollection
            {
                new LineSeries
                {
                    Values = pEVOList,
                    ScalesYAt = 0,
                    Stroke = new SolidColorBrush(Color.FromRgb(75, 75, 255)),
                    Fill = Brushes.Transparent,
                }
            };

            TCollection = new SeriesCollection
            {
                new LineSeries
                {
                    Values = TEVOList,
                    ScalesYAt = 0,
                    Stroke = new SolidColorBrush(Color.FromRgb(75, 75, 255)),
                    Fill = Brushes.Transparent,
                }
            };

            MCollection = new SeriesCollection
            {
                new LineSeries
                {
                    Values = MEVOList,
                    ScalesYAt = 0,
                    Stroke = new SolidColorBrush(Color.FromRgb(75, 75, 255)),
                    Fill = Brushes.Transparent,
                }
            };
            uChart.Series = uCollection;
            vChart.Series = vCollection;
            roChart.Series = roCOllection;
            pChart.Series = pCollection;
            Tchart.Series = TCollection;
            Mchart.Series = MCollection;
        }
        private void OpenVideo_Click(object sender, RoutedEventArgs e)
        {
            string url = "https://youtu.be/Qgto2vXkQaY";
            var psi = new ProcessStartInfo();
            psi.UseShellExecute = true;
            psi.FileName = url;
            Process.Start(psi);
        }
    }
}