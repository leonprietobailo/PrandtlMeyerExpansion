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

namespace PradntlMeyerExpansion
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        Grid mesh;
        Rules r;

        //valores iniciales
        double angulo, E,x1,y1;
        double dEta, dXi;
        Canvas plano = new Canvas();
        int divisionesy;
        Polygon[,] Polygons;

        //Gráficas
        ChartValues<double> Values = new ChartValues<double>();
        ChartValues<double> BoundaryValues = new ChartValues<double>();

        //Colores
        SolidColorBrush color = new SolidColorBrush();

        public SeriesCollection SeriesCollection;

        DataTable UTable=new DataTable();
        DataTable AndersonTable = new DataTable();
        public MainWindow()
        {
            InitializeComponent();

            color.Color = Color.FromRgb(82, 222, 197);

            Introduction.Visibility = Visibility.Visible;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
            EstudioAvanzado.Visibility= Visibility.Hidden;

            Info.Visibility = Visibility.Hidden;

            topic.BorderBrush = color;
            topic.Foreground = color;
            simulation.BorderBrush = Brushes.White;
            simulation.Foreground = Brushes.White;
            ValidationButton.BorderBrush = Brushes.White;
            ValidationButton.Foreground = Brushes.White;
            vTutorial.BorderBrush = Brushes.White;
            vTutorial.Foreground = Brushes.White;
            aboutUs.BorderBrush = Brushes.White;
            aboutUs.Foreground = Brushes.White;
            EstudioAvanzada.BorderBrush = Brushes.White;
            EstudioAvanzada.Foreground = Brushes.White;
            EstudioAvanzada.BorderBrush = Brushes.White;
            EstudioAvanzada.Foreground = Brushes.White;

            Simulation.Children.Add(plano);

            //Establecemos los valores del gráfico //Canbiaaar
            Values.Add(1.3);
            Values.Add(3);
            BoundaryValues.Add(3);
            BoundaryValues.Add(1.3);
            setChartNumbers();
            //r = new Rules(678, 0, 1.23, 0.101e6, 286.1, 2, 0.5, 1.4, 287, 10, 5.352 * Math.PI / 180, 41, 65, 40, 0.5);
        }

        private void lntro_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Visible;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
            EstudioAvanzado.Visibility = Visibility.Hidden;

            topic.BorderBrush = color;
            topic.Foreground = color;
            simulation.BorderBrush = Brushes.White;
            simulation.Foreground = Brushes.White;
            ValidationButton.BorderBrush = Brushes.White;
            ValidationButton.Foreground = Brushes.White;
            vTutorial.BorderBrush = Brushes.White;
            vTutorial.Foreground = Brushes.White;
            aboutUs.BorderBrush = Brushes.White;
            aboutUs.Foreground = Brushes.White;
            EstudioAvanzada.BorderBrush = Brushes.White;
            EstudioAvanzada.Foreground = Brushes.White;
        }

        private void Simualtion_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Visible;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
            EstudioAvanzado.Visibility = Visibility.Hidden;

            topic.BorderBrush = Brushes.White;
            topic.Foreground = Brushes.White;
            simulation.BorderBrush = color;
            simulation.Foreground = color;
            ValidationButton.BorderBrush = Brushes.White;
            ValidationButton.Foreground = Brushes.White;
            vTutorial.BorderBrush = Brushes.White;
            vTutorial.Foreground = Brushes.White;
            aboutUs.BorderBrush = Brushes.White;
            aboutUs.Foreground = Brushes.White;
            EstudioAvanzada.BorderBrush = Brushes.White;
            EstudioAvanzada.Foreground = Brushes.White;
        }

        private void Validation_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Visible;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;

            topic.BorderBrush = Brushes.White;
            topic.Foreground = Brushes.White;
            simulation.BorderBrush = Brushes.White;
            simulation.Foreground = Brushes.White;
            ValidationButton.BorderBrush = color;
            ValidationButton.Foreground = color;
            vTutorial.BorderBrush = Brushes.White;
            vTutorial.Foreground = Brushes.White;
            aboutUs.BorderBrush = Brushes.White;
            aboutUs.Foreground = Brushes.White;
            EstudioAvanzada.BorderBrush = Brushes.White;
            EstudioAvanzada.Foreground = Brushes.White;
            EstudioAvanzado.Visibility = Visibility.Hidden;

            Utable.IsChecked = true;
           
        }

        private void vTutorial_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Visible;
            AboutUs.Visibility = Visibility.Hidden;
            EstudioAvanzado.Visibility = Visibility.Hidden;

            topic.BorderBrush = Brushes.White;
            topic.Foreground = Brushes.White;
            simulation.BorderBrush = Brushes.White;
            simulation.Foreground = Brushes.White;
            ValidationButton.BorderBrush = Brushes.White;
            ValidationButton.Foreground = Brushes.White;
            vTutorial.BorderBrush = color;
            vTutorial.Foreground = color;
            aboutUs.BorderBrush = Brushes.White;
            aboutUs.Foreground = Brushes.White;
            EstudioAvanzada.BorderBrush = Brushes.White;
            EstudioAvanzada.Foreground = Brushes.White;
        }

        private void aboutUs_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Visible;
            EstudioAvanzado.Visibility = Visibility.Hidden;

            topic.BorderBrush = Brushes.White;
            topic.Foreground = Brushes.White;
            simulation.BorderBrush = Brushes.White;
            simulation.Foreground = Brushes.White;
            ValidationButton.BorderBrush = Brushes.White;
            ValidationButton.Foreground = Brushes.White;
            vTutorial.BorderBrush = Brushes.White;
            vTutorial.Foreground = Brushes.White;
            aboutUs.BorderBrush = color;
            aboutUs.Foreground = color;
            EstudioAvanzada.BorderBrush = Brushes.White;
            EstudioAvanzada.Foreground = Brushes.White;
        }

        private void estudio_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
            EstudioAvanzado.Visibility = Visibility.Visible;

            topic.BorderBrush = Brushes.White;
            topic.Foreground = Brushes.White;
            simulation.BorderBrush = Brushes.White;
            simulation.Foreground = Brushes.White;
            ValidationButton.BorderBrush = Brushes.White;
            ValidationButton.Foreground = Brushes.White;
            vTutorial.BorderBrush = Brushes.White;
            vTutorial.Foreground = Brushes.White;
            aboutUs.BorderBrush = Brushes.White;
            aboutUs.Foreground = Brushes.White;
            EstudioAvanzada.BorderBrush = color;
            EstudioAvanzada.Foreground = color;
        }
        private void LoadGrid(object sender, RoutedEventArgs e)
        {
            plano.Children.Clear();

            double u = Convert.ToDouble(u1.Text);
            double v = Convert.ToDouble(v1.Text);
            double rho = Convert.ToDouble(rho1.Text);
            double p = Convert.ToDouble(p1.Text);
            double T = Convert.ToDouble(T1.Text);
            double M = Convert.ToDouble(M1.Text);
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
            r = new Rules(u, v, rho, p, T, M, Cy, gamma, R, E, theta * Math.PI / 180, j1, xmax, H, C);

            //Cambiar cuando hagamos el ajuste al mostrar
            E = r.getE() * 11;
            x1 = r.getxMax() * 11;
            y1 = r.getH() * 11;
            angulo = r.getTheta();
            divisionesy = r.getJ() - 1;

            //Incremento y1
            double ay1;
            ay1 = y1 / divisionesy;

            //Incremento y2
            double ay2;
            ay2 = (y1 + Math.Tan(angulo) * (x1 - E)) / divisionesy;

            for (int j = 0; j < divisionesy; j++)
            {
                Polygon polygon = new Polygon();
                //Características de los rectángulos del grid
                System.Windows.Point Point11 = new System.Windows.Point(0, j * ay1);
                System.Windows.Point Point21 = new System.Windows.Point(E, j * ay1);
                System.Windows.Point Point31 = new System.Windows.Point(x1, j * ay2);

                System.Windows.Point Point41 = new System.Windows.Point(x1, (j + 1) * ay2);
                System.Windows.Point Point51 = new System.Windows.Point(E, (j + 1) * ay1);
                System.Windows.Point Point61 = new System.Windows.Point(0, (j + 1) * ay1);


                PointCollection polygonPoints1 = new PointCollection();
                polygonPoints1.Add(Point11);
                polygonPoints1.Add(Point21);
                polygonPoints1.Add(Point31);
                polygonPoints1.Add(Point41);
                polygonPoints1.Add(Point51);
                polygonPoints1.Add(Point61);
                polygon.Points = polygonPoints1;

                SolidColorBrush whiteBrush = new SolidColorBrush();
                whiteBrush.Color = Colors.White;
                polygon.Stroke = whiteBrush;
                polygon.StrokeThickness = 0.1;
                System.Windows.Thickness planoPM = new System.Windows.Thickness(520, 51, 92, 292);
                plano.Margin = planoPM;
                plano.Children.Add(polygon);
                Canvas.SetTop(polygon, 0);
                Canvas.SetLeft(polygon, 0);
            }
        }
        private void Run_Click(object sender, RoutedEventArgs e)
        {
            plano.Children.Clear();
            mesh = new Grid(r);
            mesh.PrandtlMeyerExpansion();

            E = r.getE() * 11;
            x1 = r.getxMax() * 11;
            y1 = r.getH() * 11;
            angulo = r.getTheta();
            divisionesy = r.getJ() - 1;

            Polygons = new Polygon[mesh.GetXP().Count - 1, divisionesy];

            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                //for (int j = 0; j < divisionesy; j++)
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    //Grid fase
                    Polygon polygon = new Polygon();
                    //Características de los rectángulos del grid
                    System.Windows.Point Point11 = new System.Windows.Point(mesh.GetXP()[i] * 11, (mesh.GetYP()[0][40] - mesh.GetYP()[i][j]) * 11);
                    System.Windows.Point Point21 = new System.Windows.Point(mesh.GetXP()[i + 1] * 11, (mesh.GetYP()[0][40] - mesh.GetYP()[i + 1][j]) * 11);
                    System.Windows.Point Point31 = new System.Windows.Point(mesh.GetXP()[i + 1] * 11, (mesh.GetYP()[0][40] - mesh.GetYP()[i + 1][j + 1]) * 11);
                    System.Windows.Point Point41 = new System.Windows.Point(mesh.GetXP()[i] * 11, (mesh.GetYP()[0][40] - mesh.GetYP()[i][j + 1]) * 11);
                    PointCollection polygonPoints1 = new PointCollection();
                    polygonPoints1.Add(Point11);
                    polygonPoints1.Add(Point21);
                    polygonPoints1.Add(Point31);
                    polygonPoints1.Add(Point41);
                    polygon.Points = polygonPoints1;

                    SolidColorBrush whiteBrush = new SolidColorBrush();
                    whiteBrush.Color = Colors.White;
                    polygon.Stroke = whiteBrush;
                    polygon.StrokeThickness = 0.1;
                    System.Windows.Thickness planoPM = new System.Windows.Thickness(520, 51, 92, 292);

                    plano.Margin = planoPM;
                    plano.Children.Add(polygon);

                    polygon.Tag = new Point(i, j);
                    polygon.MouseEnter += new MouseEventHandler(polygon_MouseEnter);
                    polygon.MouseLeave += new MouseEventHandler(polygon_MouseLeave);
                    Polygons[i, j] = polygon;
                }
            }
            uRB.IsChecked = true;
        }

        private void CreateUTable()
        {
            UTable.Columns.Add("");
            for (int j = 0; j < divisionesy; j++)
            {
                DataRow row = UTable.NewRow();
                row[0] = j+1;
                UTable.Rows.Add(row);
            }
        }

        private void CreateUTableAnderson()
        {
            List<double> columna19 = new List<double> {678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678, 678,678, 678, 679,683,691,701,707};
            List<double> columna88 = new List<double> { 678, 678, 678, 679, 679, 680, 681, 683, 685, 688, 690, 693, 696, 699, 702, 705, 707, 709, 711, 713, 713, 713, 712, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 711, 710,705 };
            AndersonTable.Columns.Add("");
            for (int j = 0; j < divisionesy; j++)
            {
                DataRow row = AndersonTable.NewRow();
                row[0] = j+1;
                AndersonTable.Rows.Add(row);
            }
            AndersonTable.Columns.Add("19");
            AndersonTable.Columns.Add("88");
            for (int j = 0; j < divisionesy; j++)
            {
                DataRow row = AndersonTable.Rows[j];
                row["19"] = columna88[j];
                row["88"] = columna88[j];
            }
        }

        private void u_Checked(object sender, RoutedEventArgs e)
        {
            CreateUTable();
            CreateUTableAnderson();
            double maxvalue = -100000000;
            double minvalue = 100000000;

            List<double> test = new List<double>();
            test = mesh.GetXP();


            double uValue;
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                for (int j = 0; j < divisionesy; j++)
                {
                    uValue = mesh.GetCell(j, i).getU();
                    if (uValue > maxvalue)
                    {
                        maxvalue = uValue;
                    }

                    if (uValue < minvalue)
                    {
                        minvalue = uValue;
                    }
                }
            }
            int y = 0;
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                y = 0;
                UTable.Columns.Add(Convert.ToString(i+1));
                //for (int j = 0; j < divisionesy; j++)
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    double mediaM = (mesh.GetCell(j, i).getU() + mesh.GetCell(j, i + 1).getU() + mesh.GetCell(j + 1, i + 1).getU() + mesh.GetCell(j + 1, i).getU()) / 4.0;

                    byte first = 0;
                    byte second = 0;
                    byte third = 0;

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
                    second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));

                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    Polygons[i,j].Fill = mySolidColorBrush;

                    DataRow row = UTable.Rows[y];
                    row[Convert.ToString(i+1)] = Math.Round(mesh.GetCell(j, i).getU(), 4);
                    y++;
                }
            }
        }
        private void uTable_Checked(object sender, RoutedEventArgs e)
        {
            gridData.DataContext = UTable.DefaultView;
            gridAndersonData.DataContext = AndersonTable.DefaultView;
        }


        private void v_Checked(object sender, RoutedEventArgs e)
        {
            double maxvalue = -100000000;
            double minvalue = 100000000;

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

            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                //for (int j = 0; j < divisionesy; j++)
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    double mediaM = (mesh.GetCell(j, i).getV() + mesh.GetCell(j, i + 1).getV() + mesh.GetCell(j + 1, i + 1).getV() + mesh.GetCell(j + 1, i).getV()) / 4.0;

                    byte first = 0;
                    byte second = 0;
                    byte third = 0;

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
                    second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));

                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    Polygons[i, j].Fill = mySolidColorBrush;
                }
            }
        }

        private void rho_Checked(object sender, RoutedEventArgs e)
        {
            double maxvalue = -100000000;
            double minvalue = 100000000;

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

            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                //for (int j = 0; j < divisionesy; j++)
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    double mediaM = (mesh.GetCell(j, i).getRO() + mesh.GetCell(j, i + 1).getRO() + mesh.GetCell(j + 1, i + 1).getRO() + mesh.GetCell(j + 1, i).getRO()) / 4.0;

                    byte first = 0;
                    byte second = 0;
                    byte third = 0;

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
                    second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));

                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    Polygons[i,j].Fill = mySolidColorBrush;
                }
            }
        }


        private void p_Checked(object sender, RoutedEventArgs e)
        {
            double maxvalue = -100000000;
            double minvalue = 100000000;

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

            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                //for (int j = 0; j < divisionesy; j++)
                for (int j = (divisionesy - 1); j > -1; j--)
                {

                    double mediaM = (mesh.GetCell(j, i).getP() + mesh.GetCell(j, i + 1).getP() + mesh.GetCell(j + 1, i + 1).getP() + mesh.GetCell(j + 1, i).getP()) / 4.0;

                    byte first = 0;
                    byte second = 0;
                    byte third = 0;

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
                    second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));

                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    Polygons[i,j].Fill = mySolidColorBrush;
                }
            }
        }

        private void T_Checked(object sender, RoutedEventArgs e)
        {
            double maxvalue = -100000000;
            double minvalue = 100000000;

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

            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                //for (int j = 0; j < divisionesy; j++)
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    double mediaM = (mesh.GetCell(j, i).getT() + mesh.GetCell(j, i + 1).getT() + mesh.GetCell(j + 1, i + 1).getT() + mesh.GetCell(j + 1, i).getT()) / 4.0;

                    byte first = 0;
                    byte second = 0;
                    byte third = 0;

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
                    second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));

                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    Polygons[i, j].Fill = mySolidColorBrush;
                }
            }
        }

        private void M_Checked(object sender, RoutedEventArgs e)
        {
            double maxvalue = -100000000;
            double minvalue = 100000000;

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

            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                //for (int j = 0; j < divisionesy; j++)
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    double mediaM = (mesh.GetCell(j, i).getM() + mesh.GetCell(j, i + 1).getM() + mesh.GetCell(j + 1, i + 1).getM() + mesh.GetCell(j + 1, i).getM()) / 4.0;

                    byte first = 0;
                    byte second = 0;
                    byte third = 0;

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
                    second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue)) / (maxvalue - minvalue)));
                    third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));

                    SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                    mySolidColorBrush.Color = Color.FromRgb(first, second, third);
                    Polygons[i, j].Fill = mySolidColorBrush;
                }
            }
        }



        private void vTable_Checked(object sender, RoutedEventArgs e)
        {

        }

        private void rhoTable_Checked(object sender, RoutedEventArgs e)
        {

        }

        private void pTable_Checked(object sender, RoutedEventArgs e)
        {

        }

        private void TTable_Checked(object sender, RoutedEventArgs e)
        {

        }

        private void MTable_Checked(object sender, RoutedEventArgs e)
        {

        }
        private void polygon_MouseEnter(object sender, MouseEventArgs e)
        {
            //Obtnemos la ubicación del puntero
            Polygon pol = (Polygon)sender;
            Point p = (Point)pol.Tag;
            //Mostramos las etiquetas correspondientes a los estados de las celdas
            Info.Visibility = Visibility.Visible;
            //Actualizamos los valores de fase y temperatura
            uData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getU(), 4);
            vData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getV(), 4);
            rhoData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getRO(), 4);
            pData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getP(), 4);
            TData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getT(), 4);
            MData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getM(), 4);
            xData.Content = p.X+1;
            yData.Content = p.Y+1;
        }

        private void polygon_MouseLeave(object sender, MouseEventArgs e)
        {
            Info.Visibility = Visibility.Hidden;
        }

        private void setChartNumbers()
        {
            SeriesCollection = new SeriesCollection
            {
                //Gráfica de fase, con sus propiedades
                new LineSeries
                {
                    Values = Values,
                    ScalesYAt = 0,
                    Stroke = new SolidColorBrush(Color.FromRgb(75, 75, 255)),
                    Fill = Brushes.Transparent,
                    Title = "Inside"
                },

                new LineSeries
                {
                    Values = BoundaryValues,
                    ScalesYAt = 0,
                    Stroke = new SolidColorBrush(Color.FromRgb(255, 75, 75)),
                    Fill = Brushes.Transparent,
                    Title = "Boundary"
                }
            };
            Chart1.Series = SeriesCollection;
        }
    }
}
