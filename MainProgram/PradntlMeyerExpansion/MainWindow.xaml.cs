using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
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
       
        //Gráficas
        ChartValues<double> Values = new ChartValues<double>();
        ChartValues<double> BoundaryValues = new ChartValues<double>();

        public SeriesCollection SeriesCollection;
        public MainWindow()
        {
            InitializeComponent();
            Introduction.Visibility = Visibility.Visible;

            Simulation.Children.Add(plano);

            //Establecemos los valores del gráfico //Canbiaaar
            Values.Add(1.3);
            Values.Add(3);
            BoundaryValues.Add(3);
            BoundaryValues.Add(1.3);
            setChartNumbers();
            r = new Rules(678, 0, 1.23, 0.101e6, 286.1, 2, 0.5, 1.4, 287, 10, 5.352 * Math.PI / 180, 41, 65, 40, 0.5);
            mesh = new Grid(r);
            mesh.PrandtlMeyerExpansion();
        }

        private void LoadGrid(object sender, RoutedEventArgs e)
        {
            plano.Children.Clear();

            E = r.getE()*11;
            x1 = r.getxMax() * 11;
            y1 = r.getH() * 11;
            angulo = r.getTheta();
            divisionesy = r.getJ()-1;

            //if (TabControl1.SelectedIndex == 0)
            //{
            //    E = Convert.ToInt32(E1.Text) * 11;
            //    x1 = Convert.ToInt32(x11.Text) * 11;
            //    y1 = Convert.ToInt32(y11.Text) * 11;
            //    angulo = Convert.ToDouble(theta1.Text) * Math.PI / 180;
            //    dEta= Convert.ToDouble(dEta1.Text);
            //}

            //else if (TabControl1.SelectedIndex == 1)
            //{
            //    E = Convert.ToInt32(E2.Text) * 6;
            //    x1 = Convert.ToInt32(x22.Text) * 6;
            //    y1 = Convert.ToInt32(y22.Text) * 6;
            //    angulo = Convert.ToDouble(theta2.Text);
            //    dEta = Convert.ToDouble(dEta2.Text);
            //}

            //divisionesy = 1.0 / dEta;

            //#region Frontera del PLano de SImulación
            //System.Windows.Point Point1 = new System.Windows.Point(0, 0);
            //System.Windows.Point Point2 = new System.Windows.Point(x1, 0);
            //System.Windows.Point Point3 = new System.Windows.Point(x1, y1 + (Math.Tan(angulo) * (x1 - E)));
            //System.Windows.Point Point4 = new System.Windows.Point(E, y1);
            //System.Windows.Point Point5 = new System.Windows.Point(0, y1);

            //PointCollection polygonPoints = new PointCollection();
            //polygonPoints.Add(Point1);
            //polygonPoints.Add(Point2);
            //polygonPoints.Add(Point3);
            //polygonPoints.Add(Point4);
            //polygonPoints.Add(Point5);

            //Polygon forma = new Polygon();
            //forma.Points = polygonPoints;
            //SolidColorBrush whiteBrush = new SolidColorBrush();
            //whiteBrush.Color = Colors.White;
            //forma.Stroke = whiteBrush;
            //forma.StrokeThickness = 0.1;
            //System.Windows.Thickness planoPM = new System.Windows.Thickness(475, 51, 92, 292);

            //plano.Width = x1;
            //plano.Height = y1 + (Math.Tan(angulo) * (x1 - E));
            //plano.Margin = planoPM;

            //plano.Children.Add(forma);
            //#endregion

            ////Validppppppppppppppppppp
            ////Incremento y1
            //double ay1;
            //ay1 = y1 / divisionesy;

            ////Incremento y2
            //double ay2;
            //ay2 = (y1 + Math.Tan(angulo) * (x1 - E)) / divisionesy;

            //for (int j = 0; j < divisionesy; j++)
            //{
            //    Polygon polygon = new Polygon();
            //    //Características de los rectángulos del grid
            //    System.Windows.Point Point11 = new System.Windows.Point(0, j * ay1);
            //    System.Windows.Point Point21 = new System.Windows.Point(E, j * ay1);
            //    System.Windows.Point Point31 = new System.Windows.Point(x1, j * ay2);

            //    System.Windows.Point Point41 = new System.Windows.Point(x1, (j + 1) * ay2);
            //    System.Windows.Point Point51 = new System.Windows.Point(E, (j + 1) * ay1);
            //    System.Windows.Point Point61 = new System.Windows.Point(0, (j + 1) * ay1);


            //    PointCollection polygonPoints1 = new PointCollection();
            //    polygonPoints1.Add(Point11);
            //    polygonPoints1.Add(Point21);
            //    polygonPoints1.Add(Point31);
            //    polygonPoints1.Add(Point41);
            //    polygonPoints1.Add(Point51);
            //    polygonPoints1.Add(Point61);
            //    polygon.Points = polygonPoints1;

            //    SolidColorBrush whiteBrush = new SolidColorBrush();
            //    whiteBrush.Color = Colors.White;
            //    polygon.Stroke = whiteBrush;
            //    polygon.StrokeThickness = 0.1;
            //    System.Windows.Thickness planoPM = new System.Windows.Thickness(475, 51, 92, 292);
            //    plano.Margin = planoPM;
            //    plano.Children.Add(polygon);
            //    Canvas.SetTop(polygon, 0);
            //    Canvas.SetLeft(polygon, 0);
            //}
            //LinearGradientBrush myLinearGradientBrush =new LinearGradientBrush();
            //myLinearGradientBrush.StartPoint = new Point(2, 2);
            //myLinearGradientBrush.EndPoint = new Point(2.2, 2.2);
            //myLinearGradientBrush.GradientStops.Add(new GradientStop(Colors.Yellow, 2.0));
            //myLinearGradientBrush.GradientStops.Add(new GradientStop(Colors.Red, 2.05));
            //myLinearGradientBrush.GradientStops.Add(new GradientStop(Colors.Blue, 2.1));
            //myLinearGradientBrush.GradientStops.Add(new GradientStop(Colors.LimeGreen, 2.2));
            double maxvalue=0;
            double minvalue=5;

            for (int i = 0; i < (mesh.GetXP().Count - 2); i++)
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

            for (int i = 0; i<(mesh.GetXP().Count-2); i++)
            {
                //for (int j = 0; j < divisionesy; j++)
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    if (i == 0)
                    {
                        //Grid fase
                        Polygon polygon = new Polygon();
                        //Características de los rectángulos del grid
                        System.Windows.Point Point11 = new System.Windows.Point(0, (mesh.GetYP()[0][40] - mesh.GetYP()[i][j]) * 11);
                        System.Windows.Point Point21 = new System.Windows.Point(mesh.GetXP()[i] * 11, (mesh.GetYP()[0][40] - mesh.GetYP()[i + 1][j]) * 11);
                        System.Windows.Point Point31 = new System.Windows.Point(mesh.GetXP()[i] * 11, (mesh.GetYP()[0][40] - mesh.GetYP()[i + 1][j + 1]) * 11);
                        System.Windows.Point Point41 = new System.Windows.Point(0, (mesh.GetYP()[0][40] - mesh.GetYP()[i][j + 1]) * 11);

                 
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
                        System.Windows.Thickness planoPM = new System.Windows.Thickness(500, 50, 92, 292);


                        double mediaM = (mesh.GetCell(j, i).getM()+ mesh.GetCell(j, i+1).getM()+ mesh.GetCell(j+1, i + 1).getM()+ mesh.GetCell(j+1, i).getM())/4.0;
                        byte first = 0;
                        byte second = 0;
                        byte third = 0;

                        //if (mediaM < ((maxvalue + minvalue) / 2) | mediaM == ((maxvalue + minvalue) / 2))
                        //{
                        //    first = 0;
                        //}
                        //else if (mediaM > ((maxvalue + minvalue) / 2))
                        //{
                        //    first = Convert.ToByte(Math.Round((255 * (mediaM - (maxvalue + minvalue) / 2)) / (maxvalue - (maxvalue + minvalue) / 2)));
                        //}
                        //if (mediaM < ((maxvalue + minvalue) / 2))
                        //{
                        //    second = Convert.ToByte(Math.Round((255 * (mediaM - minvalue) / 2) / (((maxvalue + minvalue) / 2) - minvalue)));
                        //}
                        //else if (mediaM > ((maxvalue + minvalue) / 2) || mediaM == ((maxvalue + minvalue) / 2))
                        //{
                        //    second = 255;
                        //}
                        //if (mediaM < ((maxvalue + minvalue) / 2))
                        //{
                        //    third = Convert.ToByte(Math.Round((255 * (mediaM - ((maxvalue + minvalue) / 2)) / (minvalue - ((maxvalue + minvalue) / 2)))));
                        //}
                        //else if (mediaM > ((maxvalue + minvalue) / 2) || mediaM == ((maxvalue + minvalue) / 2))
                        //{
                        //    third = 0;
                        //}
                        if (mediaM < ((maxvalue + minvalue) / 2))
                        {
                            first = Convert.ToByte(Math.Round(76 + (83 * (mediaM - (maxvalue + minvalue) / 2)) / (minvalue - (maxvalue + minvalue) / 2)));
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
                        polygon.Fill = mySolidColorBrush;

                        plano.Margin = planoPM;
                        plano.Children.Add(polygon);
                    }

                    else if (i!=0)
                    {
                        //Grid fase
                        Polygon polygon = new Polygon();
                        //Características de los rectángulos del grid
                        System.Windows.Point Point11 = new System.Windows.Point(mesh.GetXP()[i-1] * 11, (mesh.GetYP()[0][40] - mesh.GetYP()[i][j]) * 11);
                        System.Windows.Point Point21 = new System.Windows.Point(mesh.GetXP()[i] * 11,( mesh.GetYP()[0][40] - mesh.GetYP()[i + 1][j]) * 11);
                        System.Windows.Point Point31 = new System.Windows.Point(mesh.GetXP()[i] * 11, (mesh.GetYP()[0][40] - mesh.GetYP()[i + 1][j + 1]) * 11);
                        System.Windows.Point Point41 = new System.Windows.Point(mesh.GetXP()[i-1] * 11, (mesh.GetYP()[0][40] - mesh.GetYP()[i][j + 1]) * 11);
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
                        System.Windows.Thickness planoPM = new System.Windows.Thickness(500, 50, 92, 292);

                        double mediaM = (mesh.GetCell(j, i).getM() + mesh.GetCell(j, i + 1).getM() + mesh.GetCell(j + 1, i + 1).getM() + mesh.GetCell(j + 1, i).getM()) / 4.0;

                        byte first = 0 ;
                        byte second = 0;
                        byte third = 0;

                        //if (mediaM < ((maxvalue + minvalue) / 2) | mediaM== ((maxvalue +minvalue) / 2))
                        //{
                        //    first = 0;
                        //}
                        //else if (mediaM > ((maxvalue + minvalue) / 2))
                        //{
                        //    first = Convert.ToByte(Math.Round((255*(mediaM-(maxvalue+minvalue)/2))/(maxvalue-(maxvalue+minvalue)/2)));
                        //}
                        //if (mediaM < ((maxvalue + minvalue) / 2))
                        //{
                        //    second = Convert.ToByte(Math.Round((255 * (mediaM - minvalue) / 2) / ( ((maxvalue + minvalue) / 2) - minvalue)));
                        //}
                        //else if (mediaM > ((maxvalue + minvalue) / 2)|| mediaM == ((maxvalue + minvalue) / 2))
                        //{
                        //    second = 255;
                        //}
                        //if (mediaM < ((maxvalue + minvalue) / 2))
                        //{
                        //    third = Convert.ToByte(Math.Round((255 * (mediaM - ((maxvalue + minvalue) / 2)) / (minvalue - ((maxvalue + minvalue) / 2)))));
                        //}
                        //else if (mediaM > ((maxvalue + minvalue) / 2) || mediaM == ((maxvalue + minvalue) / 2))
                        //{
                        //    third = 255;
                        //}

                        if (mediaM < ((maxvalue + minvalue) / 2))
                        {
                            first = Convert.ToByte(Math.Round(76 + (83 * (mediaM - (maxvalue + minvalue) / 2)) / (minvalue- (maxvalue + minvalue) / 2)));
                        }
                        else if (mediaM > ((maxvalue + minvalue) / 2))
                        {
                            first = Convert.ToByte(Math.Round(76 + (164 * (mediaM - (maxvalue + minvalue) / 2)) / (maxvalue - (maxvalue + minvalue) / 2)));
                        }
                        else if (mediaM == ((maxvalue + minvalue) / 2))
                        {
                            first= Convert.ToByte(76);
                        }
                        second = Convert.ToByte(Math.Round(95 + (160 * (mediaM - minvalue) ) / (maxvalue - minvalue)));
                        third = Convert.ToByte(Math.Round((160 * (mediaM - maxvalue) / (minvalue - maxvalue))));

                        SolidColorBrush mySolidColorBrush = new SolidColorBrush();
                        mySolidColorBrush.Color = Color.FromRgb(first,second,third);
                        polygon.Fill = mySolidColorBrush;

                        plano.Margin = planoPM;
                        plano.Children.Add(polygon);
                    }
                }
            }
        }


        private void insideGrid()
        {

            //for (int i = 0; i < mesh.GetXP().Count; i++)
            //{
            //    for (int j = 0; j < divisionesy; j++)
            //    {
            //        //Grid fase
            //        Polygon polygon = new Polygon();
            //        //Características de los rectángulos del grid
            //        System.Windows.Point Point11 = new System.Windows.Point(mesh.GetXP()[i], mesh.GetYP()[i][j]);
            //        System.Windows.Point Point21 = new System.Windows.Point(mesh.GetXP()[i + 1], mesh.GetYP()[i+1][j]);
            //        System.Windows.Point Point31 = new System.Windows.Point(mesh.GetXP()[i + 1], mesh.GetYP()[i + 1][j+1]);
            //        System.Windows.Point Point41 = new System.Windows.Point(mesh.GetXP()[i], mesh.GetYP()[i][j + 1]);
            //        PointCollection polygonPoints1 = new PointCollection();
            //        polygonPoints1.Add(Point11);
            //        polygonPoints1.Add(Point21);
            //        polygonPoints1.Add(Point31);
            //        polygonPoints1.Add(Point41);
            //        polygon.Points = polygonPoints1;

            //        SolidColorBrush whiteBrush = new SolidColorBrush();
            //        whiteBrush.Color = Colors.White;
            //        polygon.Stroke = whiteBrush;
            //        polygon.StrokeThickness = 0.1;
            //        plano.Children.Add(polygon);
            //        Canvas.SetBottom(polygon, 0);
            //        Canvas.SetLeft(polygon, E - 10);
            //    }
            //}



            //#region Interior del PLano de SImulación
            ////Rectangulo
            //for (int i = 0; i < 10; i++)
            //{
            //    for (int j = 0; j < divisionesy + 1; j++)
            //    {
            //        //Grid fase
            //        Rectangle pRectangle = new Rectangle();
            //        //Características de los rectángulos del grid
            //        pRectangle.Width = E / 10;
            //        pRectangle.Height = y1 / divisionesy;
            //        pRectangle.Fill = new SolidColorBrush(Colors.Blue);
            //        pRectangle.Fill.Opacity = 0;
            //        pRectangle.StrokeThickness = 0.1;
            //        pRectangle.Stroke = Brushes.White;

            //        plano.Children.Add(pRectangle);

            //        //Posición del rectángulo dentro del grid
            //        Canvas.SetLeft(pRectangle, i * pRectangle.Width);
            //        Canvas.SetTop(pRectangle, j * pRectangle.Height);
            //    }
            //}

            ////Poligono
            //for (int i = 0; i < 10; i++)
            //{
            //    for (int j = 0; j < divisionesy; j++)
            //    {
            //        //Incremento x
            //        double ax;
            //        ax = (x1 - E) / 10;

            //        double ay1;
            //        if (i == 0)
            //        {
            //            //Incremento y1
            //            ay1 = y1 / divisionesy;
            //        }

            //        else
            //        {
            //            //Incremento y1
            //            ay1 = (y1 + Math.Tan(angulo) * (ax * (i))) / divisionesy;
            //        }

            //        //Incremento y2
            //        double ay2;
            //        ay2 = (y1 + Math.Tan(angulo) * (ax * (i + 1))) / divisionesy;

            //        //Grid fase
            //        Polygon polygon = new Polygon();
            //        //Características de los rectángulos del grid
            //        System.Windows.Point Point11 = new System.Windows.Point(10 + i * ax, j * ay1);
            //        System.Windows.Point Point21 = new System.Windows.Point(10 + (i + 1) * ax, j * ay2);
            //        System.Windows.Point Point31 = new System.Windows.Point(10 + (i + 1) * ax, (j + 1) * ay2);
            //        System.Windows.Point Point41 = new System.Windows.Point(10 + i * ax, (j + 1) * ay1);
            //        PointCollection polygonPoints1 = new PointCollection();
            //        polygonPoints1.Add(Point11);
            //        polygonPoints1.Add(Point21);
            //        polygonPoints1.Add(Point31);
            //        polygonPoints1.Add(Point41);
            //        polygon.Points = polygonPoints1;

            //        SolidColorBrush whiteBrush = new SolidColorBrush();
            //        whiteBrush.Color = Colors.White;
            //        polygon.Stroke = whiteBrush;
            //        polygon.StrokeThickness = 0.1;
            //        plano.Children.Add(polygon);
            //        Canvas.SetTop(polygon, 0);
            //        Canvas.SetLeft(polygon, E - 10);
            //    }
            //}
            //#endregion
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

        private void lntro_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Visible;
            Simulation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
        }

        private void Simualtion_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Visible;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
        }

        private void vTutorial_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Visible;
            AboutUs.Visibility = Visibility.Hidden;
        }

        private void aboutUs_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Visible;
        }
    }
}
