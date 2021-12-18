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
using System.IO;
using System.Diagnostics;

namespace PradntlMeyerExpansion
{
    public partial class MainWindow : Window
    {
        Grid mesh,meshCte;
        Rules r,rCte;

        //valores iniciales
        double angulo, E,x1,y1,dEta,dXi;
        Canvas plano = new Canvas();
        int divisionesy;
        Polygon[,] Polygons;

        // Colores
        SolidColorBrush mySolidColorBrush = new SolidColorBrush();
        SolidColorBrush color = new SolidColorBrush();

        double umax, umin, vmax, vmin, rhomax, rhomin, pmax, pmin, Tmax, Tmin, Mmax, Mmin;

        //Gráficas
        //ChartValues<double> Values = new ChartValues<double>();
        //ChartValues<double> BoundaryValues = new ChartValues<double>();

        //Estudio avanzado
        ChartValues<ObservablePoint> uEVOList = new ChartValues<ObservablePoint>();
        ChartValues<ObservablePoint> vEVOList = new ChartValues<ObservablePoint>();
        ChartValues<ObservablePoint> roEVOList = new ChartValues<ObservablePoint>();
        ChartValues<ObservablePoint> pEVOList = new ChartValues<ObservablePoint>();
        ChartValues<ObservablePoint> TEVOList = new ChartValues<ObservablePoint>();
        ChartValues<ObservablePoint> MEVOList = new ChartValues<ObservablePoint>();

        public SeriesCollection SeriesCollection, uCollection, vCollection, roCOllection, pCollection, TCollection, MCollection;

        //Tables
        DataTable UTable = new DataTable();
        DataTable VTable = new DataTable();
        DataTable rhoTable = new DataTable();
        DataTable pTable = new DataTable();
        DataTable TTable = new DataTable();
        DataTable MTable = new DataTable();

        DataTable AndersonUTable = new DataTable();
        DataTable AndersonVTable = new DataTable();
        DataTable AndersonRhoTable = new DataTable();
        DataTable AndersonpTable = new DataTable();
        DataTable AndersonTTable = new DataTable();
        DataTable AndersonMTable = new DataTable();

        public MainWindow()
        {
            InitializeComponent();

            mySolidColorBrush.Color = Color.FromRgb(75, 70, 70);

            color.Color = Color.FromRgb(82, 222, 197);

            Introduction.Visibility = Visibility.Visible;
            Introduction.Background = mySolidColorBrush;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
            EstudioAvanzado.Visibility= Visibility.Hidden;

            Info.Visibility = Visibility.Hidden;

            topic.Background = mySolidColorBrush;
            simulation.Background = Brushes.Transparent;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = Brushes.Transparent;

            EstudioAvanzada.IsEnabled = false;

            Simulation.Children.Add(plano);

            //Establecemos los valores del gráfico //Canbiaaar
            //Values.Add(1.3);
            //Values.Add(3);
            //BoundaryValues.Add(3);
            //BoundaryValues.Add(1.3);
            //setChartNumbers();

            //Sirve para cargar la tabla
            rCte = new Rules(678, 0, 1.23, 0.101e6, 286.1, 0.5, 1.4, 287, 10, 5.352 * Math.PI / 180, 41, 65, 40, 0.5);
            meshCte = new Grid(rCte);
            meshCte.PrandtlMeyerExpansion();


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

            double relU = Math.Round((uEVO - uReal) / uReal * 100, 4);
            double relV = Math.Round((vEVO - vReal) / vReal * 100,4);
            double relRO = Math.Round((roEVO - roReal) / roReal * 100,4);
            double relP = Math.Round((pEVO - pReal) / pReal * 100,4);
            double relT = Math.Round((TEVO - TReal) / TReal * 100,4);
            double relM = Math.Round((MEVO - MReal) / MReal * 100,4);

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


            gridData.DataContext = UTable.DefaultView;

        }
        private void CreateTables()
        {
            int divisionesycte = rCte.getJ();

            // u Table
            UTable.Columns.Add("y-x");
            // u Anderson Table
            AndersonUTable.Columns.Add("y-x");
            AndersonUTable.Columns.Add("18");
            AndersonUTable.Columns.Add("89");
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
                DataRow uRow = UTable.NewRow();
                uRow["y-x"] = Convert.ToString(j + 1);
                UTable.Rows.Add(uRow);
                // u Anderson Table
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

        private void uTable_Checked(object sender, RoutedEventArgs e)
        {
            gridData.Visibility = Visibility.Visible;
            gridData.DataContext = UTable.DefaultView;
            gridAndersonData.DataContext = AndersonUTable.DefaultView;
        }

        private void vTable_Checked(object sender, RoutedEventArgs e)
        {
            gridData.DataContext = VTable.DefaultView;
            gridAndersonData.DataContext = AndersonVTable.DefaultView;
        }

        private void rhoTable_Checked(object sender, RoutedEventArgs e)
        {
            gridData.DataContext = rhoTable.DefaultView;
            gridAndersonData.DataContext = AndersonRhoTable.DefaultView;
        }

        private void pTable_Checked(object sender, RoutedEventArgs e)
        {
            gridData.DataContext = pTable.DefaultView;
            gridAndersonData.DataContext = AndersonpTable.DefaultView;
        }

        private void TTable_Checked(object sender, RoutedEventArgs e)
        {
            gridData.DataContext = TTable.DefaultView;
            gridAndersonData.DataContext = AndersonTTable.DefaultView;
        }

        private void Default_Click_1(object sender, RoutedEventArgs e)
        {
            E1.Text ="10";
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

        private void MTable_Checked(object sender, RoutedEventArgs e)
        {
            gridData.DataContext = MTable.DefaultView;
            gridAndersonData.DataContext = AndersonMTable.DefaultView;
        }
        

        private void lntro_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Visible;
            Introduction.Background = mySolidColorBrush;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;

            topic.Background = mySolidColorBrush;
            simulation.Background = Brushes.Transparent;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = Brushes.Transparent;
        }

        private void Simualtion_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Visible;
            Simulation.Background = mySolidColorBrush;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
            Grids.Visibility = Visibility.Hidden;

            topic.Background = Brushes.Transparent;
            simulation.Background = mySolidColorBrush;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = Brushes.Transparent;

            EstudioAvanzada.IsEnabled = true;
            Grid.IsEnabled = false;

            
        }

        private void Validation_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Visible;
            Validation.Background = mySolidColorBrush;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;

            topic.Background = Brushes.Transparent;
            simulation.Background = Brushes.Transparent;
            ValidationButton.Background = mySolidColorBrush;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = Brushes.Transparent;

            Utable.IsChecked = true;
        }

        private void vTutorial_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Visible;
            VideoTutorial.Background = mySolidColorBrush;
            AboutUs.Visibility = Visibility.Hidden;

            topic.Background = Brushes.Transparent;
            simulation.Background = Brushes.Transparent;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = mySolidColorBrush;
            aboutUs.Background = Brushes.Transparent;
        }

        private void aboutUs_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Hidden;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Visible;
            AboutUs.Background = mySolidColorBrush;

            topic.Background = Brushes.Transparent;
            simulation.Background = Brushes.Transparent;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = mySolidColorBrush;
        }

        private void estudio_Click(object sender, RoutedEventArgs e)
        {
            Introduction.Visibility = Visibility.Hidden;
            Simulation.Visibility = Visibility.Visible;
            Simulation.Background = mySolidColorBrush;
            EstudioAvanzado.Visibility = Visibility.Visible;
            Validation.Visibility = Visibility.Hidden;
            VideoTutorial.Visibility = Visibility.Hidden;
            AboutUs.Visibility = Visibility.Hidden;
            Grids.Visibility= Visibility.Hidden;
            plano.Visibility = Visibility.Hidden;

            topic.Background = Brushes.Transparent;
            simulation.Background = mySolidColorBrush;
            ValidationButton.Background = Brushes.Transparent;
            vTutorial.Background = Brushes.Transparent;
            aboutUs.Background = Brushes.Transparent;

            EstudioAvanzada.IsEnabled = false;
            Grid.IsEnabled = true;

            computeEvolutionChange();
            setChartNumbers();
        }

        private void Grid_Click(object sender, RoutedEventArgs e)
        {
            Grids.Visibility = Visibility.Visible;
            plano.Visibility = Visibility.Visible;
            EstudioAvanzado.Visibility = Visibility.Hidden;

            EstudioAvanzada.IsEnabled = true;
            Grid.IsEnabled = false;
        }
        private void LoadGrid(object sender, RoutedEventArgs e)
        {
            RunSim.IsEnabled = true;
            Grid.IsEnabled = false;
            EstudioAvanzada.IsEnabled = false;
            Grids.Visibility = Visibility.Hidden;
            plano.Visibility = Visibility.Visible;
            EstudioAvanzado.Visibility = Visibility.Hidden;


            plano.Children.Clear();

            double u = Convert.ToDouble(u1.Text);
            double v = Convert.ToDouble(v1.Text);
            double rho = Convert.ToDouble(rho1.Text);
            double p = Convert.ToDouble(p1.Text);
            double T = Convert.ToDouble(T1.Text);
            //double M = Convert.ToDouble(M1.Text);
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

            //Cambiar cuando hagamos el ajuste al mostrar
            E = r.getE();
            x1 = r.getxMax();
            y1 = r.getH();
            angulo = r.getTheta();
            divisionesy = r.getJ() - 1;

            angulo = r.getTheta();
            divisionesy = r.getJ() - 1;

            //Incremento y1
            double ay1;
            ay1 = y1 / divisionesy;

            //Incremento y2
            double ay2;
            ay2 = (y1 + Math.Tan(angulo) * (x1 - E)) / divisionesy;

            double MaxMax1 = y1+Math.Tan(angulo)*(x1-E);

            for (int j = 0; j < divisionesy; j++)
            {

                Polygon polygon = new Polygon();
                //Características de los rectángulos del grid
                System.Windows.Point Point11 = new System.Windows.Point(0, j * ay1 * 11 * 45.1028 / MaxMax1);
                System.Windows.Point Point21 = new System.Windows.Point(E * 715 / x1, j * ay1 * 11 * 45.1028 / MaxMax1);
                System.Windows.Point Point31 = new System.Windows.Point(x1 * 715 / x1, j * ay2 * 11 * 45.1028 / MaxMax1);

                System.Windows.Point Point41 = new System.Windows.Point(x1 * 715 / x1, (j + 1) * ay2 * 11 * 45.1028 / MaxMax1);
                System.Windows.Point Point51 = new System.Windows.Point(E * 715 / x1, (j + 1) * ay1 * 11 * 45.1028 / MaxMax1);
                System.Windows.Point Point61 = new System.Windows.Point(0, (j + 1) * ay1 * 11 * 45.1028 / MaxMax1);


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
                System.Windows.Thickness planoPM = new System.Windows.Thickness(520, 100, 92, 292);
                plano.Margin = planoPM;
                plano.Children.Add(polygon);
                Canvas.SetTop(polygon, 0);
                Canvas.SetLeft(polygon, 0);
            }
        }
        private void Run_Click(object sender, RoutedEventArgs e)
        { 
            Grid.IsEnabled = false;
            EstudioAvanzada.IsEnabled = true;
            Grids.Visibility = Visibility.Visible;

            plano.Children.Clear();
            mesh = new Grid(r);
            mesh.PrandtlMeyerExpansion();
            mesh.data();


            Polygons = new Polygon[mesh.GetXP().Count - 1, divisionesy];

            int Max =mesh.GetXP().Count;
            double maxmax = mesh.GetXP()[Max-1];

            int MaxY = mesh.GetYP()[0].Length;
            double MaxMax= (mesh.GetYP()[Max-1][MaxY-1]- mesh.GetYP()[Max - 1][0]);

            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
                //for (int j = 0; j < divisionesy; j++)
                for (int j = (divisionesy - 1); j > -1; j--)
                {
                    //Grid fase
                    Polygon polygon = new Polygon();
                    //Características de los rectángulos del grid
                    System.Windows.Point Point11 = new System.Windows.Point(mesh.GetXP()[i] * 715/ maxmax, (mesh.GetYP()[0][40] - mesh.GetYP()[i][j]) * 11 * 45.1028/MaxMax);
                    System.Windows.Point Point21 = new System.Windows.Point(mesh.GetXP()[i + 1] * 715 / maxmax, (mesh.GetYP()[0][40] - mesh.GetYP()[i + 1][j]) * 11 * 45.1028 / MaxMax);
                    System.Windows.Point Point31 = new System.Windows.Point(mesh.GetXP()[i + 1] * 715 / maxmax, (mesh.GetYP()[0][40] - mesh.GetYP()[i + 1][j + 1]) * 11 * 45.1028 / MaxMax);
                    System.Windows.Point Point41 = new System.Windows.Point(mesh.GetXP()[i] * 715 / maxmax, (mesh.GetYP()[0][40] - mesh.GetYP()[i][j + 1]) * 11 * 45.1028 / MaxMax);
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
                    System.Windows.Thickness planoPM = new System.Windows.Thickness(520, 100, 92, 292);

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



        private void u_Checked(object sender, RoutedEventArgs e)
        {
            //CreateUTable();
            //CreateUTableAnderson();
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
            umax = maxvalue;
            umin = minvalue;
            int y = 0;
            for (int i = 0; i < (mesh.GetXP().Count - 1); i++)
            {
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


                }
            }
            LinearGradientBrush myHorizontalGradient = new LinearGradientBrush();
            myHorizontalGradient.StartPoint = new Point(0, 0.3);
            myHorizontalGradient.EndPoint = new Point(1, 0.3);
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(240, 255, 0), 1));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(76, 175, 80), 0.5));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(117, 95, 160), 0));
            Gradient.Fill = myHorizontalGradient;
            minValueGra.Content = Convert.ToString(Math.Round(umax)+" m/s");
            maxValueGra.Content = Convert.ToString(Math.Round(umin)+ " m/s");
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
            vmax = maxvalue;
            vmin = minvalue;
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
            LinearGradientBrush myHorizontalGradient = new LinearGradientBrush();
            myHorizontalGradient.StartPoint = new Point(0, 0.3);
            myHorizontalGradient.EndPoint = new Point(1, 0.3);
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(240, 255, 0), 1));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(76, 175, 80), 0.5));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(117, 95, 160), 0));
            Gradient.Fill = myHorizontalGradient;
            minValueGra.Content = Convert.ToString(Math.Round(vmax) + " m/s");
            maxValueGra.Content = Convert.ToString(Math.Round(vmin) + " m/s");
        }

        private void OpenVideo_Click(object sender, RoutedEventArgs e)
        {
            string url = "https://youtu.be/Qgto2vXkQaY";
            var psi = new ProcessStartInfo();
            psi.UseShellExecute = true;
            psi.FileName = url;
            Process.Start(psi);
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
            rhomax = maxvalue;
            rhomin = minvalue;
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
            LinearGradientBrush myHorizontalGradient = new LinearGradientBrush();
            myHorizontalGradient.StartPoint = new Point(0, 0.3);
            myHorizontalGradient.EndPoint = new Point(1, 0.3);
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(240, 255, 0), 1));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(76, 175, 80), 0.5));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(117, 95, 160), 0));
            Gradient.Fill = myHorizontalGradient;
            minValueGra.Content = Convert.ToString(Math.Round(rhomax,4) + " kg/m3");
            maxValueGra.Content = Convert.ToString(Math.Round(rhomin,4) + " kg/m3");
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
            pmax = maxvalue;
            pmin = minvalue;
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
            LinearGradientBrush myHorizontalGradient = new LinearGradientBrush();
            myHorizontalGradient.StartPoint = new Point(0, 0.3);
            myHorizontalGradient.EndPoint = new Point(1, 0.3);
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(240, 255, 0), 1));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(76, 175, 80), 0.5));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(117, 95, 160), 0));
            Gradient.Fill = myHorizontalGradient;
            minValueGra.Content = Convert.ToString(Math.Round(pmax) + " Pa");
            maxValueGra.Content = Convert.ToString(Math.Round(pmin) + " Pa");
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
            Tmax = maxvalue;
            Tmin = minvalue;
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
            LinearGradientBrush myHorizontalGradient = new LinearGradientBrush();
            myHorizontalGradient.StartPoint = new Point(0, 0.3);
            myHorizontalGradient.EndPoint = new Point(1, 0.3);
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(240, 255, 0), 1));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(76, 175, 80), 0.5));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(117, 95, 160), 0));
            Gradient.Fill = myHorizontalGradient;
            maxValueGra.Content = Convert.ToString(Math.Round(Tmin)+" K");
            minValueGra.Content = Convert.ToString(Math.Round(Tmax) + " K");
        }

        private void ShowOnGM_Click(object sender, RoutedEventArgs e)
        {
            r = new Rules();
            int result = r.loadRules();
            if (result == 0)
            {
                Run_Click(null, null);
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
            Mmax = maxvalue;
            Mmin = minvalue;
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
            LinearGradientBrush myHorizontalGradient = new LinearGradientBrush();
            myHorizontalGradient.StartPoint = new Point(0, 0.3);
            myHorizontalGradient.EndPoint = new Point(1, 0.3);
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(240, 255, 0), 1));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(76, 175, 80), 0.5));
            myHorizontalGradient.GradientStops.Add(new GradientStop(Color.FromRgb(117, 95, 160), 0));
            Gradient.Fill = myHorizontalGradient;
            maxValueGra.Content = Convert.ToString(Math.Round(Mmin,4));
            minValueGra.Content = Convert.ToString(Math.Round(Mmax,4));
        }



        
        private void polygon_MouseEnter(object sender, MouseEventArgs e)
        {
            //Obtnemos la ubicación del puntero
            Polygon pol = (Polygon)sender;
            Point p = (Point)pol.Tag;
            //Mostramos las etiquetas correspondientes a los estados de las celdas
            Info.Visibility = Visibility.Visible;
            //Actualizamos los valores de fase y temperatura
            uData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getU(), 2)+" m/s";
            vData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getV(), 2)+ " m/s";
            rhoData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getRO(), 2)+ " kg/m3";
            pData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getP(), 2) + " Pa";
            TData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getT(), 2) + " K";
            MData.Content = Math.Round(mesh.GetCell(Convert.ToInt32(p.Y), Convert.ToInt32(p.X)).getM(), 2);
            xData.Content = p.X+1;
            yData.Content = p.Y+1;
        }

        private void polygon_MouseLeave(object sender, MouseEventArgs e)
        {
            Info.Visibility = Visibility.Hidden;
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
    }
}