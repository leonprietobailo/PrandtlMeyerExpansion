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

namespace PradntlMeyerExpansion
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
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
