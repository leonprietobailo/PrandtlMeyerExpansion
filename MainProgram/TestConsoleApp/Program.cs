using PradntlMeyerExpansion;
using System;


namespace TestConsoleApp
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Hello World!");
            Rules r = new Rules(678, 0, 1.23, 0.101e6, 286.1, 2, 0.5, 1.4, 287, 10, 5.352 * Math.PI / 180, 41, 65, 40, 0.5);
            Grid g = new Grid(r);
            g.PrandtlMeyerExpansion();



        }
    }
}
