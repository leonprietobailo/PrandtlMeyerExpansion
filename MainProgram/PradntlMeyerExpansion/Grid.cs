using System;
using System.Collections.Generic;

namespace PradntlMeyerExpansion
{
    public class Grid
    {
        List<Cell[]> Mesh = new List<Cell[]>();

        public Cell GetCell(int ve, int ho)
        {
            return Mesh[ve][ho];
        }

    }
}

