namespace FastDC
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;
    using System.Threading.Tasks;

    internal class VoxelInfo
    {
        public SortedSet<uint> Voxels = new();

        public SortedDictionary<uint, EdgeInfo> Edges = new();

        public SortedDictionary<uint, int> Indices = new();
    }
}
