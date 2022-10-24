namespace FastDC
{
    using System;

    internal class MeshBuffer
    {
        public readonly List<MeshVertex> Vertices = new List<MeshVertex>();

        public int NumVertices => Vertices.Count;

        public readonly List<MeshTriangle> Triangles = new List<MeshTriangle>();

        public int NumTriangles => Triangles.Count;
    }
}
