namespace FastDC
{
    using GlmSharp;

    using OpenTK.Graphics;
    using OpenTK.Graphics.OpenGL;

    using System;
    using System.Diagnostics;

    internal class Mesh
    {
        public VertexArrayHandle VertexArrayObj;
        public BufferHandle VertexBuffer;
        public BufferHandle IndexBuffer;
        public int NumIndices = 0;

        internal void Dispose()
        {
            GL.DeleteBuffer(IndexBuffer);
            CheckGLError();
            GL.DeleteBuffer(VertexBuffer);
            CheckGLError();
            GL.DeleteVertexArray(VertexArrayObj);
            CheckGLError();
        }

        internal void Initialise()
        {
            VertexArrayObj = GL.GenVertexArray();
            CheckGLError();
            VertexBuffer = GL.GenBuffer();
            CheckGLError();
            IndexBuffer = GL.GenBuffer();
            CheckGLError();

            GL.BindVertexArray(VertexArrayObj);
            CheckGLError();
            GL.BindBuffer(BufferTargetARB.ArrayBuffer, VertexBuffer);
            CheckGLError();
            GL.BindBuffer(BufferTargetARB.ElementArrayBuffer, IndexBuffer);
            CheckGLError();

            GL.EnableVertexAttribArray(0);
            CheckGLError();
            GL.VertexAttribPointer(0, 3, VertexAttribPointerType.Float, false, 24, 0);
            CheckGLError();

            GL.EnableVertexAttribArray(1);
            CheckGLError();
            GL.VertexAttribPointer(1, 3, VertexAttribPointerType.Float, false, 24, 12);
            CheckGLError();

            //GL.EnableVertexAttribArray(2);
            //CheckGLError();
            //GL.VertexAttribPointer(2, 4, VertexAttribPointerType.Float, false, 48, 32);
            //CheckGLError();

            GL.BindVertexArray(VertexArrayHandle.Zero);
            CheckGLError();
        }

        internal void UploadData(MeshBuffer buffer)
        {
            GL.BindVertexArray(this.VertexArrayObj);
            CheckGLError();

            GL.BindBuffer(BufferTargetARB.ArrayBuffer, this.VertexBuffer);
            CheckGLError();
            GL.BufferData(BufferTargetARB.ArrayBuffer, buffer.Vertices.ToArray(), BufferUsageARB.StaticDraw);
            CheckGLError();

            GL.BindBuffer(BufferTargetARB.ElementArrayBuffer, this.IndexBuffer);
            CheckGLError();
            GL.BufferData(BufferTargetARB.ElementArrayBuffer, buffer.Triangles.ToArray(), BufferUsageARB.StaticDraw);
            CheckGLError();

            Console.WriteLine($"Mesh: {buffer.NumVertices} vertices");
            this.NumIndices = 3 * buffer.NumTriangles;

            GL.BindVertexArray(VertexArrayHandle.Zero);
            CheckGLError();
        }

        private static void CheckGLError()
        {
            var error = GL.GetError();
            if (error != ErrorCode.NoError)
            {
                Debug.WriteLine(error);
            }
        }
    }
}
