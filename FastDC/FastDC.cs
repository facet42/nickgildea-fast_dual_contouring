namespace FastDC
{
    using Facet42;

    using GlmSharp;

    using ImGuiNET;

    //using MathNet.Numerics.LinearAlgebra;
    //using MathNet.Numerics.LinearAlgebra.Single;

    using OpenTK.Graphics;
    using OpenTK.Graphics.OpenGL;

    //using Qef;

    using System;
    using System.Collections.Generic;
    using System.Diagnostics;
    using System.Numerics;
    using System.Runtime.Intrinsics;
    //using System.Text;
    //using System.Threading.Tasks;

    internal class FastDC
    {
        const float Scale = 8f;
        const float Distance = 128f;

        private Log? log;

        private ProgramHandle shader;

        //private const int OctreeSize = 64;
        private const int VoxelGridSize = 32;
        private const float VoxelGridOffset = VoxelGridSize / 2;
        private readonly string logPath = "FastDC.log";

        private MeshSimplificationOptions options;
        private ViewerOptions viewerOpts;
        private SuperPrimitiveConfig primConfig = new();

        private readonly List<Mesh> Meshes = new();

        private MeshBuffer meshBuffer = new();

        public bool Initialise()
        {
            this.log = new Log(this.logPath);

            var vertexSource = File.ReadAllText("voxel.vert");
            var fragmentSource = File.ReadAllText("voxel.frag");
            this.shader = Graphics.CreateProgram("FastDC", vertexSource, fragmentSource);

            this.viewerOpts.UseCms = true;
            this.viewerOpts.MeshScale = 3.5f;
            this.options.MaxEdgeSize = 2.5f;

            this.primConfig = ConfigForShape(Primitive.Cylinder);
            InitialiseMeshes();

            return true;
        }

        private void InitialiseMeshes()
        {
            this.meshBuffer = DualContour.GenerateMesh(primConfig, VoxelGridSize, log);
            var mesh = CreateGLMesh(this.meshBuffer, viewerOpts.MeshScale, options);

            this.Meshes.Clear();
            this.Meshes.Add(mesh);
        }

        private static Mesh CreateGLMesh(MeshBuffer buffer, float meshScale, MeshSimplificationOptions options)
        {
            Console.WriteLine($"Simplify iteration: error={options.MaxError}");

            var simplifiedMesh = new MeshBuffer();
            for (int i = 0; i < buffer.NumVertices; i++)
            {
                var vertex = buffer.Vertices[i];
                vertex.Position *= meshScale;
                simplifiedMesh.Vertices.Add(vertex);
            }

            simplifiedMesh.Triangles.AddRange(buffer.Triangles);

            //MeshSimplifier(simplifiedMesh, vec4.Zero, options);

            var mesh = new Mesh();
            mesh.Initialise();
            mesh.UploadData(simplifiedMesh);

            return mesh;
        }

        //private static void MeshSimplifier(MeshBuffer simplifiedMesh, vec4 zero, MeshSimplificationOptions options)
        //{
        //    if (simplifiedMesh is null)
        //    {
        //        throw new ArgumentNullException(nameof(simplifiedMesh));
        //    }

        //    Debug.WriteLine("MeshSimplifier called.");
        //}

        private void Dump(SuperPrimitiveConfig config, ViewerOptions viewer, VoxelInfo voxelInfo, MeshBuffer buffer)
        {
            if (this.log == null)
            {
                return;
            }

            this.log.WriteLine($"R:{config.R} S:{config.S}");

            this.log.WriteLine($"Scale:{viewer.MeshScale}");

            this.log.WriteLine($"Voxels{Environment.NewLine}======");

            foreach (var voxel in voxelInfo.Voxels)
            {
                this.log.WriteLine(voxel.ToString());
            }

            this.log.WriteLine($"Edges{Environment.NewLine}=====");

            foreach (var edge in voxelInfo.Edges)
            {
                this.log.WriteLine($"{edge.Key} P:{edge.Value.Position.ToString(",", "G9")} N:{edge.Value.Normal.ToString(",", "G9")} W:{edge.Value.Winding.ToString().ToLowerInvariant()}");
            }

            this.log.WriteLine($"Indices{Environment.NewLine}=======");

            foreach (var index in voxelInfo.Indices)
            {
                this.log.WriteLine($"{index.Key},{index.Value}");
            }

            this.log.WriteLine($"Vertices{Environment.NewLine}========");

            foreach (var vertex in buffer.Vertices)
            {
                this.log.WriteLine($"{vertex.Position.x:G9},{vertex.Position.y:G9},{vertex.Position.z:G9}");
                this.log.WriteLine($"{vertex.Normal.x:G9},{vertex.Normal.y:G9},{vertex.Normal.z:G9}");
            }

            this.log.WriteLine($"Triangles{Environment.NewLine}=========");

            for (int i = 0; i < buffer.NumTriangles; i++)
            {
                log.WriteLine($"{buffer.Triangles[i].V0}");
                log.WriteLine($"{buffer.Triangles[i].V1}");
                log.WriteLine($"{buffer.Triangles[i].V2}");
            }
        }

        //public static vec3 CalculateCubeQEF(vec4[] normals, vec4[] positions, int count, vec4 meanPoint)
        //{
        //    var n = new vec4[count];
        //    var p = new vec4[count];

        //    for (int i = 0; i < count; i++)
        //    {
        //        n[i] = normals[i];
        //        p[i] = positions[i];
        //    }

        //    var A = DenseMatrix.OfRowArrays(n.Select(e => new[] { e.x, e.y, e.z }).ToArray());
        //    var b = DenseVector.OfArray(normals.Zip(p.Select(v => v - meanPoint), vec4.Dot).ToArray());

        //    var pseudo = A.PseudoInverse();
        //    var leastsquares = pseudo.Multiply(b);

        //    var result = leastsquares + DenseVector.OfArray(new[] { meanPoint.x, meanPoint.y, meanPoint.z });

        //    return new vec3(result[0], result[1], result[2]);
        //}

        //private vec3 SolveQEF(vec4[] p, vec4[] n, int idx)
        //{
        //    var error = -1f;
        //    var position = new Vec3();

        //    if (idx > 2)
        //    {
        //        var ata = new Mat3(0);
        //        var atb = new Vec3(0);
        //        var accum = new Vec4(0);

        //        for (int i = 0; i < idx; i++)
        //        {
        //            var nn = new Vec3(n[i].x, n[i].y, n[i].z);
        //            var pp = new Vec3(p[i].x, p[i].y, p[i].z);
        //            Qef.Add(nn, pp, ref ata, ref atb, ref accum);
        //        }

        //        error = Qef.Solve(ata, atb, accum, out position);
        //    }

        //    if (idx <= 2 || Math.Abs(error) > 0f)
        //    {
        //        var pos = vec3.Zero;
        //        for (int i = 0; i < idx; i++)
        //        {
        //            pos += p[i].xyz;
        //        }

        //        position = new Vec3(pos.x, pos.y, pos.z) / idx;
        //    }

        //    return new vec3(position.x, position.y, position.z);
        //}

        private static string ToString(Vector4 v)
        {
            return $"{v.X:#0.000000000},{v.Y:#0.000000000},{v.Z:#0.000000000},{v.W:#0.000000000}";
        }

        private static string ToString(Vector128<float> v)
        {
            return $"{v.GetElement(0):#0.000000000},{v.GetElement(1):#0.000000000},{v.GetElement(2):#0.000000000},{v.GetElement(3):#0.000000000}";
        }

        private static byte MM_Shuffle(byte fp3, byte fp2, byte fp1, byte fp0)
        {
            return ((byte)(((fp3) << 6) | ((fp2) << 4) | ((fp1) << 2) | ((fp0))));
        }

        /// <summary>
        /// Solve least square solution using Ax = b
        /// </summary>
        /// <param name="component"></param>
        /// <param name="corner"></param>
        /// <returns></returns>
        //private Vector3? CalculateJunction(IList<vec4> vertices, IList<vec4> norms, int count, bool corner)
        //{
        //    var normals = new double[count * 3];
        //    var dots = new double[count];

        //    for (var i = 0; i < count; i++)
        //    {
        //        normals[count * 0 + i] = norms[i].x;
        //        normals[count * 1 + i] = norms[i].y;
        //        normals[count * 2 + i] = norms[i].z;

        //        // M* dot(location, normal)
        //        dots[i] = vertices[i].x * norms[i].x +
        //                  vertices[i].y * norms[i].y +
        //                  vertices[i].z * norms[i].z;
        //    }

        //    var A = Matrix<double>.Build.Dense(count, 3, normals);
        //    var b = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(dots);

        //    var x = A.Transpose().Multiply(A).Inverse().Multiply(A.Transpose().Multiply(b));
        //    if (double.IsNaN(x[0]) || double.IsNaN(x[1]) || double.IsNaN(x[2]))
        //    {
        //        return null;
        //    }

        //    return new Vector3((float)x[0], (float)x[1], (float)x[2]);
        //}

        //private float FindIntersection(float d0, float d1)
        //{
        //    var t = -d0 / (d1 - d0);

        //    return t;
        //}

        private static SuperPrimitiveConfig ConfigForShape(Primitive primitive)
        {
            var config = new SuperPrimitiveConfig();
            config.Scale = Scale;

            switch (primitive)
            {
                case Primitive.Cylinder:
                    config.S = new vec4(1f);
                    config.R = new vec2(1f, 0f);
                    break;

                case Primitive.Pill:
                    config.S = new vec4(1f, 1f, 2f, 1f);
                    config.R = new vec2(1f);
                    break;

                case Primitive.Corridor:
                    config.S = new vec4(1f, 1f, 1f, 0.25f);
                    config.R = new vec2(0.1f);
                    break;

                case Primitive.Torus:
                    config.S = new vec4(1f, 1f, 0.25f, 0.25f);
                    config.R = new vec2(1f, 0.25f);
                    break;

                case Primitive.Cube:
                default:
                    config.S = new vec4(1f);
                    config.R = new vec2(0f);
                    break;
            }

            return config;
        }

        public void Dispose()
        {
            foreach (var mesh in this.Meshes)
            {
                mesh.Dispose();
            }

            Graphics.DeleteProgram(this.shader);

            this.log?.Dispose();
        }

        public void Update(double deltaTime)
        {
            if (this.viewerOpts.RefreshModel == true)
            {
                this.viewerOpts.RefreshModel = false;

                InitialiseMeshes();
            }
        }

        public void Render()
        {
            var dir = new vec3(0f, 0f, 1f);
            // TODO: Rotate direction

            var position = dir * Distance;
            //TestFrame(this.Meshes[0]);
            DrawFrame(this.shader, this.Meshes, position, -dir, viewerOpts.DrawWireframe, viewerOpts.MeshScale);
        }

        private void DrawFrame(ProgramHandle shader, List<Mesh> meshes, vec3 position, vec3 forward, bool drawWireframe, float meshScale)
        {
            GL.Enable(EnableCap.DepthTest);
            CheckGLError();
            GL.Enable(EnableCap.CullFace);
            CheckGLError();
            GL.CullFace(CullFaceMode.Back);
            CheckGLError();

            var projection = mat4.Perspective(80f, 16f / 9f, 0.1f, 500f);
            var view = mat4.LookAt(position + forward, vec3.Zero, vec3.UnitY);
            var model = mat4.RotateX(glm.Radians(-30f)) * mat4.RotateY(glm.Radians(45f));

            GL.UseProgram(this.shader);
            CheckGLError();

            var offset = -1f * (meshes.Count / 2);
            var o = (((meshes.Count & 1) == 0) ? 0.5f : 0f);
            offset += o;
            for (int i = 0; i < meshes.Count; i++, offset++)
            {
                var mesh = meshes[i];
                //model = mat4.Translate((float)offset * new vec3(meshScale / 2f, 0f, 0f));
                var m = projection * view * (mat4.Translate(offset * VoxelGridSize * meshScale, 0, 1) * model);
                var mvp = new OpenTK.Mathematics.Matrix4(m.m00, m.m01, m.m02, m.m03, m.m10, m.m11, m.m12, m.m13, m.m20, m.m21, m.m22, m.m23, m.m30, m.m31, m.m32, m.m33);

                var mat = GL.GetUniformLocation(this.shader, "MVP");
                CheckGLError();
                GL.UniformMatrix4f(mat, false, mvp);
                CheckGLError();

                var solid = GL.GetUniformLocation(this.shader, "solidColour");
                CheckGLError();
                GL.Uniform1i(solid, 0);
                CheckGLError();

                GL.BindVertexArray(mesh.VertexArrayObj);
                CheckGLError();
                GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);
                CheckGLError();
                GL.DrawElements(PrimitiveType.Triangles, mesh.NumIndices, DrawElementsType.UnsignedInt, 0);
                CheckGLError();

                if (drawWireframe == true)
                {
                    GL.Uniform1i(solid, 1);
                    CheckGLError();

                    GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Line);
                    CheckGLError();
                    GL.DrawElements(PrimitiveType.Triangles, mesh.NumIndices, DrawElementsType.UnsignedInt, 0);
                    CheckGLError();
                    GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);
                    CheckGLError();
                }
            }

            GL.BindVertexArray(VertexArrayHandle.Zero);
            CheckGLError();
            GL.UseProgram(ProgramHandle.Zero);
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

        public void RenderUserInterface()
        {
            GuiDrawFrame(ref viewerOpts, ref options, ref primConfig);
        }

        private static void GuiDrawFrame(ref ViewerOptions viewerOpts, ref MeshSimplificationOptions options, ref SuperPrimitiveConfig primConfig)
        {
            ImGui.GetIO().MouseDrawCursor = false;

            ImGui.Begin("Options");

            if (ImGui.CollapsingHeader("Mesh Simplification Options"))
            {
                ImGui.SliderFloat("Random Edge Fraction", ref options.EdgeFraction, 0f, 1f);
                ImGui.SliderInt("Max Iterations", ref options.MaxIterations, 1, 100);
                ImGui.SliderFloat("Target Triangle Percentage", ref options.TargetPercentage, 0f, 1f, "%.3f", ImGuiSliderFlags.Logarithmic/* 1.5f*/);
                ImGui.SliderFloat("Max QEF Error", ref options.MaxError, 0f, 10f, "%.3f", ImGuiSliderFlags.Logarithmic/*1.5f*/);
                ImGui.SliderFloat("Max Edge Size", ref options.MaxEdgeSize, 0f, 10f);
                ImGui.SliderFloat("Min Angle Cosine", ref options.MinAngleCosine, 0f, 1f);
            }

            if (ImGui.CollapsingHeader("Super Primitive Config"))
            {
                if (ImGui.Button("Cube"))
                {
                    primConfig = ConfigForShape(Primitive.Cube);
                    viewerOpts.RefreshModel = true;
                }

                ImGui.SameLine();
                if (ImGui.Button("Torus"))
                {
                    primConfig = ConfigForShape(Primitive.Torus);
                    viewerOpts.RefreshModel = true;
                }

                ImGui.SameLine();
                if (ImGui.Button("Cylinder"))
                {
                    primConfig = ConfigForShape(Primitive.Cylinder);
                    viewerOpts.RefreshModel = true;
                }

                ImGui.SameLine();
                if (ImGui.Button("Pill"))
                {
                    primConfig = ConfigForShape(Primitive.Pill);
                    viewerOpts.RefreshModel = true;
                }

                ImGui.SameLine();
                if (ImGui.Button("Corridor"))
                {
                    primConfig = ConfigForShape(Primitive.Corridor);
                    viewerOpts.RefreshModel = true;
                }

                //ImGui.SliderFloat4("S", ref primConfig.S, 0f, 2f);
                //ImGui.SliderFloat2("R", ref primConfig.R, 0f, 1f);
            }

            if (ImGui.CollapsingHeader("Viewer Options"))
            {
                ImGui.SliderFloat("Mesh Scale", ref viewerOpts.MeshScale, 1f, 5f, "%.3f", ImGuiSliderFlags.Logarithmic/*1.2f*/);

                if (ImGui.RadioButton("Draw Wireframe", viewerOpts.DrawWireframe))
                {
                    viewerOpts.DrawWireframe = !viewerOpts.DrawWireframe;
                }

                if (ImGui.RadioButton("Use CMS", viewerOpts.UseCms))
                {
                    viewerOpts.UseCms = !viewerOpts.UseCms;
                    viewerOpts.RefreshModel = true;
                }
            }

            if (ImGui.Button("Refresh"))
            {
                viewerOpts.RefreshModel = true;
            }

            ImGui.End();
        }
    }
}
