namespace FastDC
{
    using GlmSharp;

    using ImGuiNET;

    using MathNet.Numerics.LinearAlgebra.Single;

    //using MathNet.Numerics.LinearAlgebra;
    //using MathNet.Numerics.LinearAlgebra.Single;

    using OpenTK.Graphics;
    using OpenTK.Graphics.OpenGL;

    //using Qef;

    using System;
    using System.Collections.Generic;
    using System.Diagnostics;
    using System.Diagnostics.CodeAnalysis;
    using System.Linq;
    using System.Numerics;
    using System.Runtime.InteropServices;
    using System.Runtime.Intrinsics;
    using System.Runtime.Intrinsics.X86;
    //using System.Text;
    //using System.Threading.Tasks;

    internal partial class FastDC
    {
        const int SvdNumSweeps = 5;
        const float PseudoInverseThreshold = 0.001f;

        const float scale = 32f;    // TODO: Move / remove this
        const float Distance = 250f;

        private Log? log;

        private ProgramHandle shader;

        //private const int OctreeSize = 64;
        private const int VoxelGridSize = 128;
        private const float VoxelGridOffset = 64f;
        private readonly string logPath = "FastDC.log";
        private static readonly vec4[] AxisOffset =
        {
            new vec4(1f, 0f, 0f, 0f),
            new vec4(0f, 1f, 0f, 0f),
            new vec4(0f, 0f, 1f, 0f)
        };

        private static readonly ivec4[][] EdgeNodeOffsets =
        {
            new[] { new ivec4(0), new ivec4(0, 0, 1, 0), new ivec4(0, 1, 0, 0), new ivec4(0, 1, 1, 0) },
            new[] { new ivec4(0), new ivec4(1, 0, 0, 0), new ivec4(0, 0, 1, 0), new ivec4(1, 0, 1, 0) },
            new[] { new ivec4(0), new ivec4(0, 1, 0, 0), new ivec4(1, 0, 0, 0), new ivec4(1, 1, 0, 0) },
        };

        private static readonly uint[] EncodedEdgeNodeOffsets =
        {
            0x00000000,
            0x00100000,
            0x00000400,
            0x00100400,
            0x00000000,
            0x00000001,
            0x00100000,
            0x00100001,
            0x00000000,
            0x00000400,
            0x00000001,
            0x00000401,
        };

        private static readonly uint[] EncodedEdgeOffsets =
        {
            0x00000000,
            0x00100000,
            0x00000400,
            0x00100400,
            0x40000000,
            0x40100000,
            0x40000001,
            0x40100001,
            0x80000000,
            0x80000400,
            0x80000001,
            0x80000401,
        };

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

            this.viewerOpts.MeshScale = 3.5f;
            this.options.MaxEdgeSize = 2.5f;

            this.primConfig = ConfigForShape(Primitive.Cylinder);
            InitialiseMeshes();

            return true;
        }

        private void InitialiseMeshes()
        {
            this.meshBuffer = GenerateMesh(primConfig);
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

        private MeshBuffer GenerateMesh(SuperPrimitiveConfig config)
        {
            var voxelInfo = FindActiveVoxels(config);

            var buffer = new MeshBuffer();

            GenerateVertexData(voxelInfo, ref buffer);
            GenerateTriangles(voxelInfo, ref buffer);

            this.Dump(config, viewerOpts, voxelInfo, buffer);

            Console.WriteLine($"Mesh: {buffer.NumVertices} {buffer.NumTriangles}");

            return buffer;
        }

        private void Dump(SuperPrimitiveConfig config, ViewerOptions viewer, VoxelInfo voxelInfo, MeshBuffer buffer)
        {
            if (this.log == null)
            {
                return;
            }

            this.log.WriteLine($"R:{config.R} S:{config.S}");

            this.log.WriteLine($"Scale:{viewer.MeshScale}");

            this.log.WriteLine("Voxels\n======");

            foreach (var voxel in voxelInfo.Voxels)
            {
                this.log.WriteLine(voxel.ToString());
            }

            this.log.WriteLine("Edges\n=====");

            foreach (var edge in voxelInfo.Edges)
            {
                this.log.WriteLine($"{edge.Key} P:{edge.Value.Position.ToString(",", "G9")} N:{edge.Value.Normal.ToString(",", "G9")} W:{edge.Value.Winding.ToString().ToLowerInvariant()}");
            }

            this.log.WriteLine("Indices\n=======");

            foreach (var index in voxelInfo.Indices)
            {
                this.log.WriteLine($"{index.Key},{index.Value}");
            }

            this.log.WriteLine("Vertices\n========");

            foreach (var vertex in buffer.Vertices)
            {
                this.log.WriteLine($"{vertex.Position.x:G9},{vertex.Position.y:G9},{vertex.Position.z:G9}");
                this.log.WriteLine($"{vertex.Normal.x:G9},{vertex.Normal.y:G9},{vertex.Normal.z:G9}");
            }

            this.log.WriteLine("Triangles\n=========");

            for (int i = 0; i < buffer.NumTriangles; i++)
            {
                log.WriteLine($"{buffer.Triangles[i].V0}");
                log.WriteLine($"{buffer.Triangles[i].V1}");
                log.WriteLine($"{buffer.Triangles[i].V2}");
            }
        }

        private void GenerateVertexData(VoxelInfo info, ref MeshBuffer buffer)
        {
            log?.WriteLine("GenerateVertexData\n==================");
            log?.WriteLine(info.Voxels.Count.ToString());

            var idxCounter = 0;
            foreach (var voxelId in info.Voxels)
            {
                log?.WriteLine($"VID:{voxelId}");

                var p = new vec4[12];
                var n = new vec4[12];

                var idx = 0;
                for (int i = 0; i < 12; i++)
                {
                    var edgeId = voxelId + EncodedEdgeOffsets[i];
                    log?.WriteLine($"EID:{edgeId}");

                    if (info.Edges.ContainsKey(edgeId) == true)
                    {
                        p[idx] = info.Edges[edgeId].Position;
                        n[idx] = info.Edges[edgeId].Normal;
                        idx++;

                        log?.WriteLine($"P:{p[idx].x:G9},{p[idx].y:G9},{p[idx].z:G9},{p[idx].w:G9}");
                        log?.WriteLine($"N:{n[idx].x:G9},{n[idx].y:G9},{n[idx].z:G9},{n[idx].w:G9}");
                    }
                }

                var error = QefSolveFromPoints4d(p, n, idx, out var nodePos);
                var nodeNormal = vec3.Zero;
                for (int i = 0; i < idx; i++)
                {
                    nodeNormal += n[i].xyz;
                }
                nodeNormal *= (1f / (float)idx);

                log?.WriteLine($"NP:{nodePos.x:G9},{nodePos.y:G9},{nodePos.z:G9}");
                log?.WriteLine($"NN:{nodeNormal.x:G9},{nodeNormal.y:G9},{nodeNormal.z:G9}");

                info.Indices.Add(voxelId, idxCounter++);

                buffer.Vertices.Add(new MeshVertex() { Position = nodePos, Normal = nodeNormal });
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

        private float QefSolveFromPoints4d(vec4[] positions, vec4[] normals, int count, out vec3 position)
        {
            if (count < 2 || count > 12)
            {
                position = vec3.Zero;
                return 0f;
            }

            var p = new Vector4[count];
            var n = new Vector4[count];
            for (int i = 0; i < count; i++)
            {
                p[i] = new Vector4(positions[i].x, positions[i].y, positions[i].z, positions[i].w);
                n[i] = new Vector4(normals[i].x, normals[i].y, normals[i].z, normals[i].w);
            }

            var error = QefSolveFromPoints(p, n, count, out var solved);
            position = new vec3(solved.X, solved.Y, solved.Z);
            return error;
        }

        private float QefSolveFromPoints(Vector4[] positions, Vector4[] normals, int count, out Vector4 solvedPosition)
        {
            var pointAccum = new Vector4();
            var ATb = new Vector4();

            var ATA = new Mat4x4();

            for (int i = 0; i < count && i < 12; i++)
            {
                QefSimdAdd(positions[i], normals[i], ref ATA, ref ATb, ref pointAccum);
            }

            log?.WriteLine($"ATb:{ATb.X:0.000000000},{ATb.Y:0.000000000},{ATb.Z:0.000000000},{ATb.W:0.000000000}");

            return QefSimdSolve(ATA, ATb, pointAccum, out solvedPosition);
        }

        private static void QefSimdAdd(Vector4 p, Vector4 n, ref Mat4x4 ATA, ref Vector4 ATb, ref Vector4 pointAccum)
        {
            //var nX = Sse.Multiply(Sse.Shuffle(n, n, MM_Shuffle(0, 0, 0, 0)), n);
            //var nY = Sse.Multiply(Sse.Shuffle(n, n, MM_Shuffle(1, 1, 1, 1)), n);
            //var nZ = Sse.Multiply(Sse.Shuffle(n, n, MM_Shuffle(2, 2, 2, 2)), n);

            var nX = new Vector4(n.X, n.X, n.X, n.X) * n;
            var nY = new Vector4(n.Y, n.Y, n.Y, n.Y) * n;
            var nZ = new Vector4(n.Z, n.Z, n.Z, n.Z) * n;

            //log?.WriteLine($"n: {n.X:0.000000000},{n.Y:0.000000000},{n.Z:0.000000000},{n.W:0.000000000}");
            //log?.WriteLine($"nX: {nX.X:0.000000000},{nX.Y:0.000000000},{nX.Z:0.000000000},{nX.W:0.000000000}");
            //log?.WriteLine($"nY: {nY.X:0.000000000},{nY.Y:0.000000000},{nY.Z:0.000000000},{nY.W:0.000000000}");
            //log?.WriteLine($"nZ: {nZ.X:0.000000000},{nZ.Y:0.000000000},{nZ.Z:0.000000000},{nZ.W:0.000000000}");

            ATA.row[0] += nX;
            ATA.row[1] += nY;
            ATA.row[2] += nZ;

            var d = Vector4.Dot(p, n);
            var x = new Vector4(d, d, d, 0f);
            x *= n;
            ATb += x;
            pointAccum += p;
        }

        private float QefSimdSolve(Mat4x4 ATA, Vector4 ATb, Vector4 pointAccum, out Vector4 solved)
        {
            var massPoint = pointAccum / pointAccum.W;

            log?.WriteLine($"pointAccum: {pointAccum.X:0.000000000},{pointAccum.Y:0.000000000},{pointAccum.Z:0.000000000},{pointAccum.W:0.000000000}");
            log?.WriteLine($"MassPoint: {massPoint.X:0.000000000},{massPoint.Y:0.000000000},{massPoint.Z:0.000000000},{massPoint.W:0.000000000}");

            var p = ATb - Vector4.Transform(massPoint, ATA.AsMatrix4x4());

            SvdSolveATAATb(ATA, p, out solved);

            var error = QefSimdCalcError(ATA, solved, ATb);

            solved += massPoint;

            log?.WriteLine($"Solved: {solved.X:0.000000000},{solved.Y:0.000000000},{solved.Z:0.000000000},{solved.W:0.000000000}");

            return error;
        }

        private static float QefSimdCalcError(Mat4x4 A, Vector4 x, Vector4 b)
        {
            var tmp = b - Vector4.Transform(x, A.AsMatrix4x4());
            return Vector4.Dot(tmp, tmp);
        }

        private void SvdSolveATAATb(Mat4x4 ATA, Vector4 ATb, out Vector4 x)
        {
            var V = new Mat4x4(
                1f, 0, 0, 0,
                0, 1f, 0, 0,
                0, 0, 1f, 0,
                0, 0, 0, 0
                );

            var sigma = SvdSolveSym(ref V, ATA);

            log?.WriteLine($"Sigma: {sigma.X:#0.000000000},{sigma.Y:#0.000000000},{sigma.Z:#0.000000000},{sigma.W:#0.000000000}");

            log?.WriteLine($"V[0]:{V[0, 0]:#0.000000000},{V[0, 1]:#0.000000000},{V[0, 2]:#0.000000000},{V[0, 3]:#0.000000000}");
            log?.WriteLine($"V[1]:{V[1, 0]:#0.000000000},{V[1, 1]:#0.000000000},{V[1, 2]:#0.000000000},{V[1, 3]:#0.000000000}");
            log?.WriteLine($"V[2]:{V[2, 0]:#0.000000000},{V[2, 1]:#0.000000000},{V[2, 2]:#0.000000000},{V[2, 3]:#0.000000000}");
            log?.WriteLine($"V[3]:{V[3, 0]:#0.000000000},{V[3, 1]:#0.000000000},{V[3, 2]:#0.000000000},{V[3, 3]:#0.000000000}");

            SvdPseudoInverse(out var Vinv, sigma, V);

            //log?.WriteLine($"Vinv[0]:{Vinv.M11:#0.000000000},{Vinv.M12:#0.000000000},{Vinv.M13:#0.000000000},{Vinv.M14:#0.000000000}");
            //log?.WriteLine($"Vinv[1]:{Vinv.M21:#0.000000000},{Vinv.M22:#0.000000000},{Vinv.M23:#0.000000000},{Vinv.M24:#0.000000000}");
            //log?.WriteLine($"Vinv[2]:{Vinv.M31:#0.000000000},{Vinv.M32:#0.000000000},{Vinv.M33:#0.000000000},{Vinv.M34:#0.000000000}");
            //log?.WriteLine($"Vinv[3]:{Vinv.M41:#0.000000000},{Vinv.M42:#0.000000000},{Vinv.M43:#0.000000000},{Vinv.M44:#0.000000000}");
            log?.WriteLine($"Vinv:{Vinv}");

            x = Vector4.Transform(ATb, Vinv.AsMatrix4x4());
        }

        private void SvdPseudoInverse(out Mat4x4 o, Vector4 sigma, Mat4x4 v)
        {
            var invDet = SvdInvDet(sigma);
            log?.WriteLine($"InvDet:{invDet}");

            //log?.WriteLine($"v[0]:{v.M11:#0.000000000},{v.M12:#0.000000000},{v.M13:#0.000000000},{v.M14:#0.000000000}");
            //log?.WriteLine($"v[1]:{v.M21:#0.000000000},{v.M22:#0.000000000},{v.M23:#0.000000000},{v.M24:#0.000000000}");
            //log?.WriteLine($"v[2]:{v.M31:#0.000000000},{v.M32:#0.000000000},{v.M33:#0.000000000},{v.M34:#0.000000000}");
            //log?.WriteLine($"v[3]:{v.M41:#0.000000000},{v.M42:#0.000000000},{v.M43:#0.000000000},{v.M44:#0.000000000}");
            log?.WriteLine($"v:{v}");

            var m = new Mat4x4();
            m.row[0] = v.row[0] * invDet;
            m.row[1] = v.row[1] * invDet;
            m.row[2] = v.row[2] * invDet;
            m.row[3] = new Vector4();

            log?.WriteLine($"m[0]:{ToString(m.row[0].AsVector128())}");
            log?.WriteLine($"m[1]:{ToString(m.row[1].AsVector128())}");
            log?.WriteLine($"m[2]:{ToString(m.row[2].AsVector128())}");
            log?.WriteLine($"m[3]:{ToString(m.row[3].AsVector128())}");

            o = new Mat4x4();

            o[0, 0] = Vector4.Dot(m.row[0], v.row[0]);
            o[0, 1] = Vector4.Dot(m.row[1], v.row[0]);
            o[0, 2] = Vector4.Dot(m.row[2], v.row[0]);
            o[0, 3] = 0f;

            o[1, 0] = Vector4.Dot(m.row[0], v.row[1]);
            o[1, 1] = Vector4.Dot(m.row[1], v.row[1]);
            o[1, 2] = Vector4.Dot(m.row[2], v.row[1]);
            o[1, 3] = 0f;

            o[2, 0] = Vector4.Dot(m.row[0], v.row[2]);
            o[2, 1] = Vector4.Dot(m.row[1], v.row[2]);
            o[2, 2] = Vector4.Dot(m.row[2], v.row[2]);
            o[2, 3] = 0f;

            o[3, 0] = m.row[3].X;
            o[3, 1] = m.row[3].Y;
            o[3, 2] = m.row[3].Z;
            o[3, 3] = m.row[3].W;

            //log?.WriteLine($"o[0]:{o.M11:#0.000000000},{o.M12:#0.000000000},{o.M13:#0.000000000},{o.M14:#0.000000000}");
            //log?.WriteLine($"o[1]:{o.M21:#0.000000000},{o.M22:#0.000000000},{o.M23:#0.000000000},{o.M24:#0.000000000}");
            //log?.WriteLine($"o[2]:{o.M31:#0.000000000},{o.M32:#0.000000000},{o.M33:#0.000000000},{o.M34:#0.000000000}");
            //log?.WriteLine($"o[3]:{o.M41:#0.000000000},{o.M42:#0.000000000},{o.M43:#0.000000000},{o.M44:#0.000000000}");
            log?.WriteLine($"o:{o}");
        }

        private static string ToString(Vector128<float> v)
        {
            return $"{v.GetElement(0):#0.000000000},{v.GetElement(1):#0.000000000},{v.GetElement(2):#0.000000000},{v.GetElement(3):#0.000000000}";
        }

        private Vector4 SvdInvDet(Vector4 v)
        {
            var ones = new Vector4(1f);
            var tol = new Vector4(PseudoInverseThreshold);

            var abs_x = Vector4.Abs(v);
            var one_over_x = ones / v;
            var abs_one_over_x = Vector4.Abs(one_over_x);
            var min_abs = Vector4.Min(abs_x, abs_one_over_x);

            return new Vector4(
                min_abs.X >= tol.X ? one_over_x.X : 0,
                min_abs.Y >= tol.X ? one_over_x.Y : 0,
                min_abs.Z >= tol.X ? one_over_x.Z : 0,
                min_abs.W >= tol.X ? one_over_x.W : 0);
        }

        private Vector4 SvdSolveSym(ref Mat4x4 v, Mat4x4 a)
        {
            var vtav = a;

            log?.WriteLine("svd_solve_sym");
            log?.WriteLine($"vtav[0]: {vtav.row[0].X:#0.000000000},{vtav.row[0].Y:#0.000000000},{vtav.row[0].Z:#0.000000000},{vtav.row[0].W:#0.000000000}");
            log?.WriteLine($"vtav[1]: {vtav.row[1].X:#0.000000000},{vtav.row[1].Y:#0.000000000},{vtav.row[1].Z:#0.000000000},{vtav.row[1].W:#0.000000000}");
            log?.WriteLine($"vtav[2]: {vtav.row[2].X:#0.000000000},{vtav.row[2].Y:#0.000000000},{vtav.row[2].Z:#0.000000000},{vtav.row[2].W:#0.000000000}");
            log?.WriteLine($"vtav[3]: {vtav.row[3].X:#0.000000000},{vtav.row[3].Y:#0.000000000},{vtav.row[3].Z:#0.000000000},{vtav.row[3].W:#0.000000000}");

            for (int i = 0; i < SvdNumSweeps; ++i)
            {
                var c = new Vector4();
                var s = new Vector4();

                if (vtav.row[0].Y != 0f)
                {
                    GivensCoeffsSym(ref c, ref s, vtav, 0, 1);
                    RotateQXY(ref vtav, c, s, 0, 1);
                    RotateXY(ref vtav, ref v, c.Y, s.Y, 0, 1);
                    vtav.row[0].Y = 0f;
                }

                if (vtav.row[0].Z != 0f)
                {
                    GivensCoeffsSym(ref c, ref s, vtav, 0, 2);
                    RotateQXY(ref vtav, c, s, 0, 2);
                    RotateXY(ref vtav, ref v, c.Y, s.Y, 0, 2);
                    vtav.row[0].Z = 0f;
                }

                if (vtav.row[1].Z != 0f)
                {
                    GivensCoeffsSym(ref c, ref s, vtav, 1, 2);
                    RotateQXY(ref vtav, c, s, 1, 2);
                    RotateXY(ref vtav, ref v, c.Z, s.Z, 1, 2);
                    vtav.row[1].Z = 0f;
                }

                log?.WriteLine($"c:{c}");
                log?.WriteLine($"s:{s}");
            }

            log?.WriteLine($"vtav[0]: {vtav.row[0].X:#0.000000000},{vtav.row[0].Y:#0.000000000},{vtav.row[0].Z:#0.000000000},{vtav.row[0].W:#0.000000000}");
            log?.WriteLine($"vtav[1]: {vtav.row[1].X:#0.000000000},{vtav.row[1].Y:#0.000000000},{vtav.row[1].Z:#0.000000000},{vtav.row[1].W:#0.000000000}");
            log?.WriteLine($"vtav[2]: {vtav.row[2].X:#0.000000000},{vtav.row[2].Y:#0.000000000},{vtav.row[2].Z:#0.000000000},{vtav.row[2].W:#0.000000000}");
            log?.WriteLine($"vtav[3]: {vtav.row[3].X:#0.000000000},{vtav.row[3].Y:#0.000000000},{vtav.row[3].Z:#0.000000000},{vtav.row[3].W:#0.000000000}");

            return new Vector4(vtav.row[0].X, vtav.row[1].Y, vtav.row[2].Z, 0f);
        }

        public static float Item(Mat4x4 m, int row, int column)
        {
            return m[row, column];
            //switch (row)
            //{
            //    case 0:
            //        switch (column)
            //        {
            //            case 0:
            //                return m.M11;
            //            case 1:
            //                return m.M12;
            //            case 2:
            //                return m.M13;
            //            case 3:
            //                return m.M14;
            //        }
            //        break;

            //    case 1:
            //        switch (column)
            //        {
            //            case 0:
            //                return m.M21;
            //            case 1:
            //                return m.M22;
            //            case 2:
            //                return m.M23;
            //            case 3:
            //                return m.M24;
            //        }
            //        break;

            //    case 2:
            //        switch (column)
            //        {
            //            case 0:
            //                return m.M31;
            //            case 1:
            //                return m.M32;
            //            case 2:
            //                return m.M33;
            //            case 3:
            //                return m.M34;
            //        }
            //        break;

            //    case 3:
            //        switch (column)
            //        {
            //            case 0:
            //                return m.M41;
            //            case 1:
            //                return m.M42;
            //            case 2:
            //                return m.M43;
            //            case 3:
            //                return m.M44;
            //        }
            //        break;
            //}

            //throw new IndexOutOfRangeException();
        }

        public static void Item(ref Mat4x4 m, int row, int column, float value)
        {
            m[row, column] = value;

            //switch (row)
            //{
            //    case 0:
            //        switch (column)
            //        {
            //            case 0:
            //                m.M11 = value;
            //                return;
            //            case 1:
            //                m.M12 = value;
            //                return;
            //            case 2:
            //                m.M13 = value;
            //                return;
            //            case 3:
            //                m.M14 = value;
            //                return;
            //        }
            //        break;

            //    case 1:
            //        switch (column)
            //        {
            //            case 0:
            //                m.M21 = value;
            //                return;
            //            case 1:
            //                m.M22 = value;
            //                return;
            //            case 2:
            //                m.M23 = value;
            //                return;
            //            case 3:
            //                m.M24 = value;
            //                return;
            //        }
            //        break;

            //    case 2:
            //        switch (column)
            //        {
            //            case 0:
            //                m.M31 = value;
            //                return;
            //            case 1:
            //                m.M32 = value;
            //                return;
            //            case 2:
            //                m.M33 = value;
            //                return;
            //            case 3:
            //                m.M34 = value;
            //                return;
            //        }
            //        break;

            //    case 3:
            //        switch (column)
            //        {
            //            case 0:
            //                m.M41 = value;
            //                return;
            //            case 1:
            //                m.M42 = value;
            //                return;
            //            case 2:
            //                m.M43 = value;
            //                return;
            //            case 3:
            //                m.M44 = value;
            //                return;
            //        }
            //        break;
            //}

            //throw new IndexOutOfRangeException();
        }

        private void RotateXY(ref Mat4x4 vtav, ref Mat4x4 v, float c, float s, int a, int b)
        {
            var simd_u = new Vector4(Item(v, 0, a), Item(v, 1, a), Item(v, 2, a), vtav[0, 3 - b]);
            var simd_v = new Vector4(Item(v, 0, b), Item(v, 1, b), Item(v, 2, b), vtav[1 - a, 2]);

            log?.WriteLine("rotate_xy");
            log?.WriteLine($"a:{a} b:{b}");
            log?.WriteLine($"vtav[0]: {vtav.row[0].X:#0.000000000},{vtav.row[0].Y:#0.000000000},{vtav.row[0].Z:#0.000000000},{vtav.row[0].W:#0.000000000}");
            log?.WriteLine($"vtav[1]: {vtav.row[1].X:#0.000000000},{vtav.row[1].Y:#0.000000000},{vtav.row[1].Z:#0.000000000},{vtav.row[1].W:#0.000000000}");
            log?.WriteLine($"vtav[2]: {vtav.row[2].X:#0.000000000},{vtav.row[2].Y:#0.000000000},{vtav.row[2].Z:#0.000000000},{vtav.row[2].W:#0.000000000}");
            log?.WriteLine($"vtav[3]: {vtav.row[3].X:#0.000000000},{vtav.row[3].Y:#0.000000000},{vtav.row[3].Z:#0.000000000},{vtav.row[3].W:#0.000000000}");

            //log?.WriteLine($"v[0]: {v.M11:#0.000000000},{v.M12:#0.000000000},{v.M13:#0.000000000},{v.M14:#0.000000000}");
            //log?.WriteLine($"v[1]: {v.M21:#0.000000000},{v.M22:#0.000000000},{v.M23:#0.000000000},{v.M24:#0.000000000}");
            //log?.WriteLine($"v[2]: {v.M31:#0.000000000},{v.M32:#0.000000000},{v.M33:#0.000000000},{v.M34:#0.000000000}");
            //log?.WriteLine($"v[3]: {v.M41:#0.000000000},{v.M42:#0.000000000},{v.M43:#0.000000000},{v.M44:#0.000000000}");
            log?.WriteLine($"v:{v}");

            log?.WriteLine($"simd_u: {simd_u.X:#0.000000000},{simd_u.Y:#0.000000000},{simd_u.Z:#0.000000000},{simd_u.W:#0.000000000}");
            log?.WriteLine($"simd_v: {simd_v.X:#0.000000000},{simd_v.Y:#0.000000000},{simd_v.Z:#0.000000000},{simd_v.W:#0.000000000}");

            var simd_c = new Vector4(c);
            var simd_s = new Vector4(s);

            var x0 = simd_c * simd_u;
            var x1 = simd_s * simd_v;
            var x = x0 - x1;

            var y0 = simd_s * simd_u;
            var y1 = simd_c * simd_v;
            var y = y0 + y1;

            Item(ref v, 0, a, x.X);
            Item(ref v, 1, a, x.Y);
            Item(ref v, 2, a, x.Z);
            vtav[0, 3 - b] = x.W;

            Item(ref v, 0, b, y.X);
            Item(ref v, 1, b, y.Y);
            Item(ref v, 2, b, y.Z);
            vtav[1 - a, 2] = y.W;

            vtav[a, b] = 0f;
        }

        private void RotateQXY(ref Mat4x4 vtav, Vector4 c, Vector4 s, int a, int b)
        {
            log?.WriteLine("rotateq_xy");
            log?.WriteLine($"vtav[0]: {vtav.row[0].X:#0.000000000},{vtav.row[0].Y:#0.000000000},{vtav.row[0].Z:#0.000000000},{vtav.row[0].W:#0.000000000}");
            log?.WriteLine($"vtav[1]: {vtav.row[1].X:#0.000000000},{vtav.row[1].Y:#0.000000000},{vtav.row[1].Z:#0.000000000},{vtav.row[1].W:#0.000000000}");
            log?.WriteLine($"vtav[2]: {vtav.row[2].X:#0.000000000},{vtav.row[2].Y:#0.000000000},{vtav.row[2].Z:#0.000000000},{vtav.row[2].W:#0.000000000}");
            log?.WriteLine($"vtav[3]: {vtav.row[3].X:#0.000000000},{vtav.row[3].Y:#0.000000000},{vtav.row[3].Z:#0.000000000},{vtav.row[3].W:#0.000000000}");

            var u = new Vector4(vtav[a, a], vtav[a, a], vtav[a, a], 0f);
            var v = new Vector4(vtav[b, b], vtav[b, b], vtav[b, b], 0f);
            var A = new Vector4(vtav[a, b], vtav[a, b], vtav[a, b], 0f);

            log?.WriteLine($"u: {u.X:#0.000000000},{u.Y:#0.000000000},{u.Z:#0.000000000},{u.W:#0.000000000}");
            log?.WriteLine($"v: {v.X:#0.000000000},{v.Y:#0.000000000},{v.Z:#0.000000000},{v.W:#0.000000000}");
            log?.WriteLine($"A: {A.X:#0.000000000},{A.Y:#0.000000000},{A.Z:#0.000000000},{A.W:#0.000000000}");

            var twos = new Vector4(2f);
            var cc = c * c; // Sse.Multiply(c, c);
            var ss = s * s; // Sse.Multiply(s, s);

            // mx = 2.0 * c * s * A;
            var c2 = twos * c; //Sse.Multiply(twos, c);
            var c2s = c2 * s; //Sse.Multiply(c2, s);
            var mx = c2s * A; // Sse.Multiply(c2s, A);

            // x = cc * u - mx + ss * v;
            var x0 = cc * u; // Sse.Multiply(cc, u);
            var x1 = x0 - mx; // Sse.Subtract(x0, mx);
            var x2 = ss * v; // Sse.Multiply(ss, v);
            var x = x1 + x2;// Sse.Add(x1, x2);

            // y = ss * u + mx + cc * v;
            var y0 = ss * u; // Sse.Multiply(ss, u);
            var y1 = y0 + mx;
            var y2 = cc * v; // Sse.Multiply(cc, v);
            var y = y1 + y2; // Sse.Add(y1, y2);

            vtav[a, a] = x.X;
            vtav[b, b] = y.X;
        }

        private void GivensCoeffsSym(ref Vector4 c_result, ref Vector4 s_result, Mat4x4 vtav, int a, int b)
        {
            var simd_pp = new Vector4(vtav[a, a], vtav[a, a], vtav[a, a], 0f);
            var simd_pq = new Vector4(vtav[a, b], vtav[a, b], vtav[a, b], 0f);
            var simd_qq = new Vector4(vtav[b, b], vtav[b, b], vtav[b, b], 0f);

            // tau = (a_qq - a_pp) / (2.f * a_pq);
            var tau = (simd_qq - simd_pp) / (2f * simd_pq);

            // stt = sqrt(1.f + tau * tau);
            var stt = Vector4.SquareRoot(Vector4.One + tau * tau);

            // tan = 1.f / ((tau >= 0.f) ? (tau + stt) : (tau - stt));
            var tan = new Vector4(
                1f / ((tau.X >= 0f) ? (tau.X + stt.X) : (tau.X - stt.X)),
                1f / ((tau.Y >= 0f) ? (tau.Y + stt.Y) : (tau.Y - stt.Y)),
                1f / ((tau.Z >= 0f) ? (tau.Z + stt.Z) : (tau.Z - stt.Z)),
                1f / ((tau.W >= 0f) ? (tau.W + stt.W) : (tau.W - stt.W)));

            // c = rsqrt(1.f + tan * tan);
            var c = Vector4.One / Vector4.SquareRoot(Vector4.One + tan * tan);

            // s = tan * c;
            var s = tan * c;

            // if pq == 0.0: c = 1.f, s = 0.f
            c_result = new Vector4(simd_pq.X == 0f ? 1 : c.X, simd_pq.Y == 0f ? 1 : c.Y, simd_pq.Z == 0f ? 1 : c.Z, simd_pq.W == 0f ? 1 : c.W);
            s_result = new Vector4(simd_pq.X == 0f ? 0 : s.X, simd_pq.Y == 0f ? 0 : s.Y, simd_pq.Z == 0f ? 0 : s.Z, simd_pq.W == 0f ? 0 : s.W);

            log?.WriteLine($"c_result:{c_result.X:G9},{c_result.Y:G9},{c_result.Z:G9},{c_result.W:G9}");
            log?.WriteLine($"s_result:{s_result.X:G9},{s_result.Y:G9},{s_result.Z:G9},{s_result.W:G9}");
        }

        private static byte MM_Shuffle(byte fp3, byte fp2, byte fp1, byte fp0)
        {
            return ((byte)(((fp3) << 6) | ((fp2) << 4) | ((fp1) << 2) | ((fp0))));
        }

        private static void GenerateTriangles(VoxelInfo voxelInfo, ref MeshBuffer buffer)
        {
            foreach (var pair in voxelInfo.Edges)
            {
                var edge = pair.Key;
                var info = pair.Value;

                //var basePos = DecodeVoxelUniqueID(edge);
                var axis = (edge >> 30) & 0xff;

                var nodeId = edge & ~0xc0000000;
                var voxelIds = new[]
                {
                    nodeId - EncodedEdgeNodeOffsets[axis * 4 + 0],
                    nodeId - EncodedEdgeNodeOffsets[axis * 4 + 1],
                    nodeId - EncodedEdgeNodeOffsets[axis * 4 + 2],
                    nodeId - EncodedEdgeNodeOffsets[axis * 4 + 3],
                };

                var edgeVoxels = new int[4];
                var numFoundVoxels = 0;

                for (int i = 0; i < 4; i++)
                {
                    edgeVoxels[numFoundVoxels++] = voxelInfo.Indices[voxelIds[i]];
                }

                if (numFoundVoxels < 4)
                {
                    continue;
                }

                if (info.Winding == true)
                {
                    buffer.Triangles.Add(new MeshTriangle() { V0 = edgeVoxels[0], V1 = edgeVoxels[1], V2 = edgeVoxels[3] });
                    buffer.Triangles.Add(new MeshTriangle() { V0 = edgeVoxels[0], V1 = edgeVoxels[3], V2 = edgeVoxels[2] });
                }
                else
                {
                    buffer.Triangles.Add(new MeshTriangle() { V0 = edgeVoxels[0], V1 = edgeVoxels[3], V2 = edgeVoxels[1] });
                    buffer.Triangles.Add(new MeshTriangle() { V0 = edgeVoxels[0], V1 = edgeVoxels[2], V2 = edgeVoxels[3] });
                }
            }
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

        private VoxelInfo FindActiveVoxels(SuperPrimitiveConfig config, bool enableLog = false)
        {
            var info = new VoxelInfo();

            for (int x = 0; x < VoxelGridSize; x++)
            {
                for (int y = 0; y < VoxelGridSize; y++)
                {
                    for (int z = 0; z < VoxelGridSize; z++)
                    {
                        var idxPos = new ivec4(x, y, z, 0);
                        var p = new vec4(x - VoxelGridOffset, y - VoxelGridOffset, z - VoxelGridOffset, 1f);

                        for (uint axis = 0; axis < 3; axis++)
                        {
                            var q = p + AxisOffset[axis];

                            var pDensity = Density(config, p);
                            var qDensity = Density(config, q);

                            var zeroCrossing = pDensity >= 0f && qDensity < 0f || pDensity < 0f && qDensity >= 0f;

                            if (zeroCrossing == false)
                            {
                                continue;
                            }

                            var t = FindIntersection(config, p, q);
                            var pos = new vec4(glm.Mix(new vec3(p), new vec3(q), new vec3(t)), 1f);

                            if (enableLog == true && log != null)
                                log.WriteLine($"P:{pDensity:G9} Q:{qDensity:G9} T:{t:G9} POS:{pos.x:G9},{pos.y:G9},{pos.z:G9},{pos.w:G9}");

                            const float H = 0.001f;
                            var d0 = Density(config, pos + new vec4(H, 0f, 0f, 0f));
                            var d1 = Density(config, pos - new vec4(H, 0f, 0f, 0f));

                            var d2 = Density(config, pos + new vec4(0f, H, 0f, 0f));
                            var d3 = Density(config, pos - new vec4(0f, H, 0f, 0f));

                            var d4 = Density(config, pos + new vec4(0f, 0f, H, 0f));
                            var d5 = Density(config, pos - new vec4(0f, 0f, H, 0f));

                            if (enableLog == true && log != null)
                                log.WriteLine($"{d0.ToString("G9").ToLowerInvariant()},{d1.ToString("G9").ToLowerInvariant()},{d2.ToString("G9").ToLowerInvariant()},{d3.ToString("G9").ToLowerInvariant()},{d4.ToString("G9").ToLowerInvariant()},{d5.ToString("G9").ToLowerInvariant()}");

                            var n = new vec4(
                                d0 - d1,
                                d2 - d3,
                                d4 - d5,
                                //Density(config, pos + new vec4(H, 0f, 0f, 0f)) - Density(config, pos - new vec4(H, 0f, 0f, 0f)),
                                //Density(config, pos + new vec4(0f, H, 0f, 0f)) - Density(config, pos - new vec4(0f, H, 0f, 0f)),
                                //Density(config, pos + new vec4(0f, 0f, H, 0f)) - Density(config, pos - new vec4(0f, 0f, H, 0f)),
                                0f);
                            var l = (1f / (float)Math.Sqrt(vec4.Dot(n, n)));
                            var normal = n * l;

                            if (enableLog == true && log != null)
                            {
                                log.WriteLine($"L:{n.Length:g9}");
                                log.WriteLine($"N:{n.x:g9},{n.y:g9},{n.z:g9},{n.w:g9}");
                                log.WriteLine($"N:{normal.x:g9},{normal.y:g9},{normal.z:g9},{normal.w:g9}");
                            }

                            var edgeInfo = new EdgeInfo() { Position = pos, Normal = normal, Winding = pDensity >= 0f };

                            var code = EncodeAxisUniqueID(axis, x, y, z);
                            info.Edges.Add(code, edgeInfo);

                            if (enableLog == true && log != null)
                                log.WriteLine($"Edge:{code}");

                            var edgeNodes = EdgeNodeOffsets[axis];

                            for (int i = 0; i < 4; i++)
                            {
                                var nodeIdxPos = idxPos - edgeNodes[i];
                                var nodeID = EncodeVoxelUniqueID(nodeIdxPos);
                                info.Voxels.Add(nodeID);
                                if (enableLog == true && log != null)
                                    log.WriteLine($"Node:{nodeID}");
                            }
                        }
                    }
                }
            }

            return info;
        }

        //private float FindIntersection(float d0, float d1)
        //{
        //    var t = -d0 / (d1 - d0);

        //    return t;
        //}

        private uint EncodeVoxelUniqueID(ivec4 idxPos, bool logEnabled = false)
        {
            //unchecked
            {
                var i = idxPos.x | (idxPos.y << 10) | (idxPos.z << 20);

                if (logEnabled == true && log != null)
                    log.WriteLine($"{idxPos.x},{idxPos.y},{idxPos.z},{idxPos.w} -> {i}");

                return (uint)i;
            }
        }

        private static uint EncodeAxisUniqueID(uint axis, int x, int y, int z)
        {
            return (uint)x | ((uint)y << 10) | ((uint)z << 20) | (axis << 30);
        }

        private static float FindIntersection(SuperPrimitiveConfig config, vec4 p0, vec4 p1)
        {
            const int FindEdgeInfoSteps = 16;
            const float FindEdgeInfoIncrement = 1f / FindEdgeInfoSteps;

            var minValue = float.MaxValue;
            var currentT = 0f;
            var t = 0f;

            for (int i = 0; i < FindEdgeInfoSteps; i++)
            {
                var p = glm.Mix(p0, p1, new vec4(currentT));
                var d = glm.Abs(Density(config, p));
                if (d < minValue)
                {
                    t = currentT;
                    minValue = d;
                }

                currentT += FindEdgeInfoIncrement;
            }

            return t;
        }

        private static float Density(SuperPrimitiveConfig config, vec4 p)
        {
            return SuperPrimitive(new vec3(p) / scale, new vec4(config.S), new vec2(config.R)) * scale;
        }

        // The "super primitve" -- use the parameters to configure different shapes from a single function
        // see https://www.shadertoy.com/view/MsVGWG

        static float SuperPrimitive(vec3 p, vec4 s, vec2 r)
        {
            var d = glm.Abs(p) - new vec3(s);

            float q = glm.Length(new vec2(glm.Max(d.x + r.x, 0f), glm.Max(d.y + r.x, 0f)));
            q += glm.Min(-r.x, glm.Max(d.x, d.y));
            q = (glm.Abs((q + s.w)) - s.w);

            return glm.Length(new vec2(glm.Max(q + r.y, 0f), glm.Max(d.z + r.y, 0f))) + glm.Min(-r.y, glm.Max(q, d.z));
        }

        private static SuperPrimitiveConfig ConfigForShape(Primitive primitive)
        {
            var config = new SuperPrimitiveConfig();

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

        //struct rubbish
        //{
        //    public float x;
        //    public float y;
        //    public float z;

        //    public rubbish(float x, float y, float z)
        //    {
        //        this.x = x;
        //        this.y = y;
        //        this.z = z;
        //    }
        //}

        //float[] g_vertex_buffer_data = new[]
        //    {
        //       -1.0f, -1.0f, 0.0f,
        //       0f, 0f, -1f,
        //       1f, 0f, 0f,
        //       1.0f, -1.0f, 0.0f,
        //       0f, 0f, -1f,
        //       0f, 1f, 0f,
        //       0.0f,  1.0f, 0.0f,
        //       0f, 0f, -1f,
        //       0f, 0f, 1f
        //    };

        //uint[] indices = new uint[]
        //{
        //    0, 1, 2
        //};

        //private void TestFrame(Mesh mesh)
        //{
        //    var VertexArrayID = GL.GenVertexArray();
        //    GL.BindVertexArray(VertexArrayID);

        //    // Generate 1 buffer, put the resulting identifier in vertexbuffer
        //    var vertexbuffer = GL.GenBuffer();
        //    // The following commands will talk about our 'vertexbuffer' buffer
        //    GL.BindBuffer(BufferTargetARB.ArrayBuffer, vertexbuffer);
        //    // Give our vertices to OpenGL.
        //    //GL.BufferData(BufferTargetARB.ArrayBuffer, g_vertex_buffer_data, BufferUsageARB.StaticDraw);
        //    GL.BufferData(BufferTargetARB.ArrayBuffer, meshBuffer.Vertices.ToArray(), BufferUsageARB.StaticDraw);

        //    // 1st attribute buffer : vertices
        //    GL.EnableVertexAttribArray(0);
        //    GL.EnableVertexAttribArray(1);
        //    GL.EnableVertexAttribArray(2);

        //    GL.BindBuffer(BufferTargetARB.ArrayBuffer, vertexbuffer);
        //    //GL.BindBuffer(BufferTargetARB.ArrayBuffer, mesh.VertexBuffer);

        //    GL.VertexAttribPointer(
        //       0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
        //       3,                  // size
        //       VertexAttribPointerType.Float,           // type
        //       false,           // normalized?
        //       36,                  // stride
        //       0            // array buffer offset
        //    );

        //    GL.VertexAttribPointer(
        //       1,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
        //       3,                  // size
        //       VertexAttribPointerType.Float,           // type
        //       false,           // normalized?
        //       36,                  // stride
        //       12            // array buffer offset
        //    );

        //    GL.VertexAttribPointer(
        //       2,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
        //       3,                  // size
        //       VertexAttribPointerType.Float,           // type
        //       false,           // normalized?
        //       36,                  // stride
        //       24            // array buffer offset
        //    );

        //    var elementBuffer = GL.GenBuffer();
        //    GL.BindBuffer(BufferTargetARB.ElementArrayBuffer, elementBuffer);
        //    GL.BufferData(BufferTargetARB.ElementArrayBuffer, indices, BufferUsageARB.StaticDraw);

        //    var projection = mat4.Perspective(glm.Radians(45f), 1920f / 1080f, 0.1f, 500f);
        //    var view = mat4.LookAt(new vec3(4, 3, 3), vec3.Zero, vec3.UnitY);
        //    var model = mat4.Identity;
        //    var mvp = projection * view * model;

        //    GL.UseProgram(this.shader);

        //    var mat = GL.GetUniformLocation(this.shader, "MVP");
        //    CheckGLError();

        //    var m4 = new OpenTK.Mathematics.Matrix4(mvp.m00, mvp.m01, mvp.m02, mvp.m03, mvp.m10, mvp.m11, mvp.m12, mvp.m13, mvp.m20, mvp.m21, mvp.m22, mvp.m23, mvp.m30, mvp.m31, mvp.m32, mvp.m33);
        //    GL.UniformMatrix4f(mat, false, m4);
        //    CheckGLError();

        //    var useColour = GL.GetUniformLocation(this.shader, "useUniformColour");
        //    CheckGLError();
        //    GL.Uniform1i(useColour, 0);
        //    CheckGLError();

        //    //GL.BindVertexArray(mesh.VertexArrayObj);
        //    //GL.DrawElements(PrimitiveType.Triangles, mesh.NumIndices, DrawElementsType.UnsignedInt, 0);

        //    // Draw the triangle !
        //    //GL.DrawArrays(PrimitiveType.Triangles, 0, 3); // Starting from vertex 0; 3 vertices total -> 1 triangle
        //    GL.DrawElements(PrimitiveType.Triangles, indices.Length, DrawElementsType.UnsignedInt, 0);

        //    GL.DisableVertexAttribArray(2);
        //    GL.DisableVertexAttribArray(1);
        //    GL.DisableVertexAttribArray(0);

        //    GL.DeleteVertexArray(VertexArrayID);
        //    GL.DeleteBuffer(elementBuffer);
        //}

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
            }

            if (ImGui.Button("Refresh"))
            {
                viewerOpts.RefreshModel = true;
            }

            ImGui.End();
        }
    }
}
