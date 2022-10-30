//#pragma warning disable RECS0018 // Comparison of floating point numbers with equality operator

using FastDC;

using GlmSharp;

using System;
using System.Numerics;

namespace Facet42
{
    internal static class Qef
    {
        public static float SolveFromPoints4d(vec4[] positions, vec4[] normals, int count, out vec3 position)
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

            var error = SolveFromPoints(p, n, count, out var solved);
            position = new vec3(solved.X, solved.Y, solved.Z);
            return error;
        }

        private static float SolveFromPoints(Vector4[] positions, Vector4[] normals, int count, out Vector4 solvedPosition)
        {
            var pointAccum = new Vector4();
            var ATb = new Vector4();

            var ATA = new Mat4x4();

            for (int i = 0; i < count && i < 12; i++)
            {
                SimdAdd(positions[i], normals[i], ref ATA, ref ATb, ref pointAccum);
            }

            //log?.WriteLine($"ATb:{ATb.X:0.000000000},{ATb.Y:0.000000000},{ATb.Z:0.000000000},{ATb.W:0.000000000}");

            return SimdSolve(ATA, ATb, pointAccum, out solvedPosition);
        }

        private static void SimdAdd(Vector4 p, Vector4 n, ref Mat4x4 ATA, ref Vector4 ATb, ref Vector4 pointAccum)
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

        private static float SimdSolve(Mat4x4 ATA, Vector4 ATb, Vector4 pointAccum, out Vector4 solved)
        {
            var massPoint = pointAccum / pointAccum.W;

            //log?.WriteLine($"pointAccum: {pointAccum.X:0.000000000},{pointAccum.Y:0.000000000},{pointAccum.Z:0.000000000},{pointAccum.W:0.000000000}");
            //log?.WriteLine($"MassPoint: {massPoint.X:0.000000000},{massPoint.Y:0.000000000},{massPoint.Z:0.000000000},{massPoint.W:0.000000000}");

            var p = ATb - Vector4.Transform(massPoint, ATA.AsMatrix4x4());

            Svd.SolveATAATb(ATA, p, out solved);

            var error = SimdCalcError(ATA, solved, ATb);

            solved += massPoint;

            //log?.WriteLine($"Solved: {solved.X:0.000000000},{solved.Y:0.000000000},{solved.Z:0.000000000},{solved.W:0.000000000}");

            return error;
        }

        private static float SimdCalcError(Mat4x4 A, Vector4 x, Vector4 b)
        {
            var tmp = b - Vector4.Transform(x, A.AsMatrix4x4());
            return Vector4.Dot(tmp, tmp);
        }
    }
}
