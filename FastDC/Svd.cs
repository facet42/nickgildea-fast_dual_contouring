namespace Facet42
{
    using FastDC;

    using System.Numerics;

    internal static class Svd
    {
        const int SvdNumSweeps = 5;
        const float PseudoInverseThreshold = 0.001f;

        public static void SolveATAATb(Mat4x4 ATA, Vector4 ATb, out Vector4 x)
        {
            var V = new Mat4x4(
                1f, 0, 0, 0,
                0, 1f, 0, 0,
                0, 0, 1f, 0,
                0, 0, 0, 0
                );

            var sigma = SolveSym(ref V, ATA);

            //log?.WriteLine($"Sigma:{ToString(sigma)}");

            //log?.WriteLine($"V[0]:{V[0, 0]:#0.000000000},{V[0, 1]:#0.000000000},{V[0, 2]:#0.000000000},{V[0, 3]:#0.000000000}");
            //log?.WriteLine($"V[1]:{V[1, 0]:#0.000000000},{V[1, 1]:#0.000000000},{V[1, 2]:#0.000000000},{V[1, 3]:#0.000000000}");
            //log?.WriteLine($"V[2]:{V[2, 0]:#0.000000000},{V[2, 1]:#0.000000000},{V[2, 2]:#0.000000000},{V[2, 3]:#0.000000000}");
            //log?.WriteLine($"V[3]:{V[3, 0]:#0.000000000},{V[3, 1]:#0.000000000},{V[3, 2]:#0.000000000},{V[3, 3]:#0.000000000}");

            PseudoInverse(out var Vinv, sigma, V);

            //log?.WriteLine($"Vinv[0]:{Vinv.M11:#0.000000000},{Vinv.M12:#0.000000000},{Vinv.M13:#0.000000000},{Vinv.M14:#0.000000000}");
            //log?.WriteLine($"Vinv[1]:{Vinv.M21:#0.000000000},{Vinv.M22:#0.000000000},{Vinv.M23:#0.000000000},{Vinv.M24:#0.000000000}");
            //log?.WriteLine($"Vinv[2]:{Vinv.M31:#0.000000000},{Vinv.M32:#0.000000000},{Vinv.M33:#0.000000000},{Vinv.M34:#0.000000000}");
            //log?.WriteLine($"Vinv[3]:{Vinv.M41:#0.000000000},{Vinv.M42:#0.000000000},{Vinv.M43:#0.000000000},{Vinv.M44:#0.000000000}");
            //log?.WriteLine($"Vinv:{Vinv}");

            x = Vector4.Transform(ATb, Vinv.AsMatrix4x4());
        }

        private static Vector4 SolveSym(ref Mat4x4 v, Mat4x4 a)
        {
            var vtav = a;

            //log?.WriteLine("svd_solve_sym");
            //log?.WriteLine($"vtav[0]: {vtav.row[0].X:#0.000000000},{vtav.row[0].Y:#0.000000000},{vtav.row[0].Z:#0.000000000},{vtav.row[0].W:#0.000000000}");
            //log?.WriteLine($"vtav[1]: {vtav.row[1].X:#0.000000000},{vtav.row[1].Y:#0.000000000},{vtav.row[1].Z:#0.000000000},{vtav.row[1].W:#0.000000000}");
            //log?.WriteLine($"vtav[2]: {vtav.row[2].X:#0.000000000},{vtav.row[2].Y:#0.000000000},{vtav.row[2].Z:#0.000000000},{vtav.row[2].W:#0.000000000}");
            //log?.WriteLine($"vtav[3]: {vtav.row[3].X:#0.000000000},{vtav.row[3].Y:#0.000000000},{vtav.row[3].Z:#0.000000000},{vtav.row[3].W:#0.000000000}");

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

                //log?.WriteLine($"c:{c}");
                //log?.WriteLine($"s:{s}");
            }

            //log?.WriteLine($"vtav[0]: {vtav.row[0].X:#0.000000000},{vtav.row[0].Y:#0.000000000},{vtav.row[0].Z:#0.000000000},{vtav.row[0].W:#0.000000000}");
            //log?.WriteLine($"vtav[1]: {vtav.row[1].X:#0.000000000},{vtav.row[1].Y:#0.000000000},{vtav.row[1].Z:#0.000000000},{vtav.row[1].W:#0.000000000}");
            //log?.WriteLine($"vtav[2]: {vtav.row[2].X:#0.000000000},{vtav.row[2].Y:#0.000000000},{vtav.row[2].Z:#0.000000000},{vtav.row[2].W:#0.000000000}");
            //log?.WriteLine($"vtav[3]: {vtav.row[3].X:#0.000000000},{vtav.row[3].Y:#0.000000000},{vtav.row[3].Z:#0.000000000},{vtav.row[3].W:#0.000000000}");

            return new Vector4(vtav.row[0].X, vtav.row[1].Y, vtav.row[2].Z, 0f);
        }

        private static void PseudoInverse(out Mat4x4 o, Vector4 sigma, Mat4x4 v)
        {
            var invDet = InvDet(sigma);
            //log?.WriteLine($"InvDet:{invDet}");

            //log?.WriteLine($"v[0]:{v.M11:#0.000000000},{v.M12:#0.000000000},{v.M13:#0.000000000},{v.M14:#0.000000000}");
            //log?.WriteLine($"v[1]:{v.M21:#0.000000000},{v.M22:#0.000000000},{v.M23:#0.000000000},{v.M24:#0.000000000}");
            //log?.WriteLine($"v[2]:{v.M31:#0.000000000},{v.M32:#0.000000000},{v.M33:#0.000000000},{v.M34:#0.000000000}");
            //log?.WriteLine($"v[3]:{v.M41:#0.000000000},{v.M42:#0.000000000},{v.M43:#0.000000000},{v.M44:#0.000000000}");
            //log?.WriteLine($"v:{v}");

            var m = new Mat4x4();
            m.row[0] = v.row[0] * invDet;
            m.row[1] = v.row[1] * invDet;
            m.row[2] = v.row[2] * invDet;
            m.row[3] = new Vector4();

            //log?.WriteLine($"m[0]:{ToString(m.row[0])}");
            //log?.WriteLine($"m[1]:{ToString(m.row[1])}");
            //log?.WriteLine($"m[2]:{ToString(m.row[2])}");
            //log?.WriteLine($"m[3]:{ToString(m.row[3])}");

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
            //log?.WriteLine($"o:{o}");
        }

        private static Vector4 InvDet(Vector4 v)
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

        private static void GivensCoeffsSym(ref Vector4 c_result, ref Vector4 s_result, Mat4x4 vtav, int a, int b)
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

            //log?.WriteLine($"c_result:{c_result.X:G9},{c_result.Y:G9},{c_result.Z:G9},{c_result.W:G9}");
            //log?.WriteLine($"s_result:{s_result.X:G9},{s_result.Y:G9},{s_result.Z:G9},{s_result.W:G9}");
        }

        private static void RotateQXY(ref Mat4x4 vtav, Vector4 c, Vector4 s, int a, int b)
        {
            //log?.WriteLine("rotateq_xy");
            //log?.WriteLine($"vtav[0]: {vtav.row[0].X:#0.000000000},{vtav.row[0].Y:#0.000000000},{vtav.row[0].Z:#0.000000000},{vtav.row[0].W:#0.000000000}");
            //log?.WriteLine($"vtav[1]: {vtav.row[1].X:#0.000000000},{vtav.row[1].Y:#0.000000000},{vtav.row[1].Z:#0.000000000},{vtav.row[1].W:#0.000000000}");
            //log?.WriteLine($"vtav[2]: {vtav.row[2].X:#0.000000000},{vtav.row[2].Y:#0.000000000},{vtav.row[2].Z:#0.000000000},{vtav.row[2].W:#0.000000000}");
            //log?.WriteLine($"vtav[3]: {vtav.row[3].X:#0.000000000},{vtav.row[3].Y:#0.000000000},{vtav.row[3].Z:#0.000000000},{vtav.row[3].W:#0.000000000}");

            var u = new Vector4(vtav[a, a], vtav[a, a], vtav[a, a], 0f);
            var v = new Vector4(vtav[b, b], vtav[b, b], vtav[b, b], 0f);
            var A = new Vector4(vtav[a, b], vtav[a, b], vtav[a, b], 0f);

            //log?.WriteLine($"u: {u.X:#0.000000000},{u.Y:#0.000000000},{u.Z:#0.000000000},{u.W:#0.000000000}");
            //log?.WriteLine($"v: {v.X:#0.000000000},{v.Y:#0.000000000},{v.Z:#0.000000000},{v.W:#0.000000000}");
            //log?.WriteLine($"A: {A.X:#0.000000000},{A.Y:#0.000000000},{A.Z:#0.000000000},{A.W:#0.000000000}");

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

        private static void RotateXY(ref Mat4x4 vtav, ref Mat4x4 v, float c, float s, int a, int b)
        {
            var simd_u = new Vector4(v[0, a], v[1, a], v[2, a], vtav[0, 3 - b]);
            var simd_v = new Vector4(v[0, b], v[1, b], v[2, b], vtav[1 - a, 2]);

            //log?.WriteLine("rotate_xy");
            //log?.WriteLine($"a:{a} b:{b}");
            //log?.WriteLine($"vtav[0]: {vtav.row[0].X:#0.000000000},{vtav.row[0].Y:#0.000000000},{vtav.row[0].Z:#0.000000000},{vtav.row[0].W:#0.000000000}");
            //log?.WriteLine($"vtav[1]: {vtav.row[1].X:#0.000000000},{vtav.row[1].Y:#0.000000000},{vtav.row[1].Z:#0.000000000},{vtav.row[1].W:#0.000000000}");
            //log?.WriteLine($"vtav[2]: {vtav.row[2].X:#0.000000000},{vtav.row[2].Y:#0.000000000},{vtav.row[2].Z:#0.000000000},{vtav.row[2].W:#0.000000000}");
            //log?.WriteLine($"vtav[3]: {vtav.row[3].X:#0.000000000},{vtav.row[3].Y:#0.000000000},{vtav.row[3].Z:#0.000000000},{vtav.row[3].W:#0.000000000}");

            //log?.WriteLine($"v[0]: {v.M11:#0.000000000},{v.M12:#0.000000000},{v.M13:#0.000000000},{v.M14:#0.000000000}");
            //log?.WriteLine($"v[1]: {v.M21:#0.000000000},{v.M22:#0.000000000},{v.M23:#0.000000000},{v.M24:#0.000000000}");
            //log?.WriteLine($"v[2]: {v.M31:#0.000000000},{v.M32:#0.000000000},{v.M33:#0.000000000},{v.M34:#0.000000000}");
            //log?.WriteLine($"v[3]: {v.M41:#0.000000000},{v.M42:#0.000000000},{v.M43:#0.000000000},{v.M44:#0.000000000}");
            //log?.WriteLine($"v:{v}");

            //log?.WriteLine($"simd_u: {simd_u.X:#0.000000000},{simd_u.Y:#0.000000000},{simd_u.Z:#0.000000000},{simd_u.W:#0.000000000}");
            //log?.WriteLine($"simd_v: {simd_v.X:#0.000000000},{simd_v.Y:#0.000000000},{simd_v.Z:#0.000000000},{simd_v.W:#0.000000000}");

            var simd_c = new Vector4(c);
            var simd_s = new Vector4(s);

            var x0 = simd_c * simd_u;
            var x1 = simd_s * simd_v;
            var x = x0 - x1;

            var y0 = simd_s * simd_u;
            var y1 = simd_c * simd_v;
            var y = y0 + y1;

            v[0, a] = x.X;
            v[1, a] = x.Y;
            v[2, a] = x.Z;
            vtav[0, 3 - b] = x.W;

            v[0, b] = y.X;
            v[1, b] = y.Y;
            v[2, b] = y.Z;
            vtav[1 - a, 2] = y.W;

            vtav[a, b] = 0f;
        }
    }
}
