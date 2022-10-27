namespace FastDC
{
    using System;
    using System.Numerics;
    using System.Text;

    public struct Mat4x4
    {
        public Vector4[] row = new Vector4[4];

        public Mat4x4()
        {
        }

        public Mat4x4(
            float m00, float m01, float m02, float m03,
            float m10, float m11, float m12, float m13,
            float m20, float m21, float m22, float m23,
            float m30, float m31, float m32, float m33)
        {
            this[0, 0] = m00;
            this[0, 1] = m01;
            this[0, 2] = m02;
            this[0, 3] = m03;

            this[1, 0] = m10;
            this[1, 1] = m11;
            this[1, 2] = m12;
            this[1, 3] = m13;

            this[2, 0] = m20;
            this[2, 1] = m21;
            this[2, 2] = m22;
            this[2, 3] = m23;

            this[3, 0] = m30;
            this[3, 1] = m31;
            this[3, 2] = m32;
            this[3, 3] = m33;
        }

        public Matrix4x4 AsMatrix4x4()
        {
            return new Matrix4x4(
                row[0].X, row[0].Y, row[0].Z, row[0].W,
                row[1].X, row[1].Y, row[1].Z, row[1].W,
                row[2].X, row[2].Y, row[2].Z, row[2].W,
                row[3].X, row[3].Y, row[3].Z, row[3].W);
        }

        public float this[int row, int column]
        {
            get
            {
                switch (column)
                {
                    case 0:
                        return this.row[row].X;

                    case 1:
                        return this.row[row].Y;

                    case 2:
                        return this.row[row].Z;

                    case 3:
                        return this.row[row].W;

                    default:
                        throw new IndexOutOfRangeException();
                }
            }

            set
            {
                switch (column)
                {
                    case 0:
                        this.row[row].X = value;
                        return;

                    case 1:
                        this.row[row].Y = value;
                        return;

                    case 2:
                        this.row[row].Z = value;
                        return;

                    case 3:
                        this.row[row].W = value;
                        return;

                    default:
                        throw new IndexOutOfRangeException();
                }
            }
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            sb.AppendLine($"{this[0, 0]:#0.000000000},{this[0, 1]:#0.000000000},{this[0, 2]:#0.000000000},{this[0, 3]:#0.000000000}");
            sb.AppendLine($"{this[1, 0]:#0.000000000},{this[1, 1]:#0.000000000},{this[1, 2]:#0.000000000},{this[1, 3]:#0.000000000}");
            sb.AppendLine($"{this[2, 0]:#0.000000000},{this[2, 1]:#0.000000000},{this[2, 2]:#0.000000000},{this[2, 3]:#0.000000000}");
            sb.Append($"{this[3, 0]:#0.000000000},{this[3, 1]:#0.000000000},{this[3, 2]:#0.000000000},{this[3, 3]:#0.000000000}");

            return sb.ToString();
        }
    }
}
