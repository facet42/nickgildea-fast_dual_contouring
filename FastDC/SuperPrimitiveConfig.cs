namespace FastDC
{
    using GlmSharp;

    public enum Primitive
    {
        Cube,
        Cylinder,
        Pill,
        Corridor,
        Torus
    }

    internal class SuperPrimitiveConfig
    {
        public vec4 S;
        public vec2 R;
    }
}
