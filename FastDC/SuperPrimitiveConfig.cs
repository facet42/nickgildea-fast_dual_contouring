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
        public float Scale;

        public float Density(vec4 p)
        {
            return SuperPrimitive(new vec3(p) / this.Scale, new vec4(this.S), new vec2(this.R)) * this.Scale;
        }

        // The "super primitve" -- use the parameters to configure different shapes from a single function
        // see https://www.shadertoy.com/view/MsVGWG

        private float SuperPrimitive(vec3 p, vec4 s, vec2 r)
        {
            var d = glm.Abs(p) - new vec3(s);

            float q = glm.Length(new vec2(glm.Max(d.x + r.x, 0f), glm.Max(d.y + r.x, 0f)));
            q += glm.Min(-r.x, glm.Max(d.x, d.y));
            q = (glm.Abs((q + s.w)) - s.w);

            return glm.Length(new vec2(glm.Max(q + r.y, 0f), glm.Max(d.z + r.y, 0f))) + glm.Min(-r.y, glm.Max(q, d.z));
        }
    }
}
