namespace FastDC
{
    using OpenTK.Graphics;
    using OpenTK.Graphics.OpenGL;

    using System;
    using System.Collections.Generic;
    using System.Diagnostics;
    using System.Linq;
    using System.Text;
    using System.Threading.Tasks;

    internal class Graphics
    {
        public static ProgramHandle CreateProgram(string name, string vertexSource, string fragmentSource)
        {
            var program = GL.CreateProgram();

            var vertex = CompileShader(name, ShaderType.VertexShader, vertexSource);
            var fragment = CompileShader(name, ShaderType.FragmentShader, fragmentSource);

            GL.AttachShader(program, vertex);
            GL.AttachShader(program, fragment);

            GL.LinkProgram(program);

            //int success = 0;
            //GL.GetProgrami(program, ProgramPropertyARB.LinkStatus, ref success);
            //if (success == 0)
            {
                GL.GetProgramInfoLog(program, out var info);
                if (string.IsNullOrWhiteSpace(info) == false)
                {
                    Console.WriteLine($"GL.LinkProgram had info log [{name}]:\n{info}");
                    Debug.WriteLine($"GL.LinkProgram had info log [{name}]:\n{info}");
                }
            }

            GL.DetachShader(program, vertex);
            GL.DetachShader(program, fragment);

            GL.DeleteShader(vertex);
            GL.DeleteShader(fragment);

            return program;
        }

        public static void DeleteProgram(ProgramHandle shader)
        {
            GL.DeleteProgram(shader);
        }

        private static ShaderHandle CompileShader(string name, ShaderType shaderType, string source)
        {
            var shader = GL.CreateShader(shaderType);

            GL.ShaderSource(shader, source);
            GL.CompileShader(shader);

            //int success = 0;
            //GL.GetShaderi(shader, ShaderParameterName.CompileStatus, ref success);
            //if (success == 0)
            {
                GL.GetShaderInfoLog(shader, out var info);
                if (string.IsNullOrWhiteSpace(info) == false)
                {
                    Console.WriteLine($"GL.CompileShader for shader '{name}' [{shaderType}] had info log:\n{info}");
                    Debug.WriteLine($"GL.CompileShader for shader '{name}' [{shaderType}] had info log:\n{info}");
                }
            }

            return shader;
        }
    }
}
