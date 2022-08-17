namespace FastDC
{
    using ImGuiNET;

    using OpenTK.Graphics;
    using OpenTK.Graphics.OpenGL;
    using OpenTK.Mathematics;
    using OpenTK.Windowing.GraphicsLibraryFramework;

    using System;
    using System.Diagnostics;
    using System.Runtime.CompilerServices;
    using System.Xml.Linq;

    internal class ImGuiController
    {
        private int windowWidth;
        private int windowHeight;
        private bool frameBegun;
        private int vertexBufferSize;
        private int indexBufferSize;
        private VertexArrayHandle vertexArray;
        private BufferHandle vertexBuffer;
        private BufferHandle indexBuffer;
        private ProgramHandle shader;
        private int shaderProjectionMatrixLoc;
        private int shaderFontTextureLoc;
        private TextureHandle fontTexture;
        private System.Numerics.Vector2 scaleFactor = new System.Numerics.Vector2(1f, 1f);

        readonly List<char> PressedChars = new();

        public ImGuiController(int x, int y)
        {
            this.windowWidth = x;
            this.windowHeight = y;

            var context = ImGui.CreateContext();
            ImGui.SetCurrentContext(context);

            var io = ImGui.GetIO();
            io.Fonts.AddFontDefault();

            io.BackendFlags |= ImGuiBackendFlags.RendererHasVtxOffset;

            CreateDeviceResources();
            SetKeyMappings();

            SetPerFrameImGuiData(1f / 60f);

            ImGui.NewFrame();
            this.frameBegun = true;
        }

        private void CreateDeviceResources()
        {
            this.vertexBufferSize = 10000;
            this.indexBufferSize = 2000;

            int prevVAO = 0;
            GL.GetInteger(GetPName.VertexArrayBinding, ref prevVAO);

            int prevArrayBuffer = 0;
            GL.GetInteger(GetPName.ArrayBufferBinding, ref prevArrayBuffer);

            this.vertexArray = GL.GenVertexArray();
            GL.BindVertexArray(this.vertexArray);

            this.vertexBuffer = GL.GenBuffer();
            GL.BindBuffer(BufferTargetARB.ArrayBuffer, this.vertexBuffer);
            GL.BufferData(BufferTargetARB.ArrayBuffer, vertexBufferSize, IntPtr.Zero, BufferUsageARB.DynamicDraw);

            this.indexBuffer = GL.GenBuffer();
            GL.BindBuffer(BufferTargetARB.ElementArrayBuffer, this.indexBuffer);
            GL.BufferData(BufferTargetARB.ElementArrayBuffer, indexBufferSize, IntPtr.Zero, BufferUsageARB.DynamicDraw);

            RecreateFontDeviceTexture();

            var VertexSource = @"#version 330 core
                uniform mat4 projection_matrix;
                layout(location = 0) in vec2 in_position;
                layout(location = 1) in vec2 in_texCoord;
                layout(location = 2) in vec4 in_color;
                out vec4 color;
                out vec2 texCoord;
                void main()
                {
                    gl_Position = projection_matrix * vec4(in_position, 0, 1);
                    color = in_color;
                    texCoord = in_texCoord;
                }";
            var FragmentSource = @"#version 330 core
                uniform sampler2D in_fontTexture;
                in vec4 color;
                in vec2 texCoord;
                out vec4 outputColor;
                void main()
                {
                    //outputColor = vec4(1, 0, 0, 1);
                    outputColor = color * texture(in_fontTexture, texCoord);
                }";

            this.shader = CreateProgram("ImGui", VertexSource, FragmentSource);
            this.shaderProjectionMatrixLoc = GL.GetUniformLocation(this.shader, "projection_matrix");
            this.shaderFontTextureLoc = GL.GetUniformLocation(this.shader, "in_fontTexture");

            var stride = Unsafe.SizeOf<ImDrawVert>();
            GL.VertexAttribPointer(0, 2, VertexAttribPointerType.Float, false, stride, 0);
            GL.VertexAttribPointer(1, 2, VertexAttribPointerType.Float, false, stride, 8);
            GL.VertexAttribPointer(2, 4, VertexAttribPointerType.UnsignedByte, true, stride, 16);

            GL.EnableVertexAttribArray(0);
            GL.EnableVertexAttribArray(1);
            GL.EnableVertexAttribArray(2);

            GL.BindVertexArray((VertexArrayHandle)prevVAO);
            GL.BindBuffer(BufferTargetARB.ArrayBuffer, (BufferHandle)prevArrayBuffer);
        }

        private void RecreateFontDeviceTexture()
        {
            var io = ImGui.GetIO();
            io.Fonts.GetTexDataAsRGBA32(out IntPtr pixels, out var width, out var height);

            var mips = (int)Math.Floor(Math.Log(Math.Max(width, height), 2));

            int prevActiveTexture = 0;
            GL.GetInteger(GetPName.ActiveTexture, ref prevActiveTexture);
            GL.ActiveTexture(TextureUnit.Texture0);

            int prevTexture2D = 0;
            GL.GetInteger(GetPName.TextureBinding2d, ref prevTexture2D);

            this.fontTexture = GL.GenTexture();
            GL.BindTexture(TextureTarget.Texture2d, this.fontTexture);
            GL.TexStorage2D(TextureTarget.Texture2d, mips, SizedInternalFormat.Rgba8, width, height);

            GL.TexSubImage2D(TextureTarget.Texture2d, 0, 0, 0, width, height, PixelFormat.Bgra, PixelType.UnsignedByte, pixels);

            GL.GenerateMipmap(TextureTarget.Texture2d);

            GL.TexParameteri(TextureTarget.Texture2d, TextureParameterName.TextureWrapS, (int)TextureWrapMode.Repeat);
            GL.TexParameteri(TextureTarget.Texture2d, TextureParameterName.TextureWrapT, (int)TextureWrapMode.Repeat);

            GL.TexParameteri(TextureTarget.Texture2d, TextureParameterName.TextureMaxLevel, mips - 1);

            GL.TexParameteri(TextureTarget.Texture2d, TextureParameterName.TextureMagFilter, (int)TextureMagFilter.Linear);
            GL.TexParameteri(TextureTarget.Texture2d, TextureParameterName.TextureMinFilter, (int)TextureMagFilter.Linear);

            GL.BindTexture(TextureTarget.Texture2d, (TextureHandle)prevTexture2D);
            GL.ActiveTexture((TextureUnit)prevActiveTexture);

            io.Fonts.SetTexID((IntPtr)this.fontTexture.Handle);
            io.Fonts.ClearTexData();
        }

        internal void Render()
        {
            if (this.frameBegun)
            {
                this.frameBegun = false;
                ImGui.Render();

                var data = ImGui.GetDrawData();
                RenderImDrawData(data);
            }
        }

        internal void Update(AppWindow window, double time)
        {
            if (this.frameBegun)
            {
                ImGui.Render();
            }

            SetPerFrameImGuiData(time);
            UpdateImGuiInput(window);

            this.frameBegun = true;
            ImGui.NewFrame();
        }

        private void SetPerFrameImGuiData(double time)
        {
            var io = ImGui.GetIO();
            io.DisplaySize = new System.Numerics.Vector2(
                this.windowWidth / this.scaleFactor.X,
                this.windowHeight / this.scaleFactor.Y);
            io.DisplayFramebufferScale = this.scaleFactor;
            io.DeltaTime = (float)time;
        }

        private void UpdateImGuiInput(AppWindow window)
        {
            var io = ImGui.GetIO();

            var mouseState = window.MouseState;
            var keyboardState = window.KeyboardState;

            io.MouseDown[0] = mouseState[MouseButton.Left];
            io.MouseDown[1] = mouseState[MouseButton.Right];
            io.MouseDown[2] = mouseState[MouseButton.Middle];

            var screenPoint = new Vector2i((int)mouseState.X, (int)mouseState.Y);
            var point = screenPoint;
            io.MousePos = new System.Numerics.Vector2(point.X, point.Y);

            foreach (Keys key in Enum.GetValues(typeof(Keys)))
            {
                if (key == Keys.Unknown)
                {
                    continue;
                }

                io.KeysDown[(int)key] = keyboardState.IsKeyDown(key);
            }

            foreach (var c in PressedChars)
            {
                io.AddInputCharacter(c);
            }

            PressedChars.Clear();

            io.KeyCtrl = keyboardState.IsKeyDown(Keys.LeftControl) || keyboardState.IsKeyDown(Keys.RightControl);
            io.KeyAlt = keyboardState.IsKeyDown(Keys.LeftAlt) || keyboardState.IsKeyDown(Keys.RightAlt);
            io.KeyShift = keyboardState.IsKeyDown(Keys.LeftShift) || keyboardState.IsKeyDown(Keys.RightShift);
            io.KeySuper = keyboardState.IsKeyDown(Keys.LeftSuper) || keyboardState.IsKeyDown(Keys.RightSuper);
        }

        internal void PressChar(char key)
        {
            this.PressedChars.Add(key);
        }

        internal static void MouseScroll(Vector2 offset)
        {
            var io = ImGui.GetIO();

            io.MouseWheel = offset.Y;
            io.MouseWheelH = offset.X;
        }

        private static void SetKeyMappings()
        {
            var io = ImGui.GetIO();

            io.KeyMap[(int)ImGuiKey.Tab] = (int)Keys.Tab;
            io.KeyMap[(int)ImGuiKey.LeftArrow] = (int)Keys.Left;
            io.KeyMap[(int)ImGuiKey.RightArrow] = (int)Keys.Right;
            io.KeyMap[(int)ImGuiKey.UpArrow] = (int)Keys.Up;
            io.KeyMap[(int)ImGuiKey.DownArrow] = (int)Keys.Down;
            io.KeyMap[(int)ImGuiKey.PageUp] = (int)Keys.PageUp;
            io.KeyMap[(int)ImGuiKey.PageDown] = (int)Keys.PageDown;
            io.KeyMap[(int)ImGuiKey.Home] = (int)Keys.Home;
            io.KeyMap[(int)ImGuiKey.End] = (int)Keys.End;
            io.KeyMap[(int)ImGuiKey.Delete] = (int)Keys.Delete;
            io.KeyMap[(int)ImGuiKey.Backspace] = (int)Keys.Backspace;
            io.KeyMap[(int)ImGuiKey.Enter] = (int)Keys.Enter;
            io.KeyMap[(int)ImGuiKey.Escape] = (int)Keys.Escape;
            io.KeyMap[(int)ImGuiKey.A] = (int)Keys.A;
            io.KeyMap[(int)ImGuiKey.C] = (int)Keys.C;
            io.KeyMap[(int)ImGuiKey.V] = (int)Keys.V;
            io.KeyMap[(int)ImGuiKey.X] = (int)Keys.X;
            io.KeyMap[(int)ImGuiKey.Y] = (int)Keys.Y;
            io.KeyMap[(int)ImGuiKey.Z] = (int)Keys.Z;
        }

        private void Dispose()
        {
            GL.DeleteVertexArray(this.vertexArray);
            GL.DeleteBuffer(this.vertexBuffer);
            GL.DeleteBuffer(this.indexBuffer);

            GL.DeleteTexture(this.fontTexture);
            GL.DeleteProgram(this.shader);
        }

        private void RenderImDrawData(ImDrawDataPtr draw_data)
        {
            if (draw_data.CmdListsCount == 0)
            {
                return;
            }

            // Get intial state.
            int prevVAO = 0; GL.GetInteger(GetPName.VertexArrayBinding, ref prevVAO);
            int prevArrayBuffer = 0; GL.GetInteger(GetPName.ArrayBufferBinding, ref prevArrayBuffer);
            int prevProgram = 0; GL.GetInteger(GetPName.CurrentProgram, ref prevProgram);
            byte prevBlendEnabled = 0; GL.GetBoolean(GetPName.Blend, ref prevBlendEnabled);
            byte prevScissorTestEnabled = 0; GL.GetBoolean(GetPName.ScissorTest, ref prevScissorTestEnabled);
            int prevBlendEquationRgb = 0; GL.GetInteger(GetPName.BlendEquationRgb, ref prevBlendEquationRgb);
            int prevBlendEquationAlpha = 0; GL.GetInteger(GetPName.BlendEquationAlpha, ref prevBlendEquationAlpha);
            int prevBlendFuncSrcRgb = 0; GL.GetInteger(GetPName.BlendSrcRgb, ref prevBlendFuncSrcRgb);
            int prevBlendFuncSrcAlpha = 0; GL.GetInteger(GetPName.BlendSrcAlpha, ref prevBlendFuncSrcAlpha);
            int prevBlendFuncDstRgb = 0; GL.GetInteger(GetPName.BlendDstRgb, ref prevBlendFuncDstRgb);
            int prevBlendFuncDstAlpha = 0; GL.GetInteger(GetPName.BlendDstAlpha, ref prevBlendFuncDstAlpha);
            byte prevCullFaceEnabled = 0; GL.GetBoolean(GetPName.CullFace, ref prevCullFaceEnabled);
            byte prevDepthTestEnabled = 0; GL.GetBoolean(GetPName.DepthTest, ref prevDepthTestEnabled);
            int prevActiveTexture = 0; GL.GetInteger(GetPName.ActiveTexture, ref prevActiveTexture);
            GL.ActiveTexture(TextureUnit.Texture0);
            int prevTexture2D = 0; GL.GetInteger(GetPName.TextureBinding2d, ref prevTexture2D);
            Span<int> prevScissorBox = stackalloc int[4];
            GL.GetInteger(GetPName.ScissorBox, prevScissorBox);
            //unsafe
            //{
            //    fixed (int* iptr = &prevScissorBox[0])
            //    {
            //        GL.GetInteger(GetPName.ScissorBox, iptr);
            //    }
            //}

            // Bind the element buffer (thru the VAO) so that we can resize it.
            GL.BindVertexArray(this.vertexArray);
            // Bind the vertex buffer so that we can resize it.
            GL.BindBuffer(BufferTargetARB.ArrayBuffer, this.vertexBuffer);
            for (int i = 0; i < draw_data.CmdListsCount; i++)
            {
                ImDrawListPtr cmd_list = draw_data.CmdListsRange[i];

                int vertexSize = cmd_list.VtxBuffer.Size * Unsafe.SizeOf<ImDrawVert>();
                if (vertexSize > this.vertexBufferSize)
                {
                    int newSize = (int)Math.Max(this.vertexBufferSize * 1.5f, vertexSize);

                    GL.BufferData(BufferTargetARB.ArrayBuffer, newSize, IntPtr.Zero, BufferUsageARB.DynamicDraw);
                    this.vertexBufferSize = newSize;

                    Console.WriteLine($"Resized dear imgui vertex buffer to new size {this.vertexBufferSize}");
                }

                int indexSize = cmd_list.IdxBuffer.Size * sizeof(ushort);
                if (indexSize > this.indexBufferSize)
                {
                    int newSize = (int)Math.Max(this.indexBufferSize * 1.5f, indexSize);
                    GL.BufferData(BufferTargetARB.ElementArrayBuffer, newSize, IntPtr.Zero, BufferUsageARB.DynamicDraw);
                    this.indexBufferSize = newSize;

                    Console.WriteLine($"Resized dear imgui index buffer to new size {this.indexBufferSize}");
                }
            }

            // Setup orthographic projection matrix into our constant buffer
            ImGuiIOPtr io = ImGui.GetIO();
            Matrix4 mvp = Matrix4.CreateOrthographicOffCenter(
                0.0f,
                io.DisplaySize.X,
                io.DisplaySize.Y,
                0.0f,
                -1.0f,
                1.0f);

            GL.UseProgram(this.shader);
            GL.UniformMatrix4f(this.shaderProjectionMatrixLoc, false, mvp);
            GL.Uniform1i(this.shaderFontTextureLoc, 0);
            CheckGLError("Projection");

            GL.BindVertexArray(this.vertexArray);
            CheckGLError("VAO");

            draw_data.ScaleClipRects(io.DisplayFramebufferScale);

            GL.Enable(EnableCap.Blend);
            //GL.Enable(EnableCap.ScissorTest);
            GL.BlendEquation(BlendEquationModeEXT.FuncAdd);
            GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.OneMinusSrcAlpha);
            GL.Disable(EnableCap.CullFace);
            GL.Disable(EnableCap.DepthTest);

            // Render command lists
            for (int n = 0; n < draw_data.CmdListsCount; n++)
            {
                ImDrawListPtr cmd_list = draw_data.CmdListsRange[n];

                GL.BufferSubData(BufferTargetARB.ArrayBuffer, IntPtr.Zero, cmd_list.VtxBuffer.Size * Unsafe.SizeOf<ImDrawVert>(), cmd_list.VtxBuffer.Data);
                CheckGLError($"Data Vert {n}");

                GL.BufferSubData(BufferTargetARB.ElementArrayBuffer, IntPtr.Zero, cmd_list.IdxBuffer.Size * sizeof(ushort), cmd_list.IdxBuffer.Data);
                CheckGLError($"Data Idx {n}");

                for (int cmd_i = 0; cmd_i < cmd_list.CmdBuffer.Size; cmd_i++)
                {
                    ImDrawCmdPtr pcmd = cmd_list.CmdBuffer[cmd_i];
                    if (pcmd.UserCallback != IntPtr.Zero)
                    {
                        throw new NotImplementedException();
                    }
                    else
                    {
                        GL.ActiveTexture(TextureUnit.Texture0);
                        GL.BindTexture(TextureTarget.Texture2d, new TextureHandle((int)pcmd.TextureId));
                        CheckGLError("Texture");

                        // We do _windowHeight - (int)clip.W instead of (int)clip.Y because gl has flipped Y when it comes to these coordinates
                        //var clip = pcmd.ClipRect;
                        //GL.Scissor((int)clip.X, this.windowHeight - (int)clip.W, (int)(clip.Z - clip.X), (int)(clip.W - clip.Y));
                        //CheckGLError("Scissor");
                        //GL.Scissor(0, 0, windowWidth, windowHeight);

                        if ((io.BackendFlags & ImGuiBackendFlags.RendererHasVtxOffset) != 0)
                        {
                            GL.DrawElementsBaseVertex(PrimitiveType.Triangles, (int)pcmd.ElemCount, DrawElementsType.UnsignedShort, (IntPtr)(pcmd.IdxOffset * sizeof(ushort)), unchecked((int)pcmd.VtxOffset));
                        }
                        else
                        {
                            GL.DrawElements(PrimitiveType.Triangles, (int)pcmd.ElemCount, DrawElementsType.UnsignedShort, (int)pcmd.IdxOffset * sizeof(ushort));
                        }
                        CheckGLError("Draw");
                    }
                }
            }

            GL.Disable(EnableCap.Blend);
            GL.Disable(EnableCap.ScissorTest);

            // Reset state
            GL.BindTexture(TextureTarget.Texture2d, new TextureHandle(prevTexture2D));
            GL.ActiveTexture((TextureUnit)prevActiveTexture);
            GL.UseProgram(new ProgramHandle(prevProgram));
            GL.BindVertexArray(new VertexArrayHandle(prevVAO));
            GL.Scissor(prevScissorBox[0], prevScissorBox[1], prevScissorBox[2], prevScissorBox[3]);
            GL.BindBuffer(BufferTargetARB.ArrayBuffer, new BufferHandle(prevArrayBuffer));
            GL.BlendEquationSeparate((BlendEquationModeEXT)prevBlendEquationRgb, (BlendEquationModeEXT)prevBlendEquationAlpha);
            GL.BlendFuncSeparate(
                (BlendingFactor)prevBlendFuncSrcRgb,
                (BlendingFactor)prevBlendFuncDstRgb,
                (BlendingFactor)prevBlendFuncSrcAlpha,
                (BlendingFactor)prevBlendFuncDstAlpha);
            if (prevBlendEnabled != 0) GL.Enable(EnableCap.Blend); else GL.Disable(EnableCap.Blend);
            if (prevDepthTestEnabled != 0) GL.Enable(EnableCap.DepthTest); else GL.Disable(EnableCap.DepthTest);
            if (prevCullFaceEnabled != 0) GL.Enable(EnableCap.CullFace); else GL.Disable(EnableCap.CullFace);
            if (prevScissorTestEnabled != 0) GL.Enable(EnableCap.ScissorTest); else GL.Disable(EnableCap.ScissorTest);
        }

        private static void CheckGLError(string title)
        {
            OpenTK.Graphics.OpenGL.ErrorCode error;
            int i = 1;
            while ((error = GL.GetError()) != OpenTK.Graphics.OpenGL.ErrorCode.NoError)
            {
                Debug.Print($"{title} ({i++}): {error}");
            }
        }

        internal void WindowResized(int width, int height)
        {
            this.windowWidth = width;
            this.windowHeight = height;
        }

        private static ProgramHandle CreateProgram(string name, string vertexSource, string fragmentSource)
        {
            var program = GL.CreateProgram();

            var vertex = CompileShader(name, ShaderType.VertexShader, vertexSource);
            var fragment = CompileShader(name, ShaderType.FragmentShader, fragmentSource);

            GL.AttachShader(program, vertex);
            GL.AttachShader(program, fragment);

            GL.LinkProgram(program);

            int success = 0;
            GL.GetProgrami(program, ProgramPropertyARB.LinkStatus, ref success);
            //if (success == 0)
            {
                GL.GetProgramInfoLog(program, out var info);
                Debug.WriteLine($"GL.LinkProgram had info log [{name}]:\n{info}");
            }

            GL.DetachShader(program, vertex);
            GL.DetachShader(program, fragment);

            GL.DeleteShader(vertex);
            GL.DeleteShader(fragment);

            return program;
        }

        private static ShaderHandle CompileShader(string name, ShaderType shaderType, string source)
        {
            var shader = GL.CreateShader(shaderType);

            GL.ShaderSource(shader, source);
            GL.CompileShader(shader);

            int success = 0;
            GL.GetShaderi(shader, ShaderParameterName.CompileStatus, ref success);
            //if (success == 0)
            {
                GL.GetShaderInfoLog(shader, out var info);
                Console.WriteLine($"GL.CompileShader for shader '{name}' [{shaderType}] had info log:\n{info}");
            }

            return shader;
        }

        public void DestroyDeviceObjects()
        {
            this.Dispose();
        }
    }
}
