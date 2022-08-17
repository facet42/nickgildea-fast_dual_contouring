namespace FastDC
{
    using System;

    using ImGuiNET;

    using OpenTK.Graphics.OpenGL;
    using OpenTK.Mathematics;
    using OpenTK.Windowing.Common;
    using OpenTK.Windowing.Desktop;

    internal class AppWindow : GameWindow
    {
        ImGuiController controller;

        public AppWindow(int width, int height)
            : base(
                  GameWindowSettings.Default,
                  new NativeWindowSettings()
                  {
                      Size = new Vector2i(width, height),
                      APIVersion = new Version(3, 3)
                  })
        {
        }

        protected override void OnLoad()
        {
            base.OnLoad();

            this.Title += ": OpenGL Version " + GL.GetString(StringName.Version);

            this.controller = new ImGuiController(this.ClientSize.X, this.ClientSize.Y);
        }

        protected override void OnResize(ResizeEventArgs e)
        {
            base.OnResize(e);

            GL.Viewport(0, 0, this.ClientSize.X, this.ClientSize.Y);

            this.controller.WindowResized(this.ClientSize.X, this.ClientSize.Y);
        }

        protected override void OnRenderFrame(FrameEventArgs e)
        {
            base.OnRenderFrame(e);

            this.controller.Update(this, e.Time);

            GL.ClearColor(0, 32, 48, 255);
            GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit | ClearBufferMask.StencilBufferBit);

            ImGui.ShowDemoWindow();

            this.controller.Render();

            this.SwapBuffers();
        }

        protected override void OnTextInput(TextInputEventArgs e)
        {
            base.OnTextInput(e);

            this.controller.PressChar((char)e.Unicode);
        }

        protected override void OnMouseWheel(MouseWheelEventArgs e)
        {
            base.OnMouseWheel(e);

            ImGuiController.MouseScroll(e.Offset);
        }
    }
}
