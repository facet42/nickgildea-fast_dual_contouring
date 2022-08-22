namespace FastDC
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;
    using System.Threading.Tasks;
    using ImGuiNET;

    using OpenTK.Graphics;
    using OpenTK.Graphics.OpenGL;
    using OpenTK.Windowing.Desktop;

    internal class Application
    {
        private AppWindow? window;

        FastDC fastDC = new FastDC();

        public int Run(int width, int height)
        {
            if (Initialise(width, height) == false)
            {
                return -1;
            }

            this.window?.Run();

            return 0;
        }

        private bool Initialise(int width, int height)
        {
            this.window = new AppWindow(width, height);
            this.window.Render += Render;
            this.window.RenderUserInterface += RenderUserInterface;
            this.window.Exit += Exit;

            return this.fastDC.Initialise();
        }

        private void Render(object? sender, EventArgs e)
        {
            this.fastDC.Render();
        }

        private void RenderUserInterface(object? sender, EventArgs e)
        {
            this.fastDC.RenderUserInterface();
        }

        private void Exit(object? sender, EventArgs e)
        {
            this.fastDC.Dispose();
        }
    }
}
