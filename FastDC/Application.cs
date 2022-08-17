namespace FastDC
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;
    using System.Threading.Tasks;

    using OpenTK.Graphics.OpenGL;
    using OpenTK.Windowing.Desktop;

    internal class Application
    {
        private AppWindow window;

        public int Run(int width, int height)
        {
            if (Initialise(width, height) == false)
            {
                return -1;
            }

            this.window?.Run();

            Shutdown();

            return 0;
        }

        private bool Initialise(int width, int height)
        {
            this.window = new AppWindow(width, height);

            //ReportError(SDL_GL_SetAttribute(SDL_GLattr.SDL_GL_CONTEXT_MAJOR_VERSION, 3));
            //ReportError(SDL_GL_SetAttribute(SDL_GLattr.SDL_GL_CONTEXT_MINOR_VERSION, 3));
            //ReportError(SDL_GL_SetAttribute(SDL_GLattr.SDL_GL_CONTEXT_PROFILE_MASK, SDL_GLprofile.SDL_GL_CONTEXT_PROFILE_CORE));
            //ReportError(SDL_GL_SetAttribute(SDL_GLattr.SDL_GL_DOUBLEBUFFER, 1));

            //this.context = SDL_GL_CreateContext(window);
            //if (this.context == IntPtr.Zero)
            //{
            //    return false;
            //}

            //ReportError(SDL_GL_MakeCurrent(this.window, this.context));

            //GL.Viewport(0, 0, width, height);

            return true;
        }

        private void Shutdown()
        {
            //SDL_DestroyWindow(this.window);
            //SDL_Quit();
        }

        private static void ReportError(int error)
        {
            //if (error >= 0)
            //{
            //    return;
            //}

            //Console.WriteLine(SDL_GetError());
        }
    }
}
