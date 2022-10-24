namespace FastDC
{
    using System;
    using System.IO;
    using System.Linq;
    using System.Text;
    using System.Threading.Tasks;

    internal class Log : IDisposable
    {
        private StreamWriter file;

        public Log(string path)
        {
            this.file = new StreamWriter(path, false, Encoding.Default);
            this.file.AutoFlush = true;
        }

        public void WriteLine(string message, params object?[] args)
        {
            this.file.WriteLine(message, args);
        }

        public void Dispose()
        {
            this.file.Close();
            this.file.Dispose();
        }
    }
}
