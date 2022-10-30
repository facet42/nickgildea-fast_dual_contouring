namespace FastDC
{
    internal struct ViewerOptions
    {
        public float MeshScale;
        public bool DrawWireframe;
        public bool UseCms;
        public bool RefreshModel;

        public ViewerOptions()
        {
            this.MeshScale = 1.0f;
            this.DrawWireframe = false;
            this.UseCms = true;
            this.RefreshModel = false;
        }
    }
}
