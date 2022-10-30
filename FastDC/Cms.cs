namespace Facet42
{
    using FastDC;

    using GlmSharp;

    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;
    using System.Threading.Tasks;

    internal static class Cms
    {
        private static readonly uint[] EncodedEdgeOffsets =
        {
            0x00000000,
            0x00100000,
            0x00000400,
            0x00100400,
            0x40000000,
            0x40100000,
            0x40000001,
            0x40100001,
            0x80000000,
            0x80000400,
            0x80000001,
            0x80000401,
        };

        public static void GenerateVertexData(VoxelInfo info, ref MeshBuffer buffer, Log? log)
        {
            log?.WriteLine($"GenerateVertexData{Environment.NewLine}==================");

            var idxCounter = 0;
            foreach (var voxelId in info.Voxels)
            {
                log?.WriteLine($"VID:{voxelId}");

                var p = new vec4[12];
                var n = new vec4[12];

                var idx = 0;
                for (int i = 0; i < 12; i++)
                {
                    var edgeId = voxelId + EncodedEdgeOffsets[i];
                    log?.WriteLine($"EID:{edgeId}");

                    if (info.Edges.ContainsKey(edgeId) == true)
                    {
                        p[idx] = info.Edges[edgeId].Position;
                        n[idx] = info.Edges[edgeId].Normal;

                        log?.WriteLine($"P:{p[idx].x:G9},{p[idx].y:G9},{p[idx].z:G9},{p[idx].w:G9}");
                        log?.WriteLine($"N:{n[idx].x:G9},{n[idx].y:G9},{n[idx].z:G9},{n[idx].w:G9}");

                        idx++;
                    }
                }
            }
        }

        public static void GenerateTriangles(VoxelInfo voxelInfo, ref MeshBuffer buffer)
        {
        }
    }
}
