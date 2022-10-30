namespace Facet42
{
    using FastDC;
    using GlmSharp;

    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;
    using System.Threading.Tasks;

    internal static class DualContour
    {
        private static readonly vec4[] AxisOffset =
        {
            new vec4(1f, 0f, 0f, 0f),
            new vec4(0f, 1f, 0f, 0f),
            new vec4(0f, 0f, 1f, 0f)
        };

        private static readonly ivec4[][] EdgeNodeOffsets =
        {
            new[] { new ivec4(0), new ivec4(0, 0, 1, 0), new ivec4(0, 1, 0, 0), new ivec4(0, 1, 1, 0) },
            new[] { new ivec4(0), new ivec4(1, 0, 0, 0), new ivec4(0, 0, 1, 0), new ivec4(1, 0, 1, 0) },
            new[] { new ivec4(0), new ivec4(0, 1, 0, 0), new ivec4(1, 0, 0, 0), new ivec4(1, 1, 0, 0) },
        };

        private static readonly uint[] EncodedEdgeNodeOffsets =
        {
            0x00000000,
            0x00100000,
            0x00000400,
            0x00100400,
            0x00000000,
            0x00000001,
            0x00100000,
            0x00100001,
            0x00000000,
            0x00000400,
            0x00000001,
            0x00000401,
        };

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

        public static MeshBuffer GenerateMesh(SuperPrimitiveConfig config, int voxelGridSize, Log? log)
        {
            var voxelInfo = FindActiveVoxels(config, voxelGridSize, log);

            var buffer = new MeshBuffer();

            DualContour.GenerateVertexData(voxelInfo, ref buffer);
            DualContour.GenerateTriangles(voxelInfo, ref buffer);

            //this.Dump(config, viewerOpts, voxelInfo, buffer);

            Console.WriteLine($"Mesh: {buffer.NumVertices} {buffer.NumTriangles}");

            return buffer;
        }

        private static VoxelInfo FindActiveVoxels(SuperPrimitiveConfig config, int VoxelGridSize, Log? log)
        {
            var VoxelGridOffset = VoxelGridSize / 2;
            var info = new VoxelInfo();

            for (int x = 0; x < VoxelGridSize; x++)
            {
                for (int y = 0; y < VoxelGridSize; y++)
                {
                    for (int z = 0; z < VoxelGridSize; z++)
                    {
                        var idxPos = new ivec4(x, y, z, 0);
                        var p = new vec4(x - VoxelGridOffset, y - VoxelGridOffset, z - VoxelGridOffset, 1f);

                        for (uint axis = 0; axis < 3; axis++)
                        {
                            var q = p + AxisOffset[axis];

                            var pDensity = config.Density(p);
                            var qDensity = config.Density(q);

                            var zeroCrossing = pDensity >= 0f && qDensity < 0f || pDensity < 0f && qDensity >= 0f;

                            if (zeroCrossing == false)
                            {
                                continue;
                            }

                            var t = FindIntersection(config, p, q);
                            var pos = new vec4(glm.Mix(new vec3(p), new vec3(q), new vec3(t)), 1f);

                            log?.WriteLine($"P:{pDensity:G9} Q:{qDensity:G9} T:{t:G9} POS:{pos.x:G9},{pos.y:G9},{pos.z:G9},{pos.w:G9}");

                            const float H = 0.001f;
                            var d0 = config.Density(pos + new vec4(H, 0f, 0f, 0f));
                            var d1 = config.Density(pos - new vec4(H, 0f, 0f, 0f));

                            var d2 = config.Density(pos + new vec4(0f, H, 0f, 0f));
                            var d3 = config.Density(pos - new vec4(0f, H, 0f, 0f));

                            var d4 = config.Density(pos + new vec4(0f, 0f, H, 0f));
                            var d5 = config.Density(pos - new vec4(0f, 0f, H, 0f));

                            log?.WriteLine($"{d0.ToString("G9").ToLowerInvariant()},{d1.ToString("G9").ToLowerInvariant()},{d2.ToString("G9").ToLowerInvariant()},{d3.ToString("G9").ToLowerInvariant()},{d4.ToString("G9").ToLowerInvariant()},{d5.ToString("G9").ToLowerInvariant()}");

                            var n = new vec4(d0 - d1, d2 - d3, d4 - d5, 0f);
                            var l = (1f / (float)Math.Sqrt(vec4.Dot(n, n)));
                            var normal = n * l;

                            log?.WriteLine($"L:{n.Length:g9}");
                            log?.WriteLine($"N:{n.x:g9},{n.y:g9},{n.z:g9},{n.w:g9}");
                            log?.WriteLine($"N:{normal.x:g9},{normal.y:g9},{normal.z:g9},{normal.w:g9}");

                            var edgeInfo = new EdgeInfo() { Position = pos, Normal = normal, Winding = pDensity >= 0f };

                            var code = EncodeAxisUniqueID(axis, x, y, z);
                            info.Edges.Add(code, edgeInfo);

                            log?.WriteLine($"Edge:{code}");

                            var edgeNodes = EdgeNodeOffsets[axis];

                            for (int i = 0; i < 4; i++)
                            {
                                var nodeIdxPos = idxPos - edgeNodes[i];
                                var nodeID = EncodeVoxelUniqueID(nodeIdxPos);
                                info.Voxels.Add(nodeID);

                                log?.WriteLine($"Node:{nodeID}");
                            }
                        }
                    }
                }
            }

            return info;
        }

        private static uint EncodeVoxelUniqueID(ivec4 idxPos)
        {
            uint i = (uint)idxPos.x | ((uint)idxPos.y << 10) | ((uint)idxPos.z << 20);

            //log?.WriteLine($"{idxPos.x},{idxPos.y},{idxPos.z},{idxPos.w} -> {i}");

            return (uint)i;
        }

        private static uint EncodeAxisUniqueID(uint axis, int x, int y, int z)
        {
            return (uint)x | ((uint)y << 10) | ((uint)z << 20) | (axis << 30);
        }

        private static float FindIntersection(SuperPrimitiveConfig config, vec4 p0, vec4 p1)
        {
            const int FindEdgeInfoSteps = 16;
            const float FindEdgeInfoIncrement = 1f / FindEdgeInfoSteps;

            var minValue = float.MaxValue;
            var currentT = 0f;
            var t = 0f;

            for (int i = 0; i < FindEdgeInfoSteps; i++)
            {
                var p = glm.Mix(p0, p1, new vec4(currentT));
                var d = glm.Abs(config.Density(p));
                if (d < minValue)
                {
                    t = currentT;
                    minValue = d;
                }

                currentT += FindEdgeInfoIncrement;
            }

            return t;
        }

        public static void GenerateVertexData(VoxelInfo info, ref MeshBuffer buffer)
        {
            //log?.WriteLine($"GenerateVertexData{Environment.NewLine}==================");
            //log?.WriteLine($"Voxels:{info.Voxels.Count}");

            var idxCounter = 0;
            foreach (var voxelId in info.Voxels)
            {
                //log?.WriteLine($"VID:{voxelId}");

                var p = new vec4[12];
                var n = new vec4[12];

                var idx = 0;
                for (int i = 0; i < 12; i++)
                {
                    var edgeId = voxelId + EncodedEdgeOffsets[i];
                    //log?.WriteLine($"EID:{edgeId}");

                    if (info.Edges.ContainsKey(edgeId) == true)
                    {
                        p[idx] = info.Edges[edgeId].Position;
                        n[idx] = info.Edges[edgeId].Normal;

                        //log?.WriteLine($"P:{p[idx].x:G9},{p[idx].y:G9},{p[idx].z:G9},{p[idx].w:G9}");
                        //log?.WriteLine($"N:{n[idx].x:G9},{n[idx].y:G9},{n[idx].z:G9},{n[idx].w:G9}");

                        idx++;
                    }
                }

                var error = Qef.SolveFromPoints4d(p, n, idx, out var nodePos);
                var nodeNormal = vec3.Zero;
                for (int i = 0; i < idx; i++)
                {
                    nodeNormal += n[i].xyz;
                }
                nodeNormal *= (1f / (float)idx);

                //log?.WriteLine($"NP:{nodePos.x:G9},{nodePos.y:G9},{nodePos.z:G9}");
                //log?.WriteLine($"NN:{nodeNormal.x:G9},{nodeNormal.y:G9},{nodeNormal.z:G9}");

                info.Indices.Add(voxelId, idxCounter++);

                buffer.Vertices.Add(new MeshVertex() { Position = nodePos, Normal = nodeNormal });
            }
        }

        public static void GenerateTriangles(VoxelInfo voxelInfo, ref MeshBuffer buffer)
        {
            foreach (var pair in voxelInfo.Edges)
            {
                var edge = pair.Key;
                var info = pair.Value;

                //var basePos = DecodeVoxelUniqueID(edge);
                var axis = (edge >> 30) & 0xff;

                var nodeId = edge & ~0xc0000000;
                var voxelIds = new[]
                {
                    nodeId - EncodedEdgeNodeOffsets[axis * 4 + 0],
                    nodeId - EncodedEdgeNodeOffsets[axis * 4 + 1],
                    nodeId - EncodedEdgeNodeOffsets[axis * 4 + 2],
                    nodeId - EncodedEdgeNodeOffsets[axis * 4 + 3],
                };

                var edgeVoxels = new int[4];
                var numFoundVoxels = 0;

                for (int i = 0; i < 4; i++)
                {
                    edgeVoxels[numFoundVoxels++] = voxelInfo.Indices[voxelIds[i]];
                }

                if (numFoundVoxels < 4)
                {
                    continue;
                }

                if (info.Winding == true)
                {
                    buffer.Triangles.Add(new MeshTriangle() { V0 = edgeVoxels[0], V1 = edgeVoxels[1], V2 = edgeVoxels[3] });
                    buffer.Triangles.Add(new MeshTriangle() { V0 = edgeVoxels[0], V1 = edgeVoxels[3], V2 = edgeVoxels[2] });
                }
                else
                {
                    buffer.Triangles.Add(new MeshTriangle() { V0 = edgeVoxels[0], V1 = edgeVoxels[3], V2 = edgeVoxels[1] });
                    buffer.Triangles.Add(new MeshTriangle() { V0 = edgeVoxels[0], V1 = edgeVoxels[2], V2 = edgeVoxels[3] });
                }
            }
        }
    }
}
