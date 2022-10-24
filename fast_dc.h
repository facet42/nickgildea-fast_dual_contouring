#ifndef 	HAS_DC_H_BEEN_INCLUDED
#define		HAS_DC_H_BEEN_INCLUDED

#include <unordered_map>
#include <unordered_set>

#include	"ng_mesh_simplify.h"

struct ViewerOptions
{
	float meshScale = 1.f;
	bool drawWireframe = false;
	bool refreshModel = false;
};

struct SuperPrimitiveConfig
{
	enum Type
	{
		Cube,
		Cylinder,
		Pill,
		Corridor,
		Torus,
	};

	glm::vec4 s;
	glm::vec2 r;
};

SuperPrimitiveConfig ConfigForShape(const SuperPrimitiveConfig::Type& type);
MeshBuffer* GenerateMesh(const SuperPrimitiveConfig& config, const ViewerOptions& viewerOpts);

void Dump(SuperPrimitiveConfig& config, std::unordered_set<uint32_t> activeVoxels, std::unordered_map<uint32_t, struct EdgeInfo> activeEdges, MeshBuffer* buffer);

#endif //	HAS_DC_H_BEEN_INCLUDED
