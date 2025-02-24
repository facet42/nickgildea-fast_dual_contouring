#version 330

layout(location=0) in vec3 position;
layout(location=1) in vec3 normal;
//layout(location=2) in vec3 colour;

uniform mat4 MVP;
smooth out vec3 vertexColour;
smooth out vec3 vertexNormal;

void main()
{
	vertexColour = vec3(.5f, 0.5f, 0.5f) * (vec3(0.5f) + (vec3(0.5f) * normal));
	//vertexColour = colour;
	vertexNormal = normal;

	gl_Position = MVP * vec4(position, 1);
}
