#version 330

smooth in vec3 vertexColour;
smooth in vec3 vertexNormal;

uniform int solidColour;

out vec4 colour;

void main()
{
/*
	if (useUniformColour > 0)
	{
		gl_FragColor = colour;
	}
	else
	{	
		vec3 lightDir = -normalize(vec3(0, -5, 30));
		float d = dot(vertexNormal, -lightDir);
		d = max(0.6, d);
		gl_FragColor = vec4(vertexColour, 1) * d;
	}
*/

	if (solidColour > 0)
	{
		colour = vec4(1, 1, 1, 1);
	}
	else
	{	
		vec3 lightDir = -normalize(vec3(0, -5, 30));
		float d = dot(vertexNormal, -lightDir);
		d = max(0.6, d);
		colour = vec4(vertexColour, 1) * d;
	}
}
