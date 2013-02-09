#version 330
uniform mat4 MVP;
uniform mat4 MV;
in vec3 vertex;
in float data;
out float g_data;
out vec3 g_pos;
void main(void)
{
    gl_Position = MVP * vec4(vertex,1);
    g_pos = (MV * vec4(vertex,1)).xyz;
    g_data = data;
}
