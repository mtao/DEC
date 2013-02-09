#version 330
uniform mat4 MVP;
in vec3 vertex;
in float data;
out float g_data;
void main(void)
{
    gl_Position = MVP * vec4(vertex,1);
    g_data = data;
}
