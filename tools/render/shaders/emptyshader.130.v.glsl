#version 130
uniform mat4 MVP;
attribute vec3 vertex;
varying float f_data;
void main(void)
{
//    gl_Position = MVP * vec4(vertex,1);
    gl_Position = vec4(0.0,0.0,0.0,1.0);
    f_data = 1.0;
}
