#version 330
layout(triangles) in;
layout(triangle_strip, max_vertices=3) out;
out float f_data;
void main()
{
    f_data = 1.0;
    gl_Position = gl_in[0].gl_Position;
    EmitVertex();
    gl_Position = gl_in[1].gl_Position;
    EmitVertex();
    gl_Position = gl_in[2].gl_Position;
    EmitVertex();
    EndPrimitive();
}
