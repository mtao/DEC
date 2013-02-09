#version 330

in float f_data;
layout(location=0, index=0) out vec4 fragment_color;
void main(void)
{
    fragment_color = vec4(f_data,0,-f_data,1);
}
