#version 130

varying float f_data;
void main(void)
{
    gl_FragColor = vec4(f_data,0,-f_data,1);
}
