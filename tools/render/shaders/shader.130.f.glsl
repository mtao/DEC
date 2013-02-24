#version 130

varying float f_data;
void main(void)
{
    /*
    if(abs(f_data) < 0.004) {
        discard;
    }
    */
    gl_FragColor = vec4(f_data,0,-f_data,1);
}
