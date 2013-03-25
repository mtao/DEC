#include "sphere.hpp"
#include <iostream>

int main(int argc, char * argv[]) {
    if(argc < 2) {
        std::cout << "Usage: sphereGen 4 output.obj" << std::endl;
        return 1;
    }
    int depth = atoi(argv[1]);
    SphereMeshFactory<float> smf(depth);
    if(argc >= 3){
        char * filename = argv[2];
        smf.write(filename);
    } else {
        smf.write(std::cout);
    }

    return 0;
}

