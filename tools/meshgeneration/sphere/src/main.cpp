#include "sphere.hpp"
#include <iostream>

int main(int argc, char * argv[]) {
    if(argc != 3) {
        std::cout << "Usage: sphereGen 4 output.obj" << std::endl;
        return 1;
    }
    int depth = atoi(argv[1]);
    char * filename = argv[2];

    SphereMeshFactory<float> smf(depth);
    smf.write(filename);
    return 0;
}

