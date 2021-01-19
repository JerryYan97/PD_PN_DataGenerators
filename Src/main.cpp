//
// Created by jiaruiyan on 1/17/21.
//
#include "macro.h"
#include <iostream>
#include "app.h"

int main(){
    // NOTE: My eigen version is 3.3.7
    std::cout << "Eigen version:" << EIGEN_WORLD_VERSION << "."
              << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION << std::endl;
    App my_app;
    my_app.run(1001, 80);
    return 0;
}