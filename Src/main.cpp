//
// Created by jiaruiyan on 1/17/21.
//
#include "macro.h"
#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;

int main(){
    std::cout << "Hello World!" << TEST << std::endl;
    // Test Eigen
    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;
    std::cout << "Eigen version:" << EIGEN_WORLD_VERSION << "."
              << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION << std::endl;
    return 0;
}