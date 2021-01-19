//
// Created by jiaruiyan on 1/17/21.
//
#include "macro.h"
#include <iostream>
#include <time.h>
#include <igl/opengl/glfw/Viewer.h>

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

    //
    // Test MKL
    //
    // Without MKL: 18.7719s
    // With MKL-SEQ: 7.55845s
    // With MKL-TBB: 7.30708s
    srand((unsigned)time(NULL));
    clock_t start,finish;
    double totaltime;
    start=clock();

    Eigen::MatrixXf m1 = Eigen::MatrixXf::Random(7000, 7000);
    Eigen::MatrixXf m2 = Eigen::MatrixXf::Random(7000, 7000);
    Eigen::MatrixXf m3 = m1*m2;

    finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    std::cout << "The running time of codes above:" << totaltime << "s" << std::endl;

    // Inline mesh of a cube
    const Eigen::MatrixXd V= (Eigen::MatrixXd(8,3)<<
                                                  0.0,0.0,0.0,
            0.0,0.0,1.0,
            0.0,1.0,0.0,
            0.0,1.0,1.0,
            1.0,0.0,0.0,
            1.0,0.0,1.0,
            1.0,1.0,0.0,
            1.0,1.0,1.0).finished();
    const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
                                                    1,7,5,
            1,3,7,
            1,4,3,
            1,2,4,
            3,8,7,
            3,4,8,
            5,7,8,
            5,8,6,
            1,5,6,
            1,6,2,
            2,6,8,
            2,8,4).finished().array()-1;

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.launch();

    return 0;
}