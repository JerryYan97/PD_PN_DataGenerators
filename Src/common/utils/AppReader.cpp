//
// Created by jiaruiyan on 1/19/21.
//

#include "AppReader.h"
void AppReader::read_test_case(int id,
                    Eigen::MatrixXd &X,
                    Eigen::MatrixXi &Tet,
                    Eigen::VectorXi &TetTag){
    // Unused API variables
    Eigen::MatrixXi Tri;
    Eigen::VectorXi TriTag;

    // Call API
    const std::string readfile = "./Data/MeshModels/box3D_v518_t2112.msh";
    if (!igl::readMSH(readfile, X, Tri, Tet, TriTag, TetTag)){
        throw std::invalid_argument("Cannot read .msh file!");
    }
}

