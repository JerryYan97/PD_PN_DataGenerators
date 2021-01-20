//
// Created by jiaruiyan on 1/19/21.
//

#include "AppReader.h"
#include <oneapi/tbb.h>

void AppReader::read_test_case(int id,
                    Eigen::MatrixXd &X,
                    Eigen::MatrixXi &Tet,
                    Eigen::VectorXi &TetTag,
                    Eigen::VectorXi &BTriTag){
    // Unused API variables
    Eigen::MatrixXi Tri;
    Eigen::VectorXi TriTag;

    // Call API
    std::string meshStr = "box3D_v518_t2112.msh";
    const std::string readfile = this->m_default_mesh_path + meshStr;
    if (!igl::readMSH(readfile, X, Tri, Tet, TriTag, TetTag)){
        throw std::invalid_argument("Cannot read .msh file!");
    }

    // Find boundary triangles

}

