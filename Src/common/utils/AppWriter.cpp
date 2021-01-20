//
// Created by jiaruiyan on 1/19/21.
//

#include "AppWriter.h"

void AppWriter::write_anim_seq(int frame_id, std::string& filename,
                    Eigen::MatrixXd& X, Eigen::VectorXi& TetTag, Eigen::VectorXi& BoundaryTriTag){
    // Fill the boundary structure
    Eigen::VectorXi o_boundary_tri_idx;
    Eigen::MatrixXd o_boundary_tri_X;
    for(int i = 0; i < BoundaryTriTag.rows(); ++i){

    }
    // Output to .obj file
}
