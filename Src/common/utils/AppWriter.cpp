//
// Created by jiaruiyan on 1/19/21.
//

#include "AppWriter.h"
#include <igl/writeOBJ.h>

void AppWriter::write_anim_seq(int frame_id, std::string& filename,
                    Eigen::MatrixXd& X, Eigen::MatrixXi& BTri){
    // Fill the boundary structure -- it should also have a wiser way to output used X and BTri.
    // Output to .obj file
    std::string obj_path = "./Data/PNData/tmp.obj";
    igl::writeOBJ(obj_path, X, BTri);
}
