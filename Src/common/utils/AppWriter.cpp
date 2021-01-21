//
// Created by jiaruiyan on 1/19/21.
//

#include "AppWriter.h"

void AppWriter::write_anim_seq(int frame_id, std::string& filename_path,
                    Eigen::MatrixXd& X, Eigen::MatrixXi& BTri){
    // Fill the boundary structure -- it should also have a wiser way to output used X and BTri.
    // Output to .obj file
    int frame_id_len = 8;
    std::string fid_str = std::to_string(frame_id);
    std::string obj_path = filename_path + std::string(frame_id_len - fid_str.length(), '0').append(fid_str.append(".obj"));
    igl::writeOBJ(obj_path, X, BTri);
}
