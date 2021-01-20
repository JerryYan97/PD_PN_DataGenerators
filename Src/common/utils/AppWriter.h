//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_APPWRITER_H
#define PD_PN_GENERATORS_APPWRITER_H

#include <igl/writeOBJ.h>

class AppWriter {
private:
    std::string m_default_output_path;

public:
    AppWriter(){
        m_default_output_path = "./Data/PNData";
    }

    void write_anim_seq(int frame_id, std::string& filename,
                        Eigen::MatrixXd& X, Eigen::VectorXi& TetTag, Eigen::VectorXi& BoundaryTriTag);


};


#endif //PD_PN_GENERATORS_APPWRITER_H
