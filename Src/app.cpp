//
// Created by jiaruiyan on 1/18/21.
//

#include "app.h"

void App::run(int test_case_id, int frame_cnt){
    // Read test case.
    Eigen::MatrixXd X;
    Eigen::MatrixXi Tet;
    Eigen::VectorXi TetTag;
    Eigen::VectorXi BTriTag;
    m_reader->read_test_case(1001, X, Tet, TetTag, BTriTag);

    // Output Anim sequence
    // m_writer->write_anim_seq();
}
