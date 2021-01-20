//
// Created by jiaruiyan on 1/18/21.
//

#include "app.h"
#include <filesystem>

namespace fs = std::filesystem;

void delete_dir_content(const fs::path& dir_path) {
    for (auto& path: fs::directory_iterator(dir_path)) {
        fs::remove_all(path);
    }
}

void App::run(int test_case_id, int frame_cnt){
    // Read test case.
    Eigen::MatrixXd X;
    Eigen::MatrixXi Tet;
    Eigen::MatrixXi BTri;
    m_reader->read_test_case(1001, X, Tet, BTri);

    // Create and Clear output folder
    fs::create_directory("./Data/PNData/");
    delete_dir_content("./Data/PNData/");

    // Output Anim sequence
    std::string tmp = "TMP";
    m_writer->write_anim_seq(0, tmp, X, BTri);
}
