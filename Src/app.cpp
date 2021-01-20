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
    // Read test case:
    TestCaseInfo TInfo;
    m_reader->read_test_case(1001, TInfo);
    PNSimulator sim(TInfo);

    // Create and Clear output folder
    fs::create_directory("./Data/PNData/");
    delete_dir_content("./Data/PNData/");

    // main loop
    for (int i = 0; i < frame_cnt; ++i) {
        // Output Anim sequence
        m_writer->write_anim_seq(i, TInfo.name_path, sim.GetXRef(), TInfo.boundary_tri);
    }
}
