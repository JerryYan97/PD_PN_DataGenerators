//
// Created by jiaruiyan on 1/18/21.
//

#include "app.h"
#include <filesystem>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>

namespace fs = std::filesystem;

void delete_dir_content(const fs::path& dir_path) {
    for (auto& path: fs::directory_iterator(dir_path)) {
        fs::remove_all(path);
    }
}

void App::set_force_field(TestCaseInfo Tinfo, std::shared_ptr<ForceField> ff) {
    // TODO: Implement a more general force field settings.
}

void App::runVisualization(){
    igl::opengl::glfw::Viewer viewer;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOFF("./Data/MeshModels/bunny.off", V, F);
    viewer.launch_init();
    while (true){
        viewer.data().set_mesh(V, F);
        viewer.launch_rendering(false);
    }
}

void App::run(int test_case_id, int frame_cnt, bool visualization){
    // Read test case:
    TestCaseInfo TInfo;
    m_reader->read_test_case(1001, TInfo);
    // Set force field:
    ForceFieldInfo ffinfo;
    ffinfo.dir_force = Eigen::Vector3d(0.1, 0.1, 0.1);
    ffinfo.type = DIRECT_FORCE;
    TInfo.force_field->SetForceField(ffinfo);
    // Create and init simulator:
    PNSimulator sim(TInfo, TInfo.force_field);

    // Visualization
    igl::opengl::glfw::Viewer viewer;
    if (visualization){
        viewer.launch_init();
        viewer.data().set_mesh(sim.GetXRef(), TInfo.boundary_tri);
        viewer.launch_rendering(false);
    }

    // Create and Clear output folder
    fs::create_directory("./Data/PNData/");
    delete_dir_content("./Data/PNData/");

    // main loop
    for (int i = 0; i < frame_cnt; ++i) {
        std::cout << "==================== Frame: " << i << " ====================" << std::endl;
        sim.step();
        // Output Anim sequence
        m_writer->write_anim_seq(i, TInfo.name_path, sim.GetXRef(), TInfo.boundary_tri);
        if (visualization){
            viewer.data().set_mesh(sim.GetXRef(), TInfo.boundary_tri);
            viewer.launch_rendering(false);
        }
    }
    if (visualization){
        viewer.launch_shut();
    }
}

void App::runTestCase() {
    // Read test case:
    TestCaseInfo TInfo;
    m_reader->read_test_case(1001, TInfo);
    // Set force field:
    ForceFieldInfo ffinfo;
    ffinfo.dir_force = Eigen::Vector3d(0.1, 0.1, 0.1);
    ffinfo.type = DIRECT_FORCE;
    TInfo.force_field->SetForceField(ffinfo);
    // Create and init simulator:
    PNSimulator sim(TInfo, TInfo.force_field);

    sim.hessian_test_case();
}
