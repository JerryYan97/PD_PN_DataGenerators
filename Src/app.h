//
// Created by jiaruiyan on 1/18/21.
//

#ifndef PD_PN_GENERATORS_APP_H
#define PD_PN_GENERATORS_APP_H
#include <memory>
#include "common/utils/AppReader.h"
#include "common/utils/AppWriter.h"
#include "Simulators/PNSimulator.h"
#include <igl/opengl/glfw/Viewer.h>

class App{

private:
    std::unique_ptr<AppReader> m_reader;
    std::unique_ptr<AppWriter> m_writer;
    std::unique_ptr<DataGenerator> m_data_generator;

public:
    App(){
        m_reader = std::make_unique<AppReader>();
        m_writer = std::make_unique<AppWriter>();
    }
    void run(int test_case_id, int frame_cnt);
    ~App()= default;
};

#endif //PD_PN_GENERATORS_APP_H
