//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_APPREADER_H
#define PD_PN_GENERATORS_APPREADER_H

#include <igl/readMSH.h>
#include "InfoStructs.h"

class AppReader {
private:
    std::string m_default_mesh_path;

public:
    AppReader() : m_default_mesh_path("./Data/MeshModels/")
    {}

    void read_test_case(int id, TestCaseInfo& info);
};


#endif //PD_PN_GENERATORS_APPREADER_H
