//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_DATAGENERATOR_H
#define PD_PN_GENERATORS_DATAGENERATOR_H

#include <memory>
#include "ForceFields/DirectForceField.h"
#include "../common/utils/InfoStructs.h"


class DataGenerator {
protected:
    double dt;
    double la;
    double mu;
    double rho;

    Eigen::MatrixXd X;
    Eigen::MatrixXi Tet;
public:
    DataGenerator(TestCaseInfo& info)
    : X(info.init_X), Tet(info.tet), dt(info.dt)
    {}
};


#endif //PD_PN_GENERATORS_DATAGENERATOR_H
