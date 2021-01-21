//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_DATAGENERATOR_H
#define PD_PN_GENERATORS_DATAGENERATOR_H

#include <memory>
#include <iostream>
#include "../common/ForceFields/DirectForceField.h"
#include "../common/utils/InfoStructs.h"
#include <oneapi/tbb.h>

class DataGenerator {
protected:
    double dt;
    double la;
    double mu;
    double density;

    int n_elements;
    int n_verts;

    Eigen::MatrixXd X;
    Eigen::MatrixXi Tet;

    Eigen::VectorXd mass;
    std::vector<Eigen::Matrix<double, 3, 3>> restT;
    std::shared_ptr<ForceField> m_force_field;

    void compute_restT_m();

public:
    DataGenerator(TestCaseInfo& info, std::shared_ptr<ForceField> ff)
    : X(info.init_X), Tet(info.tet), dt(info.dt), density(info.density)
    {
        la = info.E * info.nu / ((1.0 + info.nu) * (1.0 - 2.0 * info.nu));
        mu = info.E / (2.0 * (1.0 + info.nu));
        m_force_field = std::move(ff);
        n_elements = info.tet.rows();
        n_verts = info.init_X.rows();

        restT.resize(n_elements);
        mass = Eigen::VectorXd::Zero(n_verts);
        compute_restT_m();
    }

    virtual void step() = 0;

    Eigen::MatrixXi& GetTetRef(){
        return Tet;
    }

    Eigen::MatrixXd& GetXRef(){
        return X;
    }
};


#endif //PD_PN_GENERATORS_DATAGENERATOR_H
