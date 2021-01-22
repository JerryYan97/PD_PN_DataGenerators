//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_DATAGENERATOR_H
#define PD_PN_GENERATORS_DATAGENERATOR_H

#ifndef EIGEN_VECTORIZE_SSE4_2
#define EIGEN_VECTORIZE_SSE4_2
#endif

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include <memory>
#include <iostream>
#include "../common/ForceFields/DirectForceField.h"
#include "../common/utils/InfoStructs.h"
#include "../common/physics/FixedCoRotMaterial.h"
#include <oneapi/tbb.h>

class DataGenerator {
protected:
    double dt;
    FixedCoRotMaterial m_material;

    int n_elements;
    int n_verts;

    Eigen::MatrixXd X;
    Eigen::MatrixXd vel;
    Eigen::MatrixXi Tet;

    Eigen::VectorXd mass;
    Eigen::VectorXd vol;
    std::vector<Eigen::Matrix<double, 3, 3>> restT;
    std::vector<Eigen::Matrix<double, 3, 3>> restTInv;
    std::shared_ptr<ForceField> m_force_field;

    void compute_restT_m();

public:
    DataGenerator(TestCaseInfo& info, std::shared_ptr<ForceField> ff)
    : X(info.init_X), Tet(info.tet), dt(info.dt)
    {
        double la = info.E * info.nu / ((1.0 + info.nu) * (1.0 - 2.0 * info.nu));
        double mu = info.E / (2.0 * (1.0 + info.nu));
        m_material = FixedCoRotMaterial(la, mu, info.density);
        m_force_field = std::move(ff);
        n_elements = info.tet.rows();
        n_verts = info.init_X.rows();

        restT.resize(n_elements);
        restTInv.resize(n_elements);
        vol = Eigen::VectorXd::Zero(n_elements);
        mass = Eigen::VectorXd::Zero(n_verts);
        vel = Eigen::MatrixXd::Zero(X.rows(), X.cols());
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
