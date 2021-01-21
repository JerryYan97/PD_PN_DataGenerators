//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_INFOSTRUCTS_H
#define PD_PN_GENERATORS_INFOSTRUCTS_H
#include <string>
#include <memory>
#include <Eigen/Dense>
#include "../ForceFields/DirectForceField.h"

// Here, the concept of test case is a more general idea than the testcase in python.
struct TestCaseInfo{
    std::string name;
    std::string name_path;
    int id;

    // Small sim variables
    double dt;

    // Material properties
    double E;
    double nu;
    double density;

    // Large sim variables
    Eigen::MatrixXd init_X;
    Eigen::MatrixXi boundary_tri;
    Eigen::MatrixXi tet;
    Eigen::MatrixXi dirichlet;

    // Force field configuration info
    std::shared_ptr<ForceField> force_field;
};

#endif //PD_PN_GENERATORS_INFOSTRUCTS_H
