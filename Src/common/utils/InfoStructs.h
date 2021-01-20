//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_INFOSTRUCTS_H
#define PD_PN_GENERATORS_INFOSTRUCTS_H
#include <string>
#include <Eigen/Dense>
struct TestCaseInfo{
    std::string name;
    std::string name_path;
    int id;
    double dt;
    double E;
    double mu;
    Eigen::MatrixXd init_X;
    Eigen::MatrixXi boundary_tri;
    Eigen::MatrixXi tet;
};

#endif //PD_PN_GENERATORS_INFOSTRUCTS_H
