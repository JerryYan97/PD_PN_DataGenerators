//
// Created by jiaruiyan on 1/21/21.
//

#include "FixedCoRotMaterial.h"
#include <iostream>
double FixedCoRotMaterial::compute_energy_density(const Eigen::Vector3d& singularValues){
    double sigmam12Sum = (singularValues - Eigen::Vector3d::Ones()).squaredNorm();
    double J = singularValues.prod();
    return mu * sigmam12Sum + 0.5 * la * (J - 1.0) * (J - 1.0);
}

