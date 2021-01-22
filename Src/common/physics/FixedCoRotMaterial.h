//
// Created by jiaruiyan on 1/21/21.
//

#ifndef PD_PN_GENERATORS_FIXEDCOROTMATERIAL_H
#define PD_PN_GENERATORS_FIXEDCOROTMATERIAL_H

#include <vector>
#include <Eigen/Dense>
class FixedCoRotMaterial {
protected:
    double la;
    double mu;
    double density;

public:
    FixedCoRotMaterial() = default;
    FixedCoRotMaterial(double la, double mu, double density)
    : la(la), mu(mu), density(density)
    {}

    double compute_energy_density(const Eigen::Vector3d& singularValues);

    double GetDensity(){
        return density;
    }
};


#endif //PD_PN_GENERATORS_FIXEDCOROTMATERIAL_H
