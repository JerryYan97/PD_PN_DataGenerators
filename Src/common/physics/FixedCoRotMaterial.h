//
// Created by jiaruiyan on 1/21/21.
//

#ifndef PD_PN_GENERATORS_FIXEDCOROTMATERIAL_H
#define PD_PN_GENERATORS_FIXEDCOROTMATERIAL_H

#define EIGEN_USE_MKL_ALL
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
    void compute_dE_div_dsigma(const Eigen::Vector3d& singularValues, Eigen::Vector3d& dE_div_dsigma) const;
    void compute_d2E_div_dsigma2(const Eigen::Vector3d& singularValues, Eigen::Matrix3d& d2E_div_dsigma2) const;
    void compute_BLeftCoef(const Eigen::Vector3d& singularValues, Eigen::Vector3d& BLeftCoef) const;


    void first_piola_kirchoff_stress_derivative(const Eigen::Matrix3d& U, const Eigen::Vector3d& Sigma,
                                                const Eigen::Matrix3d& V, Eigen::Matrix<double, 9, 9>& dPdF) const;

    void first_piola_kirchoff_stress(const Eigen::Matrix3d& F, const Eigen::Matrix3d& U, const Eigen::Vector3d& Sigma,
                                     const Eigen::Matrix3d& V, Eigen::Matrix3d& P) const;

    double GetDensity(){
        return density;
    }
};


#endif //PD_PN_GENERATORS_FIXEDCOROTMATERIAL_H
