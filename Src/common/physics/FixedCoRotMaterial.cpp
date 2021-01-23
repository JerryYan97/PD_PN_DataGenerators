//
// Created by jiaruiyan on 1/21/21.
//

#include "FixedCoRotMaterial.h"
#include "../math/math_tools.h"
#include <iostream>

double FixedCoRotMaterial::compute_energy_density(const Eigen::Vector3d& singularValues){
    double sigmam12Sum = (singularValues - Eigen::Vector3d::Ones()).squaredNorm();
    double J = singularValues.prod();
    return mu * sigmam12Sum + 0.5 * la * (J - 1.0) * (J - 1.0);
}

void FixedCoRotMaterial::compute_d2E_div_dsigma2(const Eigen::Vector3d& singularValues,
                                                 Eigen::Matrix3d& d2E_div_dsigma2)
                                                 const {
    const double sigmaProd = singularValues.prod();
    Eigen::Vector3d sigmaProd_noI;
    sigmaProd_noI[0] = singularValues[1] * singularValues[2];
    sigmaProd_noI[1] = singularValues[2] * singularValues[0];
    sigmaProd_noI[2] = singularValues[0] * singularValues[1];
    double _2u = mu * 2;
    d2E_div_dsigma2(0, 0) = _2u + la * sigmaProd_noI[0] * sigmaProd_noI[0];
    d2E_div_dsigma2(1, 1) = _2u + la * sigmaProd_noI[1] * sigmaProd_noI[1];
    d2E_div_dsigma2(2, 2) = _2u + la * sigmaProd_noI[2] * sigmaProd_noI[2];
    d2E_div_dsigma2(0, 1) = d2E_div_dsigma2(1, 0) = la * (singularValues[2] * (sigmaProd - 1.0) + sigmaProd_noI[0] * sigmaProd_noI[1]);
    d2E_div_dsigma2(0, 2) = d2E_div_dsigma2(2, 0) = la * (singularValues[1] * (sigmaProd - 1.0) + sigmaProd_noI[0] * sigmaProd_noI[2]);
    d2E_div_dsigma2(2, 1) = d2E_div_dsigma2(1, 2) = la * (singularValues[0] * (sigmaProd - 1.0) + sigmaProd_noI[2] * sigmaProd_noI[1]);
}

void FixedCoRotMaterial::compute_dE_div_dsigma(const Eigen::Vector3d& singularValues,
                                               Eigen::Vector3d& dE_div_dsigma) const{
    const double sigmaProdm1lambda = la * (singularValues.prod() - 1.0);
    Eigen::Vector3d sigmaProd_noI;
    sigmaProd_noI[0] = singularValues[1] * singularValues[2];
    sigmaProd_noI[1] = singularValues[2] * singularValues[0];
    sigmaProd_noI[2] = singularValues[0] * singularValues[1];
    double _2u = mu * 2;
    dE_div_dsigma[0] = (_2u * (singularValues[0] - 1.0) + sigmaProd_noI[0] * sigmaProdm1lambda);
    dE_div_dsigma[1] = (_2u * (singularValues[1] - 1.0) + sigmaProd_noI[1] * sigmaProdm1lambda);
    dE_div_dsigma[2] = (_2u * (singularValues[2] - 1.0) + sigmaProd_noI[2] * sigmaProdm1lambda);
}

void FixedCoRotMaterial::compute_BLeftCoef(const Eigen::Vector3d& singularValues, Eigen::Vector3d& BLeftCoef) const{
    const double sigmaProd = singularValues.prod();
    const double halfLambda = la / 2.0;
    BLeftCoef[0] = mu - halfLambda * singularValues[2] * (sigmaProd - 1);
    BLeftCoef[1] = mu - halfLambda * singularValues[0] * (sigmaProd - 1);
    BLeftCoef[2] = mu - halfLambda * singularValues[1] * (sigmaProd - 1);
}

void FixedCoRotMaterial::first_piola_kirchoff_stress(
        const Eigen::Matrix3d& F, const Eigen::Matrix3d& U, const Eigen::Vector3d& Sigma,
        const Eigen::Matrix3d& V, Eigen::Matrix3d& P) const {
    Eigen::Matrix3d JFInvT;
    computeCofactorMtr(F, JFInvT);
    P = 2.0 * mu * (F - U * V.transpose()) + la * (Sigma.prod() - 1.0) * JFInvT;
}

void FixedCoRotMaterial::first_piola_kirchoff_stress_derivative(const Eigen::Matrix3d& U,
                                                                const Eigen::Vector3d& Sigma,
                                                                const Eigen::Matrix3d& V,
                                                                Eigen::Matrix<double, 9, 9>& dPdF) const{
    // compute A
    Eigen::Vector3d dE_div_dsigma;
    compute_dE_div_dsigma(Sigma, dE_div_dsigma);
    Eigen::Matrix3d d2E_div_dsigma2;
    compute_d2E_div_dsigma2(Sigma, d2E_div_dsigma2);
    makePD(d2E_div_dsigma2);

    // compute B
    Eigen::Vector3d BLeftCoef;
    compute_BLeftCoef(Sigma, BLeftCoef);
    Eigen::Matrix2d B[3];
    for (int cI = 0; cI < 3; cI++){
        int cI_post = (cI + 1) % 3;
        double rightCoef = dE_div_dsigma[cI] + dE_div_dsigma[cI_post];
        double sum_sigma = Sigma[cI] + Sigma[cI_post];
        const double eps = 1.0e-6;
        if (sum_sigma < eps) {
            rightCoef /= 2.0 * eps;
        }
        else {
            rightCoef /= 2.0 * sum_sigma;
        }
        const double& leftCoef = BLeftCoef[cI];
        B[cI](0, 0) = B[cI](1, 1) = leftCoef + rightCoef;
        B[cI](0, 1) = B[cI](1, 0) = leftCoef - rightCoef;
        makePD(B[cI]);
    }

    // compute M using A(d2E_div_dsigma2) and B
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(9, 9);
    // A
    M(0, 0) = d2E_div_dsigma2(0, 0);
    M(0, 4) = d2E_div_dsigma2(0, 1);
    M(0, 8) = d2E_div_dsigma2(0, 2);
    M(4, 0) = d2E_div_dsigma2(1, 0);
    M(4, 4) = d2E_div_dsigma2(1, 1);
    M(4, 8) = d2E_div_dsigma2(1, 2);
    M(8, 0) = d2E_div_dsigma2(2, 0);
    M(8, 4) = d2E_div_dsigma2(2, 1);
    M(8, 8) = d2E_div_dsigma2(2, 2);
    // B01
    M(1, 1) = B[0](0, 0);
    M(1, 3) = B[0](0, 1);
    M(3, 1) = B[0](1, 0);
    M(3, 3) = B[0](1, 1);
    // B12
    M(5, 5) = B[1](0, 0);
    M(5, 7) = B[1](0, 1);
    M(7, 5) = B[1](1, 0);
    M(7, 7) = B[1](1, 1);
    // B20
    M(2, 2) = B[2](1, 1);
    M(2, 6) = B[2](1, 0);
    M(6, 2) = B[2](0, 1);
    M(6, 6) = B[2](0, 0);

    // compute dP_div_dF
    for (int i = 0; i < 3; i++) {
        int _dim_i = i * 3;
        for (int j = 0; j < 3; j++) {
            int ij = _dim_i + j;
            for (int r = 0; r < 3; r++) {
                int _dim_r = r * 3;
                for (int s = 0; s < 3; s++) {
                    int rs = _dim_r + s;
                    if (ij > rs) {
                        // bottom left, same as upper right
                        continue;
                    }
                    dPdF(ij, rs) = M(0, 0) * U(i, 0) * V(j, 0) * U(r, 0) * V(s, 0) + M(0, 4) * U(i, 0) * V(j, 0) * U(r, 1) * V(s, 1) + M(0, 8) * U(i, 0) * V(j, 0) * U(r, 2) * V(s, 2) + M(4, 0) * U(i, 1) * V(j, 1) * U(r, 0) * V(s, 0) + M(4, 4) * U(i, 1) * V(j, 1) * U(r, 1) * V(s, 1) + M(4, 8) * U(i, 1) * V(j, 1) * U(r, 2) * V(s, 2) + M(8, 0) * U(i, 2) * V(j, 2) * U(r, 0) * V(s, 0) + M(8, 4) * U(i, 2) * V(j, 2) * U(r, 1) * V(s, 1) + M(8, 8) * U(i, 2) * V(j, 2) * U(r, 2) * V(s, 2) + M(1, 1) * U(i, 0) * V(j, 1) * U(r, 0) * V(s, 1) + M(1, 3) * U(i, 0) * V(j, 1) * U(r, 1) * V(s, 0) + M(3, 1) * U(i, 1) * V(j, 0) * U(r, 0) * V(s, 1) + M(3, 3) * U(i, 1) * V(j, 0) * U(r, 1) * V(s, 0) + M(5, 5) * U(i, 1) * V(j, 2) * U(r, 1) * V(s, 2) + M(5, 7) * U(i, 1) * V(j, 2) * U(r, 2) * V(s, 1) + M(7, 5) * U(i, 2) * V(j, 1) * U(r, 1) * V(s, 2) + M(7, 7) * U(i, 2) * V(j, 1) * U(r, 2) * V(s, 1) + M(2, 2) * U(i, 0) * V(j, 2) * U(r, 0) * V(s, 2) + M(2, 6) * U(i, 0) * V(j, 2) * U(r, 2) * V(s, 0) + M(6, 2) * U(i, 2) * V(j, 0) * U(r, 0) * V(s, 2) + M(6, 6) * U(i, 2) * V(j, 0) * U(r, 2) * V(s, 0);
                    if (ij < rs) {
                        dPdF(rs, ij) = dPdF(ij, rs);
                    }
                }
            }
        }
    }
}