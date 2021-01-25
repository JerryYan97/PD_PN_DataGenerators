//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_PNSIMULATOR_H
#define PD_PN_GENERATORS_PNSIMULATOR_H
#define EIGEN_USE_MKL_ALL
#include "DataGenerator.h"
#include <Eigen/SparseCholesky>

class Parallel_Hessian_Gradient_Computation_Elasticity;
class Parallel_Hessian_Gradient_Computation_Inerita;

class PNSimulator : public DataGenerator {
private:
    double compute_energy(std::vector<Eigen::Matrix3d>& F, std::vector<Eigen::Matrix3d>& U,
                          std::vector<Eigen::Vector3d>& Sigma, std::vector<Eigen::Matrix3d>& V,
                          Eigen::MatrixXd& xTilde, bool calF);
    void compute_hessian_and_gradient(Eigen::SparseMatrix<double>& H, Eigen::VectorXd& grad_E,
                                      const Eigen::MatrixXd& xTilde, const std::vector<Eigen::Matrix3d>& F,
                                      const std::vector<Eigen::Matrix3d>& U, const std::vector<Eigen::Vector3d>& Sigma,
                                      const std::vector<Eigen::Matrix3d>& V);



    friend Parallel_Hessian_Gradient_Computation_Elasticity;
    friend Parallel_Hessian_Gradient_Computation_Inerita;
public:
    PNSimulator(TestCaseInfo &info, std::shared_ptr<ForceField> ff) : DataGenerator(info, ff) {}
    virtual void step();
    void hessian_test_case();
};


#endif //PD_PN_GENERATORS_PNSIMULATOR_H
