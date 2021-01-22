//
// Created by jiaruiyan on 1/19/21.
//

#include "PNSimulator.h"

#include "../common/math/math_tools.h"

class Parallel_Compute_Inertia_Energy{
    Eigen::VectorXd* m_mass;
    Eigen::MatrixXd* m_X;
    Eigen::MatrixXd* m_xTilde;
public:
    double acc_energy;
    void operator()( const tbb::blocked_range<size_t>& r ){
        Eigen::VectorXd& mass = *m_mass;
        Eigen::MatrixXd& X = *m_X;
        Eigen::MatrixXd& xTilde = *m_xTilde;
        for(size_t i = r.begin(); i != r.end(); ++i){
            acc_energy += 0.5 * mass(i) * (X.row(i) - xTilde.row(i)).squaredNorm();
        }
    }

    Parallel_Compute_Inertia_Energy(Parallel_Compute_Inertia_Energy& x, tbb::split)
    : m_mass(x.m_mass), m_X(x.m_X), m_xTilde(x.m_xTilde), acc_energy(0.0)
    {}

    void join(const Parallel_Compute_Inertia_Energy& y){
        acc_energy += y.acc_energy;
    }

    Parallel_Compute_Inertia_Energy(Eigen::VectorXd* mass, Eigen::MatrixXd* X, Eigen::MatrixXd* xTilde)
    : m_mass(mass), m_X(X), m_xTilde(xTilde), acc_energy(0.0)
    {}
};

double PNSimulator::compute_energy(std::vector<Eigen::Matrix<double, 3, 3>>& F,
                                   Eigen::MatrixXd& xTilde, bool calF){
    double energy(0.0);
    // Compute deformation gradient and material energy
    Eigen::VectorXd tetE = Eigen::VectorXd::Zero(n_elements);
    tbb::parallel_for(size_t(0), size_t(n_elements), [&](size_t i){
        if(calF){
            Eigen::Matrix3d Xt;
            Eigen::Vector4i Teti = Tet.row(i);
            Xt.col(0) = X.row(Teti(1)) - X.row(Teti(0));
            Xt.col(1) = X.row(Teti(2)) - X.row(Teti(0));
            Xt.col(2) = X.row(Teti(3)) - X.row(Teti(0));
            F[i] = Xt * restTInv[i];
        }
        Eigen::Matrix3d Ui, Vi;
        Eigen::Vector3d Si;
        Singular_Value_Decomposition(F[i], Ui, Si, Vi);
        tetE(i) = vol[i] * dt * dt * m_material.compute_energy_density(Si);
    });
    energy += tetE.sum();
    // Compute inertia energy
    // Inertia
    Parallel_Compute_Inertia_Energy pcie(&mass, &X, &xTilde);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, n_verts), pcie);
    energy += pcie.acc_energy;
    return energy;
}

void PNSimulator::step() {
    // Compute xn and x tilde
    Eigen::MatrixXd Xprev(X);
    Eigen::MatrixXd xTilde = X + dt * vel;
    tbb::parallel_for(size_t(0), size_t(n_verts), [&](size_t i){
        xTilde.row(i) += (dt * dt * m_force_field->GetForce(X.row(i)) / mass(i));
    });

    // Compute deformation gradient and Eprev
    std::vector<Eigen::Matrix<double, 3, 3>> F(n_elements);
    double Eprev = compute_energy(F, xTilde, true);
    /*
    do {
        // Construct Hessian
        // Construct rhs
        // Calculate sol
        double alpha = 1.0;
        do {
            // Apply sol
            // alpha decrease
            alpha *= 0.5;
        }while (compute_energy(F, xTilde, true) > Eprev);
    }while (true);
    */
}
