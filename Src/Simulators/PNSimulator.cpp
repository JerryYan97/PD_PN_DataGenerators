//
// Created by jiaruiyan on 1/19/21.
//

#include "PNSimulator.h"

#include "../common/math/math_tools.h"
#include <Eigen/SparseCholesky>

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

double PNSimulator::compute_energy(std::vector<Eigen::Matrix3d>& F, std::vector<Eigen::Matrix3d>& U,
                                   std::vector<Eigen::Vector3d>& Sigma, std::vector<Eigen::Matrix3d>& V,
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
            Singular_Value_Decomposition(F[i], U[i], Sigma[i], V[i]);
        }
        tetE(i) = vol[i] * dt * dt * m_material.compute_energy_density(Sigma[i]);
    });
    energy += tetE.sum();
    // Compute inertia energy
    // Inertia
    Parallel_Compute_Inertia_Energy pcie(&mass, &X, &xTilde);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, n_verts), pcie);
    energy += pcie.acc_energy;
    return energy;
}

// NOTE: I change the rhs to positive to make everything looks like Newton's method.
class Parallel_Hessian_Gradient_Computation_Elasticity{
    const std::vector<Eigen::Matrix3d>* my_F;
    const std::vector<Eigen::Matrix3d>* my_U;
    const std::vector<Eigen::Vector3d>* my_Sigma;
    const std::vector<Eigen::Matrix3d>* my_V;
    const PNSimulator* my_sim;
public:
    Eigen::SparseMatrix<double> acc_H;
    Eigen::VectorXd acc_grad_E;
    void operator()( const tbb::blocked_range<size_t>& r ){
        const std::vector<Eigen::Matrix3d>& F = *my_F;
        const std::vector<Eigen::Matrix3d>& U = *my_U;
        const std::vector<Eigen::Vector3d>& Sigma = *my_Sigma;
        const std::vector<Eigen::Matrix3d>& V = *my_V;
        std::vector<Eigen::Triplet<double>> t_vec;
        t_vec.reserve(12 * 12 * (r.end() - r.begin()));
        
        for(size_t i = r.begin(); i != r.end(); ++i){
            Eigen::Matrix<double, 9, 9> dPdF = Eigen::Matrix<double, 9, 9>::Zero(9, 9);
            my_sim->m_material.first_piola_kirchoff_stress_derivative(U[i], Sigma[i], V[i], dPdF);
            dPdF *= (my_sim->dt * my_sim->dt * my_sim->vol(i));
            // std::cout << "dPdF:" << std::endl << dPdF << std::endl;

            Eigen::Matrix3d P = Eigen::Matrix3d::Zero();
            my_sim->m_material.first_piola_kirchoff_stress(F[i], U[i], Sigma[i], V[i], P);
            P *= (my_sim->dt * my_sim->dt * my_sim->vol(i));
            // std::cout << "P:" << std::endl << P << std::endl;

            Eigen::Matrix<double, 12, 9> intermediate = Eigen::Matrix<double, 12, 9>::Zero();
            const Eigen::Matrix3d& IB = my_sim->restTInv[i];
            for (int colI = 0; colI < 9; ++colI) {
                intermediate(3, colI) = IB(0, 0) * dPdF(0, colI) + IB(0, 1) * dPdF(3, colI) + IB(0, 2) * dPdF(6, colI);
                intermediate(4, colI) = IB(0, 0) * dPdF(1, colI) + IB(0, 1) * dPdF(4, colI) + IB(0, 2) * dPdF(7, colI);
                intermediate(5, colI) = IB(0, 0) * dPdF(2, colI) + IB(0, 1) * dPdF(5, colI) + IB(0, 2) * dPdF(8, colI);
                intermediate(6, colI) = IB(1, 0) * dPdF(0, colI) + IB(1, 1) * dPdF(3, colI) + IB(1, 2) * dPdF(6, colI);
                intermediate(7, colI) = IB(1, 0) * dPdF(1, colI) + IB(1, 1) * dPdF(4, colI) + IB(1, 2) * dPdF(7, colI);
                intermediate(8, colI) = IB(1, 0) * dPdF(2, colI) + IB(1, 1) * dPdF(5, colI) + IB(1, 2) * dPdF(8, colI);
                intermediate(9, colI) = IB(2, 0) * dPdF(0, colI) + IB(2, 1) * dPdF(3, colI) + IB(2, 2) * dPdF(6, colI);
                intermediate(10, colI) = IB(2, 0) * dPdF(1, colI) + IB(2, 1) * dPdF(4, colI) + IB(2, 2) * dPdF(7, colI);
                intermediate(11, colI) = IB(2, 0) * dPdF(2, colI) + IB(2, 1) * dPdF(5, colI) + IB(2, 2) * dPdF(8, colI);
                intermediate(0, colI) = -intermediate(3, colI) - intermediate(6, colI) - intermediate(9, colI);
                intermediate(1, colI) = -intermediate(4, colI) - intermediate(7, colI) - intermediate(10, colI);
                intermediate(2, colI) = -intermediate(5, colI) - intermediate(8, colI) - intermediate(11, colI);
            }

            std::vector<int> indMap(12);
            indMap[0] = my_sim->Tet(i, 0) * 3;
            indMap[1] = my_sim->Tet(i, 0) * 3 + 1;
            indMap[2] = my_sim->Tet(i, 0) * 3 + 2;
            indMap[3] = my_sim->Tet(i, 1) * 3;
            indMap[4] = my_sim->Tet(i, 1) * 3 + 1;
            indMap[5] = my_sim->Tet(i, 1) * 3 + 2;
            indMap[6] = my_sim->Tet(i, 2) * 3;
            indMap[7] = my_sim->Tet(i, 2) * 3 + 1;
            indMap[8] = my_sim->Tet(i, 2) * 3 + 2;
            indMap[9] = my_sim->Tet(i, 3) * 3;
            indMap[10] = my_sim->Tet(i, 3) * 3 + 1;
            indMap[11] = my_sim->Tet(i, 3) * 3 + 2;

            for(int rowI = 0; rowI < 12; ++rowI){
                std::vector<Eigen::Triplet<double>> tmp_vec(12);
                tmp_vec[0] = Eigen::Triplet(indMap[rowI], indMap[3], IB(0, 0) * intermediate(rowI, 0) + IB(0, 1) * intermediate(rowI, 3) + IB(0, 2) * intermediate(rowI, 6));
                tmp_vec[1] = Eigen::Triplet(indMap[rowI], indMap[4], IB(0, 0) * intermediate(rowI, 1) + IB(0, 1) * intermediate(rowI, 4) + IB(0, 2) * intermediate(rowI, 7));
                tmp_vec[2] = Eigen::Triplet(indMap[rowI], indMap[5], IB(0, 0) * intermediate(rowI, 2) + IB(0, 1) * intermediate(rowI, 5) + IB(0, 2) * intermediate(rowI, 8));
                tmp_vec[3] = Eigen::Triplet(indMap[rowI], indMap[6], IB(1, 0) * intermediate(rowI, 0) + IB(1, 1) * intermediate(rowI, 3) + IB(1, 2) * intermediate(rowI, 6));
                tmp_vec[4] = Eigen::Triplet(indMap[rowI], indMap[7], IB(1, 0) * intermediate(rowI, 1) + IB(1, 1) * intermediate(rowI, 4) + IB(1, 2) * intermediate(rowI, 7));
                tmp_vec[5] = Eigen::Triplet(indMap[rowI], indMap[8], IB(1, 0) * intermediate(rowI, 2) + IB(1, 1) * intermediate(rowI, 5) + IB(1, 2) * intermediate(rowI, 8));
                tmp_vec[6] = Eigen::Triplet(indMap[rowI], indMap[9], IB(2, 0) * intermediate(rowI, 0) + IB(2, 1) * intermediate(rowI, 3) + IB(2, 2) * intermediate(rowI, 6));
                tmp_vec[7] = Eigen::Triplet(indMap[rowI], indMap[10], IB(2, 0) * intermediate(rowI, 1) + IB(2, 1) * intermediate(rowI, 4) + IB(2, 2) * intermediate(rowI, 7));
                tmp_vec[8] = Eigen::Triplet(indMap[rowI], indMap[11], IB(2, 0) * intermediate(rowI, 2) + IB(2, 1) * intermediate(rowI, 5) + IB(2, 2) * intermediate(rowI, 8));
                tmp_vec[9] = Eigen::Triplet(indMap[rowI], indMap[0], -tmp_vec[0].value() - tmp_vec[3].value() - tmp_vec[6].value());
                tmp_vec[10] = Eigen::Triplet(indMap[rowI], indMap[1], -tmp_vec[1].value() - tmp_vec[4].value() - tmp_vec[7].value());
                tmp_vec[11] = Eigen::Triplet(indMap[rowI], indMap[2], -tmp_vec[2].value() - tmp_vec[5].value() - tmp_vec[8].value());
                for (int j = 0; j < 12; ++j) {
                    int row_idx = tmp_vec[j].row();
                    int col_idx = tmp_vec[j].col();
                    if(!my_sim->m_fixed_label[row_idx] && !my_sim->m_fixed_label[col_idx]){
                        t_vec.push_back(tmp_vec[j]);
                    }
                }
            }
            double R10 = IB(0, 0) * P(0, 0) + IB(0, 1) * P(0, 1) + IB(0, 2) * P(0, 2);
            double R11 = IB(0, 0) * P(1, 0) + IB(0, 1) * P(1, 1) + IB(0, 2) * P(1, 2);
            double R12 = IB(0, 0) * P(2, 0) + IB(0, 1) * P(2, 1) + IB(0, 2) * P(2, 2);
            double R20 = IB(1, 0) * P(0, 0) + IB(1, 1) * P(0, 1) + IB(1, 2) * P(0, 2);
            double R21 = IB(1, 0) * P(1, 0) + IB(1, 1) * P(1, 1) + IB(1, 2) * P(1, 2);
            double R22 = IB(1, 0) * P(2, 0) + IB(1, 1) * P(2, 1) + IB(1, 2) * P(2, 2);
            double R30 = IB(2, 0) * P(0, 0) + IB(2, 1) * P(0, 1) + IB(2, 2) * P(0, 2);
            double R31 = IB(2, 0) * P(1, 0) + IB(2, 1) * P(1, 1) + IB(2, 2) * P(1, 2);
            double R32 = IB(2, 0) * P(2, 0) + IB(2, 1) * P(2, 1) + IB(2, 2) * P(2, 2);
            if(!my_sim->m_fixed_label[my_sim->Tet(i, 1) * 3]){
                acc_grad_E[my_sim->Tet(i, 1) * 3 + 0] += R10;
                acc_grad_E[my_sim->Tet(i, 1) * 3 + 1] += R11;
                acc_grad_E[my_sim->Tet(i, 1) * 3 + 2] += R12;
            }
            if(!my_sim->m_fixed_label[my_sim->Tet(i, 2) * 3]){
                acc_grad_E[my_sim->Tet(i, 2) * 3 + 0] += R20;
                acc_grad_E[my_sim->Tet(i, 2) * 3 + 1] += R21;
                acc_grad_E[my_sim->Tet(i, 2) * 3 + 2] += R22;
            }
            if(!my_sim->m_fixed_label[my_sim->Tet(i, 3) * 3]){
                acc_grad_E[my_sim->Tet(i, 3) * 3 + 0] += R30;
                acc_grad_E[my_sim->Tet(i, 3) * 3 + 1] += R31;
                acc_grad_E[my_sim->Tet(i, 3) * 3 + 2] += R32;
            }
            if(!my_sim->m_fixed_label[my_sim->Tet(i, 0) * 3]){
                acc_grad_E[my_sim->Tet(i, 0) * 3 + 0] += -R10 - R20 - R30;
                acc_grad_E[my_sim->Tet(i, 0) * 3 + 1] += -R11 - R21 - R31;
                acc_grad_E[my_sim->Tet(i, 0) * 3 + 2] += -R12 - R22 - R32;
            }
        }
        Eigen::SparseMatrix<double> local_H(my_sim->n_verts * 3, my_sim->n_verts * 3);
        local_H.setFromTriplets(t_vec.begin(), t_vec.end());
        acc_H += local_H;
    }

    Parallel_Hessian_Gradient_Computation_Elasticity(Parallel_Hessian_Gradient_Computation_Elasticity& x, tbb::split)
    : my_F(x.my_F), my_U(x.my_U), my_Sigma(x.my_Sigma), my_V(x.my_V), my_sim(x.my_sim),
    acc_H(Eigen::SparseMatrix<double>(x.my_sim->n_verts * 3, x.my_sim->n_verts * 3)),
    acc_grad_E(Eigen::VectorXd::Zero(x.my_sim->n_verts * 3))
    {}

    void join(const Parallel_Hessian_Gradient_Computation_Elasticity& y){
        acc_H += y.acc_H;
        acc_grad_E += y.acc_grad_E;
    }

    Parallel_Hessian_Gradient_Computation_Elasticity(const std::vector<Eigen::Matrix3d>* F,
                                                     const std::vector<Eigen::Matrix3d>* U,
                                                     const std::vector<Eigen::Vector3d>* Sigma,
                                                     const std::vector<Eigen::Matrix3d>* V,
                                                     const PNSimulator* sim)
    : my_F(F), my_U(U), my_Sigma(Sigma), my_V(V), my_sim(sim),
    acc_H(Eigen::SparseMatrix<double>(sim->n_verts * 3, sim->n_verts * 3)),
    acc_grad_E(Eigen::VectorXd::Zero(sim->n_verts * 3))
    {}
};

void PNSimulator::compute_hessian_and_gradient(Eigen::SparseMatrix<double>& H, Eigen::VectorXd& grad_E,
                                               const Eigen::MatrixXd& xTilde,
                                               const std::vector<Eigen::Matrix<double, 3, 3>>& F,
                                               const std::vector<Eigen::Matrix3d>& U,
                                               const std::vector<Eigen::Vector3d>& Sigma,
                                               const std::vector<Eigen::Matrix3d>& V) {
    // Inerita and dirichlet mat
    // std::vector<Eigen::Triplet<double>> triplet_vec;
    // triplet_vec.reserve(3 * n_verts);
    tbb::concurrent_vector<Eigen::Triplet<double>> triplet_vec;
    triplet_vec.reserve(3 * n_verts);
    tbb::parallel_for(size_t(0), size_t(n_verts), [&](size_t i){
        if(!m_fixed_label[3 * i]){
            triplet_vec.push_back(Eigen::Triplet<double>(3 * i, 3 * i, mass[i]));
            triplet_vec.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 1, mass[i]));
            triplet_vec.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 2, mass[i]));
            grad_E[3 * i] += (mass[i] * (X(i, 0) - xTilde(i, 0)));
            grad_E[3 * i + 1] += (mass[i] * (X(i, 1) - xTilde(i, 1)));
            grad_E[3 * i + 2] += (mass[i] * (X(i, 2) - xTilde(i, 2)));
        }else{
            triplet_vec.push_back(Eigen::Triplet<double>(3 * i, 3 * i, 1.0));
            triplet_vec.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 1, 1.0));
            triplet_vec.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 2, 1.0));
        }
    });
    H.setFromTriplets(triplet_vec.begin(), triplet_vec.end());
    // std::cout << "Hessian:" << std::endl << H << std::endl;
    // std::cout << "Grad_E:" << std::endl << grad_E << std::endl;
    // Elasticity
    Parallel_Hessian_Gradient_Computation_Elasticity phgce(&F, &U, &Sigma, &V, this);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, n_elements), phgce);
    // std::cout << "phgce.acc_H:" << std::endl << phgce.acc_H << std::endl;
    H += phgce.acc_H;
    grad_E += phgce.acc_grad_E;
}

void PNSimulator::step() {
    std::cout.precision(17);
    // Compute xn and x tilde
    Eigen::MatrixXd Xprev(X);
    Eigen::MatrixXd Xt(X);
    Eigen::MatrixXd xTilde(n_verts, 3);
    tbb::parallel_for(size_t(0), size_t(n_verts), [&](size_t i){
        Eigen::Vector3d f = m_force_field->GetForce(X.row(i));
        xTilde(i, 0) = X(i, 0) + dt * vel(i, 0) + (dt * dt * f[0] / mass(i));
        xTilde(i, 1) = X(i, 1) + dt * vel(i, 1) + (dt * dt * f[1] / mass(i));
        xTilde(i, 2) = X(i, 2) + dt * vel(i, 2) + (dt * dt * f[2] / mass(i));
    });
    // std::cout << "m:" << mass << std::endl;
    // Compute deformation gradient and Eprev
    std::vector<Eigen::Matrix3d> F(n_elements);
    std::vector<Eigen::Matrix3d> U_Vec(n_elements);
    std::vector<Eigen::Matrix3d> V_Vec(n_elements);
    std::vector<Eigen::Vector3d> Sigma_Vec(n_elements);
    double residual = 0.0;
    double one_over_dt = 1.0 / dt;
    double Eprev = compute_energy(F, U_Vec, Sigma_Vec, V_Vec, xTilde, true);
    std::cout << "Eprev:" << Eprev << std::endl;
    do {
        // Construct Hessian and gradient
        Eigen::SparseMatrix<double> H(n_verts * 3, n_verts * 3);
        Eigen::VectorXd grad_E = Eigen::VectorXd::Zero(n_verts * 3);
        compute_hessian_and_gradient(H, grad_E, xTilde, F, U_Vec, Sigma_Vec, V_Vec);
        // std::cout << "Hessian:" << std::endl << H << std::endl;
        // std::cout << "Grad_E:" << std::endl << grad_E << std::endl;
        // Calculate sol
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> LLT_Solver;
        LLT_Solver.compute(H);
        Eigen::VectorXd p = -LLT_Solver.solve(grad_E);
        /*
        std::cout << "p == 1 idx:" << std::endl;
        for (int i = 0; i < 3 * n_verts; ++i) {
            if(p[i] == -1.0){
                std::cout << i << std::endl;
            }
        }
         */
        // std::cout << p << std::endl;
        double alpha = 1.0;
        double Ex = 0.0;
        do {
            // Apply sol
            std::cout << "Applied alpha:" << alpha << std::endl;
            tbb::parallel_for(size_t(0), size_t(n_verts), [&](size_t i){
                X(i, 0) = Xprev(i, 0) + alpha * p[3 * i];
                X(i, 1) = Xprev(i, 1) + alpha * p[3 * i + 1];
                X(i, 2) = Xprev(i, 2) + alpha * p[3 * i + 2];
            });
            // alpha decrease
            alpha *= 0.5;
            Ex = compute_energy(F, U_Vec, Sigma_Vec, V_Vec, xTilde, true);
        }while (Ex > Eprev);
        Eprev = Ex;
        Xprev = X;
        residual = p.cwiseAbs().maxCoeff();
        /*
        std::cout << "p abs min:" << p.cwiseAbs().minCoeff() << std::endl;
        std::cout << "p abs max:" << p.cwiseAbs().maxCoeff() << std::endl;
        std::cout << "p == 1 idx:" << p.cwiseAbs().maxCoeff() << std::endl;
        for (int i = 0; i < 3 * n_verts; ++i) {
            if(p[i] == 1.0){
                std::cout << i << std::endl;
            }
        }
        std::cout << "grad_E == 0 idx:" << p.cwiseAbs().maxCoeff() << std::endl;
        for (int i = 0; i < 3 * n_verts; ++i) {
            if(grad_E[i] == 0.0){
                std::cout << i << std::endl;
            }
        }
        */
        std::cout << "grad E abs min:" << grad_E.cwiseAbs().minCoeff() << std::endl;
        std::cout << "Search Direction Residual : " << residual * one_over_dt << std::endl;
        std::cout << "Ex : " << Ex << std::endl;
    }while ((residual * one_over_dt) > 0.0001);
    // Calculate velocity
    vel = (X - Xt) * one_over_dt;
    // std::cout << vel.row(0) << std::endl;
}
