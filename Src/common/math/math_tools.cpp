//
// Created by jiaruiyan on 1/22/21.
//
#include "math_tools.h"

#include "SVD/ImplicitQRSVD.h"

// IQRSVD
void Singular_Value_Decomposition(const Eigen::Matrix3d& A,
                                  Eigen::Matrix3d& U,
                                  Eigen::Vector3d& singular_values,
                                  Eigen::Matrix3d& V){
    JIXIE::singularValueDecomposition(A, U, singular_values, V);
}

// 3d
void makePD(Eigen::Matrix3d& symMtr)
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(symMtr);
    if (eigenSolver.eigenvalues()[0] >= 0.0) {
        return;
    }
    Eigen::DiagonalMatrix<double, 3> D(eigenSolver.eigenvalues());
    int rows = 3;
    for (int i = 0; i < rows; i++) {
        if (D.diagonal()[i] < 0.0) {
            D.diagonal()[i] = 0.0;
        }
        else {
            break;
        }
    }
    symMtr = eigenSolver.eigenvectors() * D * eigenSolver.eigenvectors().transpose();
}

// 2d
void makePD(Eigen::Matrix2d& symMtr){
    const double a = symMtr(0, 0);
    const double b = (symMtr(0, 1) + symMtr(1, 0)) / 2.0;
    const double d = symMtr(1, 1);

    double b2 = b * b;
    const double D = a * d - b2;
    const double T_div_2 = (a + d) / 2.0;
    const double sqrtTT4D = std::sqrt(T_div_2 * T_div_2 - D);
    const double L2 = T_div_2 - sqrtTT4D;
    if (L2 < 0.0) {
        const double L1 = T_div_2 + sqrtTT4D;
        if (L1 <= 0.0) {
            symMtr.setZero();
        }
        else {
            if (b2 == 0.0) {
                symMtr << L1, 0.0, 0.0, 0.0;
            }
            else {
                const double L1md = L1 - d;
                const double L1md_div_L1 = L1md / L1;
                symMtr(0, 0) = L1md_div_L1 * L1md;
                symMtr(0, 1) = symMtr(1, 0) = b * L1md_div_L1;
                symMtr(1, 1) = b2 / L1;
            }
        }
    }
}

void computeCofactorMtr(const Eigen::Matrix3d& F, Eigen::Matrix3d& A)
{
    A(0, 0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
    A(0, 1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
    A(0, 2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);
    A(1, 0) = F(0, 2) * F(2, 1) - F(0, 1) * F(2, 2);
    A(1, 1) = F(0, 0) * F(2, 2) - F(0, 2) * F(2, 0);
    A(1, 2) = F(0, 1) * F(2, 0) - F(0, 0) * F(2, 1);
    A(2, 0) = F(0, 1) * F(1, 2) - F(0, 2) * F(1, 1);
    A(2, 1) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
    A(2, 2) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
}
