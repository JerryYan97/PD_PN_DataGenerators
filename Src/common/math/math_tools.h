//
// Created by jiaruiyan on 1/21/21.
//
//

#ifndef PD_PN_GENERATORS_MATH_TOOLS_H
#define PD_PN_GENERATORS_MATH_TOOLS_H

#ifndef EIGEN_VECTORIZE_SSE4_2
#define EIGEN_VECTORIZE_SSE4_2
#endif

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include <Eigen/Dense>
#include "SVD/ImplicitQRSVD.h"

// IQRSVD
void Singular_Value_Decomposition(const Eigen::Matrix3d& A,
                                  Eigen::Matrix3d& U,
                                  Eigen::Vector3d& singular_values,
                                  Eigen::Matrix3d& V){
    JIXIE::singularValueDecomposition(A, U, singular_values, V);
}

#endif //PD_PN_GENERATORS_MATH_TOOLS_H
