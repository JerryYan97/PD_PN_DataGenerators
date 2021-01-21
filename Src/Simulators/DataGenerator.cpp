//
// Created by jiaruiyan on 1/19/21.
//

#include "DataGenerator.h"
#include <oneapi/tbb.h>

class Parallel_restT_M{
    std::vector<Eigen::Matrix<double, 3, 3>>* my_restT;
    double density;
    int n_verts;
    Eigen::MatrixXd* my_X;
    Eigen::MatrixXi* my_Tet;
public:
    Eigen::VectorXd acc_mass;
    void operator()(const tbb::blocked_range<size_t>& r){
        std::vector<Eigen::Matrix<double, 3, 3>>& restT = *my_restT;
        Eigen::MatrixXd& X = *my_X;
        Eigen::MatrixXi& Tet = *my_Tet;
        for (int i = r.begin(); i != r.end(); ++i) {
            restT[i].col(0) = X.row(Tet(i, 1)) - X.row(Tet(i, 0));
            restT[i].col(1) = X.row(Tet(i, 2)) - X.row(Tet(i, 0));
            restT[i].col(2) = X.row(Tet(i, 3)) - X.row(Tet(i, 0));
            double m = restT[i].determinant() * density / (3.0 * 2.0 * 4.0);
            assert(m > 0);
            for (int j = 0; j < 4; ++j) {
                acc_mass[Tet(i, j)] += m;
            }
        }
    }

    Parallel_restT_M(Parallel_restT_M& x, tbb::split)
    : my_restT(x.my_restT), density(x.density), n_verts(x.n_verts),
    my_X(x.my_X), my_Tet(x.my_Tet), acc_mass(Eigen::VectorXd::Zero(x.n_verts))
    {}

    void join(const Parallel_restT_M& y){
        acc_mass += y.acc_mass;
    }

    Parallel_restT_M(std::vector<Eigen::Matrix<double, 3, 3>>* my_restT, double density, int n_verts,
                     Eigen::MatrixXd* my_X, Eigen::MatrixXi* my_Tet)
                     :  my_restT(my_restT), density(density), n_verts(n_verts),
                     my_X(my_X), my_Tet(my_Tet), acc_mass(Eigen::VectorXd::Zero(n_verts))
    {}
};

void DataGenerator::compute_restT_m(){
    Parallel_restT_M prtm(&restT, density, n_verts, &X, &Tet);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, n_elements), prtm);
    mass = prtm.acc_mass;
}
