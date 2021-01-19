//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_APPREADER_H
#define PD_PN_GENERATORS_APPREADER_H

#include <igl/readMSH.h>

class AppReader {
public:
    void read_test_case(int id,
                        Eigen::MatrixXd &X,
                        Eigen::MatrixXi &Tet,
                        Eigen::VectorXi &TetTag);
};


#endif //PD_PN_GENERATORS_APPREADER_H
