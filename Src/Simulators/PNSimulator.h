//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_PNSIMULATOR_H
#define PD_PN_GENERATORS_PNSIMULATOR_H
#include "DataGenerator.h"

class PNSimulator : public DataGenerator {
private:

public:

    PNSimulator(TestCaseInfo &info) : DataGenerator(info) {}

    void step();

    Eigen::MatrixXi& GetTetRef(){
        return Tet;
    }

    Eigen::MatrixXd& GetXRef(){
        return X;
    }
};


#endif //PD_PN_GENERATORS_PNSIMULATOR_H
