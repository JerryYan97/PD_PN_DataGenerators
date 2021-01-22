//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_PNSIMULATOR_H
#define PD_PN_GENERATORS_PNSIMULATOR_H
#include "DataGenerator.h"

class PNSimulator : public DataGenerator {
private:
    double compute_energy(std::vector<Eigen::Matrix<double, 3, 3>>& F,
                          Eigen::MatrixXd& xTilde, bool calF);
public:
    PNSimulator(TestCaseInfo &info, std::shared_ptr<ForceField> ff) : DataGenerator(info, ff) {}
    virtual void step();
};


#endif //PD_PN_GENERATORS_PNSIMULATOR_H
