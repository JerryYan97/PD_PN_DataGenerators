//
// Created by jiaruiyan on 1/19/21.
//

#ifndef PD_PN_GENERATORS_PNSIMULATOR_H
#define PD_PN_GENERATORS_PNSIMULATOR_H
#include "DataGenerator.h"

class PNSimulator : public DataGenerator {
private:

public:
    PNSimulator(TestCaseInfo &info, std::shared_ptr<ForceField> ff) : DataGenerator(info, ff) {}
    virtual void step();
};


#endif //PD_PN_GENERATORS_PNSIMULATOR_H
