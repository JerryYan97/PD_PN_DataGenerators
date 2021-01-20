//
// Created by jiaruiyan on 1/20/21.
//

#ifndef PD_PN_GENERATORS_FORCEFIELD_H
#define PD_PN_GENERATORS_FORCEFIELD_H

#include <Eigen/Dense>

enum ForceType{
    DIRECT_FORCE,
    CIRCLE_FORCE
};

struct ForceFieldInfo{
    ForceType type;
    Eigen::Vector3f dir_force;
};

class ForceField {
protected:
    ForceFieldInfo m_info;
public:
    virtual Eigen::Vector3f GetForce(Eigen::Vector3f pos) = 0;
};

#endif //PD_PN_GENERATORS_FORCEFIELD_H
