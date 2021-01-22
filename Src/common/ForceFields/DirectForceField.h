//
// Created by jiaruiyan on 1/20/21.
//

#ifndef PD_PN_GENERATORS_DIRECTFORCEFIELD_H
#define PD_PN_GENERATORS_DIRECTFORCEFIELD_H

#include "ForceField.h"
class DirectForceField : public ForceField{

private:
    Eigen::Vector3d m_force;

public:
    DirectForceField() : m_force(Eigen::Vector3d(0.0, 0.0, 0.0)) {}
    void SetForceField(ForceFieldInfo info){
        m_force = info.dir_force;
    }
    virtual Eigen::Vector3d GetForce(Eigen::Vector3d pos){
        return m_force;
    }
};

#endif //PD_PN_GENERATORS_DIRECTFORCEFIELD_H
