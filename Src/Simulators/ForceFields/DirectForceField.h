//
// Created by jiaruiyan on 1/20/21.
//

#ifndef PD_PN_GENERATORS_DIRECTFORCEFIELD_H
#define PD_PN_GENERATORS_DIRECTFORCEFIELD_H

#include "ForceField.h"
class DirectForceField : public ForceField{

private:
    Eigen::Vector3f m_force;

public:
    DirectForceField(ForceFieldInfo info) {
        m_info = info;
        m_force = m_info.dir_force;
    }

    virtual Eigen::Vector3f GetForce(Eigen::Vector3f pos){
        return m_force;
    }
};

#endif //PD_PN_GENERATORS_DIRECTFORCEFIELD_H
