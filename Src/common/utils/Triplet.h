//
// Created by jiaruiyan on 1/20/21.
//
// https://github.com/ipc-sim/IPC/blob/master/src/Utils/Triplet.h

#ifndef PD_PN_GENERATORS_TRIPLET_H
#define PD_PN_GENERATORS_TRIPLET_H
#include <iostream>
class Triplet {
public:
    int key[3];

    Triplet() = default;

    Triplet(const int *p_key) {
        key[0] = p_key[0];
        key[1] = p_key[1];
        key[2] = p_key[2];
    }

    Triplet(int key0, int key1, int key2) {
        key[0] = key0;
        key[1] = key1;
        key[2] = key2;
    }

    ~Triplet() = default;

    bool operator<(const Triplet &right) const {
        if (key[0] < right.key[0]) {
            return true;
        } else if (key[0] == right.key[0]) {
            if (key[1] < right.key[1]) {
                return true;
            } else if (key[1] == right.key[1]) {
                if (key[2] < right.key[2]) {
                    return true;
                }
            }
        }
        return false;
    }

    bool operator==(const Triplet& obj) const{
        if (key[0] == obj.key[0] && key[1] == obj.key[1] && key[2] == obj.key[2]){
            return true;
        }else{
            return false;
        }
    }
};

inline std::ostream& operator<<(std::ostream& os, const Triplet& t)
{
    os << t.key[0] << ' ' << t.key[1] << ' ' << t.key[2];
    return os;
}

class TripletHashFunc{
public:
    size_t operator()(const Triplet& t) const
    {
        return t.key[0] + t.key[1] + t.key[2];
    }
};

#endif //PD_PN_GENERATORS_TRIPLET_H
