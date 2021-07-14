//
// Created by ppxjd3 on 14/07/2021.
//

#ifndef INC_3DMOLECULARDYNAMICS_BASEPLATE_H
#define INC_3DMOLECULARDYNAMICS_BASEPLATE_H

#include <cmath>


class BasePlate {
public:
    BasePlate(double z, double A, double T): _z(z), _A(A), _T(T){}

    void update(double Time);

    void set_zi(double z) {_z = z;}
    void set_A(double A) {_A = A;}
    void set_T(double T) {_T = T; _omega=2*M_PI/T;}

    double& z() {return _z;}
    double z() const {return _z;}

    double& vz() {return _vz;}
    double vz() const {return _vz;}


private:
    double _z{0};
    double _vz{0};
    double _A{0};
    double _T{0};
    double _omega{0};

};


#endif //INC_3DMOLECULARDYNAMICS_BASEPLATE_H
