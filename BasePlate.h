//
// Created by ppxjd3 on 14/07/2021.
//

#ifndef INC_3DMOLECULARDYNAMICS_BASEPLATE_H
#define INC_3DMOLECULARDYNAMICS_BASEPLATE_H


class BasePlate {
public:
    BasePlate(double z, double A, double T): _z(z), _A(A), _T(T){}

private:
    double _z;
    double _vz{0};
    double _A;
    double _T;

};


#endif //INC_3DMOLECULARDYNAMICS_BASEPLATE_H
