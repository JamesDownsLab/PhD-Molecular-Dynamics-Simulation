//
// Created by ppxjd3 on 14/07/2021.
//

#ifndef INC_3DMOLECULARDYNAMICS_PARTICLE_H
#define INC_3DMOLECULARDYNAMICS_PARTICLE_H

#include <iostream>
#include "Vector.h"
#include "BasePlate.h"

inline double normalize(double dx, double L) {
    while (dx < -L / 2) dx += L;
    while (dx >= L / 2) dx -= L;
    return dx;
}

class Particle {

    friend std::istream& operator >> (std::istream& is, Particle& p);
    friend void force(Particle& p1, Particle& p2, double lx, double ly, double lz);
    friend void force(Particle& p, BasePlate& basePlate, double force_constant);


public:
    Particle() : rtd0(null), rtd1(null), rtd2(null), rtd3(null), rtd4(null) {}
    Particle(const Particle& rhs) = default;
    Particle(Particle && rhs) = default;
    Particle& operator=(const Particle & rhs) = default;
    Particle& operator=(Particle && rhs) = default;
    virtual ~Particle() = default;

    ///////////////////////////////////////
    /// Getters
    ///////////////////////////////////////
    double& r() { return _r; }
    double r() const { return _r; }
    double m() const { return _m; }

    double& x() { return rtd0.x(); }
    double x() const { return rtd0.x(); }
    double& y() { return rtd0.y(); }
    double y() const { return rtd0.y(); }
    double& z() { return rtd0.z(); }
    double z() const { return rtd0.z(); }
    double& vx() { return rtd1.x(); }
    double vx() const { return rtd1.x(); }
    double& vy() { return rtd1.y(); }
    double vy() const { return rtd1.y(); }
    double& vz() { return rtd1.z(); }
    double vz() const { return rtd1.z(); }
    Vector& force() {return _force;}
    Vector force() const {return _force;}

    ///////////////////////////////////////
    /// Setters
    ///////////////////////////////////////
    void add_force(const Vector& f) {_force += f;}
    void periodic_bc(double x_0, double y_0, double lx, double ly);
    void set_force_to_zero() { _force = null;}

    void velocity_verlet(double dt, Vector G);

private:
    Vector rtd0{null}, rtd1{null}, rtd2{null}, rtd3{null}, rtd4{null};
    Vector rtd0_old, rtd0_new;
    Vector _force;
    double _r, _m, _youngs_modulus, _poisson, _coeff_res, _coeff_fric, _damping_constant, _force_constant;
};


#endif //INC_3DMOLECULARDYNAMICS_PARTICLE_H
