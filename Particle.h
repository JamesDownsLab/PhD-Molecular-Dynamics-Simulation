//
// Created by ppxjd3 on 14/07/2021.
//

#ifndef INC_3DMOLECULARDYNAMICS_PARTICLE_H
#define INC_3DMOLECULARDYNAMICS_PARTICLE_H

#include <iostream>
#include "Vector.h"
#include "BasePlate.h"
#include "Options.h"

inline double normalize(double dx, double L) {
    while (dx < -L / 2) dx += L;
    while (dx >= L / 2) dx -= L;
    return dx;
}

class Particle {
    //////////////////////
    /// Force calculation
    //////////////////////
    friend void force(Particle& p1, Particle& p2, double lx, double ly, double lz);
    friend void force(Particle& p, Particle& b, BasePlate& basePlate);


public:
    Particle(double x, double y, double z, ParticleProps props) :
        rtd0(x, y, z), rtd1(null_vec), rtd2(null_vec), _r(props.radius),
        _m(props.mass), _youngs_modulus(props.youngs_modulus), _poisson(props.poisson),
        _damping_constant(props.damping_factor){}
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
    void add_torque(const Vector & t) {_torque += t;}
    void periodic_bc(double x_0, double y_0, double lx, double ly);
    void set_force_to_zero() { _force = null_vec;}
    void set_torque_to_zero() {_torque = null_vec;}
    void set_z(double z) {rtd0 += Vector(0, 0, z);}

    ///////////////////////////////////////
    /// Integration
    ////////////////////////////////////
    void predict(double dt);
    void correct(double dt, Vector G);

private:
    Vector rtd0{null_vec}, rtd1{null_vec}, rtd2{null_vec}, rtd3{null_vec}, rtd4{null_vec};
    Vector _force{null_vec};
    Vector _torque{null_vec};
    double _r;
    double _m;
    double _youngs_modulus;
    double _poisson;
    double _damping_constant;
};
#endif //INC_3DMOLECULARDYNAMICS_PARTICLE_H
