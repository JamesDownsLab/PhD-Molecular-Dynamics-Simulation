//
// Created by ppxjd3 on 14/07/2021.
//

#include "Particle.h"
#include <cmath>

std::istream &operator>>(std::istream &is, Particle &p) {
    is >> p.rtd0 >> p._r >> p._m >> p._youngs_modulus >> p._poisson >> p._coeff_res >> p._coeff_fric >> p._damping_constant;

    p._force_constant = 2*p._youngs_modulus* sqrt(p._r)/(3*(1-p._poisson*p._poisson));
    return is;
}

void Particle::velocity_verlet(double dt, Vector G) {

        double a1 = dt;
        double a2 = dt * dt / 2;

        rtd0 += rtd1 * a1 + rtd2 * a2;
        Vector new_rtd2 = (_force) * (1 / _m) + G;
        rtd1 += (rtd2 + new_rtd2) * (dt / 2);
        rtd2 = new_rtd2;
}

void force(Particle &p1, Particle &p2, double lx, double ly, double lz) {
    double dx = normalize(p1.x() - p2.x(), lx);
    double dy = normalize(p1.y() - p2.y(), ly);
    double dz = normalize(p1.z() - p2.z(), lz);
    if (std::abs(dx) < p1.r() + p2.r() && std::abs(dy) < p1.r() + p2.r() && std::abs(dz) < p1.r() + p2.r()){
        double rr = sqrt(dx*dx + dy*dy + dz*dz);
        double r1 = p1.r();
        double r2 = p2.r();

        // Overlap
        double xi = r1 + r2 - rr;

        if (xi > 0) { // If overlapping
            double sqrt_xi = sqrt(xi);
            double rr_rez = 1 / rr;

            // Unit vectors
            double ex = dx * rr_rez;
            double ey = dy * rr_rez;
            double ez = dz * rr_rez;

            // Relative velocities
            double dvx = p1.vx() - p2.vx();
            double dvy = p1.vy() - p2.vy();
            double dvz = p1.vz() - p2.vz();

            // Overlap rate
            double xidot = -(ex * dvx + ey * dvy + ez*dvz);

            double elastic_force = p1._force_constant * xi * sqrt_xi;
            double dissipative_force = p1._force_constant * p1._damping_constant * sqrt_xi * xidot;
            double fn = elastic_force + dissipative_force;

            if (fn < 0) fn = 0;
            p1.add_force(Vector(fn * ex, fn*ey, fn*ez));
            p1.add_force(Vector(-fn*ex, -fn*ey, -fn*ez));


        }
    }
}

void Particle::periodic_bc(double x_0, double y_0, double lx, double ly) {
    while (rtd0.x() < x_0) rtd0.x() += lx;
    while (rtd0.x() > x_0 + lx) rtd0.x() -= lx;
    while (rtd0.y() < y_0) rtd0.y() += ly;
    while (rtd0.y() > y_0 + ly) rtd0.y() -= ly;
}

void force(Particle &p, BasePlate &basePlate, double force_constant) {
    double dz = p.z() - basePlate.z();

    // Overlap
    double xi = p.r() - dz;
    if (xi > 0){ // Overlapping
        double xidot = -(p.vz() - basePlate.vz());
        double elastic_force = force_constant * pow(xi, 1.5);
        double dissipative_force = force_constant * p._damping_constant*xidot*sqrt(xi)/2;
        double fn = elastic_force + dissipative_force;

        if (fn > 0) {
            p.add_force(Vector(0, 0, fn));
        }
    }
}
