//
// Created by ppxjd3 on 14/07/2021.
//

#include "Particle.h"
#include <cmath>


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

            double Y = p1._youngs_modulus;
            double poisson = p1._poisson;
            double force_constant = 2*Y*sqrt(r1)/(3*(1-poisson*poisson));

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

            double elastic_force = force_constant * xi * sqrt_xi;
            double dissipative_force = force_constant * p1._damping_constant * sqrt_xi * xidot;
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

void force(Particle &p, Particle &b, BasePlate &basePlate) {
    double dx = p.x() - b.x();
    double dy = p.y() - b.y();
    double dz = p.z() - basePlate.z();
    if (std::abs(dx) < p.r() + b.r() && std::abs(dy) < p.r() + b.r() && std::abs(dz) < p.r() + b.r()) {
        double rr = sqrt(dx*dx + dy*dy + dz*dz);
        double r1 = p.r();
        double r2 = b.r();

        // Overlap
        double xi = r1 + r2 - rr;
        if (xi > 0) {
            double force_constant = 4*sqrt(p.r())/3*(1/((1-p._poisson*p._poisson)/p._youngs_modulus + (1-b._poisson*b._poisson)/b._youngs_modulus));
            double damping_constant = (p._damping_constant + b._damping_constant)/2;

            double sqrt_xi = sqrt(xi);
            double rr_rez = 1 / rr;

            // Unit vectors
            double ex = dx * rr_rez;
            double ey = dy * rr_rez;
            double ez = dz * rr_rez;

            // Relative velocities
            double dvx = p.vx();
            double dvy = p.vy();
            double dvz = p.vz() - basePlate.vz();

            // Overlap rate
            double xidot = -(ex * dvx + ey * dvy + ez*dvz);

            double elastic_force = force_constant * sqrt_xi * xi;
            double dissipative_force = force_constant * damping_constant * xidot * sqrt_xi / 2;
            double fn = elastic_force + dissipative_force;


            if (fn > 0) {
                p.add_force(Vector(ex * fn, ey * fn, ez * fn));
            }
        }
    }
}

void Particle::predict(double dt) {
    double a1 = dt;
    double a2 = a1*dt/2;
    double a3 = a2*dt/3;
    double a4 = a3*dt/4;

    rtd0 += a1*rtd1 + a2*rtd2 + a3*rtd3 + a4*rtd4;
    rtd1 += a1*rtd2 + a2*rtd3 + a3*rtd4;
    rtd2 += a1*rtd3 + a2*rtd4;
    rtd3 += a1*rtd4;
}

void Particle::correct(double dt, Vector G, double mass) {

    static Vector accel, corr;

    double dtrez = 1/dt;
    const double coeff0 = double(19)/double(90) * (dt*dt/double(2));
    const double coeff1 = double(3)/double(4)*(dt/double(2));
    const double coeff3 = double(1)/double(2)*(double(3)*dtrez);
    const double coeff4 = double(1)/double(12)*(double(12)*(dtrez*dtrez));

    accel = Vector((1/mass)*_force.x()+G.x(), (1/mass)*_force.y()+G.y(), (1/mass)*_force.z()+G.z());
    corr = accel - rtd2;
    rtd0 += coeff0*corr;
    rtd1 += coeff1*corr;
    rtd2 = accel;
    rtd3 += coeff3*corr;
    rtd4 += coeff4*corr;
}
