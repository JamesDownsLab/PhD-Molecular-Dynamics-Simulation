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

        if (xi > 1e-10) { // If overlapping



            double Y = p1._youngs_modulus;
            double poisson = p1._poisson;
            double force_constant = 2*Y*sqrt(r1)/(3*(1-poisson*poisson));

            double sqrt_xi = sqrt(xi);
            double rr_rez = 1 / rr;

            // Unit vectors
            double ex = dx * rr_rez;
            double ey = dy * rr_rez;
            double ez = dz * rr_rez;

            double dvx = p1.vx() - p2.vx();
            double dvy = p1.vy() - p2.vy();
            double dvz = p1.vz() - p2.vz();

            Eigen::Vector3d n = {ex, ey, ez};

            Eigen::Vector3d v1 = p1.rtd1;
            Eigen::Vector3d v2 = p2.rtd1;
            Eigen::Vector3d omega1 = p1.rot1;
            Eigen::Vector3d omega2 = p2.rot2;
            Eigen::Vector3d vrel = v1 - v2 - (r1*omega1 + r2*omega2).cross(n);
            Eigen::Vector3d vtrel = vrel - vrel.dot(n)*n;
            Eigen::Vector3d t = vtrel * (1/vtrel.norm());

            double xidot = -(ex * dvx + ey * dvy + ez*dvz);

            double gamma = p1._tangential_damping;

            // Normal Forces
            double elastic_force = force_constant * xi * sqrt_xi;
            double dissipative_force = force_constant * p1._damping_constant * sqrt_xi * xidot;
            double fn = elastic_force + dissipative_force;
            if (fn < 0) fn = 0;

            // Tangential forces
            double mu = p1._friction;
            double ft = -gamma * xidot;
            if (ft < -mu*fn) ft = -mu*fn;
            if (ft>mu*fn) ft = mu*fn;

            // Total force
            Eigen::Vector3d force = fn*n + ft*t;

            Eigen::Vector3d torque = force.cross(n);

            p1.add_force(force);
            p2.add_force(-1*force);

            p1.add_torque(torque);
            p2.add_torque(-1*torque);


        }
    }
}


void Particle::periodic_bc(double x_0, double y_0, double lx, double ly) {
    while (rtd0.x() < x_0) rtd0.x() += lx;
    while (rtd0.x() > x_0 + lx) rtd0.x() -= lx;
    while (rtd0.y() < y_0) rtd0.y() += ly;
    while (rtd0.y() > y_0 + ly) rtd0.y() -= ly;
}

bool force(Particle &p, Particle &b, BasePlate &basePlate, double timestep) {
    double dx = p.x() - b.x();
    double dy = p.y() - b.y();
    double dz = p.z() - (basePlate.z()+b.z());
    if (std::abs(dx) < p.r() + b.r() && std::abs(dy) < p.r() + b.r() && std::abs(dz) < p.r() + b.r()) {
        double rr = sqrt(dx*dx + dy*dy + dz*dz);
        double r1 = p.r();
        double r2 = b.r();

        // Overlap
        double xi = r1 + r2 - rr;
        if (xi > 1e-10) {



            double Y = (p._youngs_modulus*b._youngs_modulus)/(p._youngs_modulus+b._youngs_modulus);
            double poisson = 0.5*(p._poisson+b._poisson);
            double A = 0.5*(p._damping_constant + b._damping_constant);
            double force_constant = 2*Y*sqrt(r1)/(3*(1-poisson*poisson));
            double sqrt_xi = sqrt(xi);

            // Unit vectors
            double rr_rez = 1/rr;
            double ex = dx * rr_rez;
            double ey = dy * rr_rez;
            double ez = dz * rr_rez;

            // Relative velocities
            double dvx = p.vx();
            double dvy = p.vy();
            double dvz = p.vz() - basePlate.vz();

            Eigen::Vector3d n = {ex, ey, ez};

            Eigen::Vector3d v1 = p.rtd1;
            Eigen::Vector3d v2 = {0, 0, basePlate.vz()};
            Eigen::Vector3d omega1 = p.rot1;
            Eigen::Vector3d omega2 = {0, 0, 0};
            Eigen::Vector3d vrel = v1 - v2 -  (r1*omega1).cross(n);
            Eigen::Vector3d vtrel = vrel - vrel.dot(n)*n;

            // Update the contacts
            if (p.base_contacts.find(b.index) == p.base_contacts.end()){
                // No contact
                p.base_contacts[b.index] = vtrel*timestep;
            }
            else {
                p.base_contacts[b.index] += vtrel*timestep;
            }


            double vtrel_size = vtrel.norm();
            Eigen::Vector3d t;
            if (vtrel_size > 0){
                t = vtrel * (1/vtrel_size);
            }
            else{
                t = {0, 0, 0};
            }

            double xidot = -(ex * dvx + ey * dvy + ez*dvz);
            double gamma = 0.5*(p._tangential_damping+b._tangential_damping);

            // Normal forces
            double elastic_force = force_constant * xi * sqrt_xi;
            double dissipative_force = force_constant * A * sqrt_xi * xidot;
            double fn = elastic_force + dissipative_force;
            if (fn < 0) fn = 0;

            // Tangential forces
            double mu = p._friction;
            double elongation = p.base_contacts[b.index].norm();
            double ft = -gamma * elongation;
            if (ft < -mu*fn) ft = -mu*fn;
            if (ft > mu*fn) ft = mu*fn;

            // Total force
            Eigen::Vector3d force = fn*n + ft*t;
            Eigen::Vector3d torque = force.cross(n);
            p.add_force(force);
            p.add_torque(torque);
            return true;
        }
        else{
            return false;
        }
    }
}

void Particle::predict(double dt) {
    double a1 = dt;
    double a2 = a1*dt/2;
    double a3 = a2*dt/3;

    rtd0 += a1*rtd1 + a2*rtd2 + a3*rtd3;
    rtd1 += a1*rtd2 + a2*rtd3;
    rtd2 += a1*rtd3;

    rot0 += a1*rot1 + a2*rot2 + a3*rot3;
    rot1 += a1*rot2 + a2*rot3;
    rot2 += a1*rot3;
}

void Particle::correct(double dt, Eigen::Vector3d G) {

    static Eigen::Vector3d accel, corr, rot_accel, rot_corr;

    double dtrez = 1/dt;
    const double coeff0 = double(1)/double(6) * (dt*dt/double(2));
    const double coeff1 = double(5)/double(6)*(dt/double(2));
    const double coeff3 = double(1)/double(3)*(double(3)*dtrez);

    accel = Eigen::Vector3d((1/_m)*_force.x()+G.x(), (1/_m)*_force.y()+G.y(), (1/_m)*_force.z()+G.z());
    corr = accel - rtd2;
    rtd0 += coeff0*corr;
    rtd1 += coeff1*corr;
    rtd2 = accel;
    rtd3 += coeff3*corr;

    rot_accel = _torque * (1/J);
    rot_corr = rot_accel - rot2;
    rot0 += coeff0*rot_corr;
    rot1 += coeff1*rot_corr;
    rot2 = rot_accel;
    rot3 += coeff3*rot_corr;
}

void Particle::update_base_contacts(std::set<size_t>& contacts){
    std::set<size_t> to_delete;
    for (const auto& [key, value] : base_contacts){
        if (!contacts.contains(key)){
            to_delete.insert(key);
        }
    }
    for (size_t k: to_delete){
        base_contacts.erase(k);
    }
}
