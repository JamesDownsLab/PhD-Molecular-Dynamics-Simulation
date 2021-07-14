//
// Created by ppxjd3 on 14/07/2021.
//

#ifndef INC_3DMOLECULARDYNAMICS_SETUP_H
#define INC_3DMOLECULARDYNAMICS_SETUP_H

#include <random>

struct SystemProperties {
    double timestep;
    double lx;
    double ly;
    double lz;
    double dimple_rad;
    double dimple_spacing;
    double dimple_depth;
    double base_height;
};

struct MaterialProperties {
    double radius;
    double mass;
    double youngs_modulus;
    double poisson;
    double coeff_res;
    double coeff_fric;
};

MaterialProperties ballProps {
    2e-3,
    0.0387e-3,
    4e-3,
    0.48,
    0.1,
    0.5
};

void dump(
        std::ofstream& os,
        double x,
        double y,
        double z,
        MaterialProperties& props){
    os.precision(8);
    os << x << "\t" << y << "\t" << z << "\t"
        << props.radius << "\t" << props.mass << "\t"
        << props.youngs_modulus << "\t" << props.poisson << "\t"
        << props.coeff_res << "\t" << props.coeff_fric << "\n";
}

void dump_preamble(std::ofstream& os, SystemProperties& p) {
    os.precision(10);
    os << "#timestep: " << p.timestep << "\n"
        << "#lx: " << p.lx << "\n"
        << "#ly: " << p.ly << "\n"
        << "#lz: " << p.lz << "\n"
        << "#dimple_rad: " << p.dimple_rad << "\n"
        << "#dimple_spacing: " << p.dimple_spacing << "\n"
        << "#dimple_depth: " << p.dimple_depth << "\n"
        << "#base_height: " << p.base_height << "\n";
}

std::vector<std::pair<double, double>> create_lattice(double lx, double ly, double L){
    double dx = L;
    double dy = L*sqrt(3.0)/2.0;
    int Nx = floor(lx / dx);
    int Ny = floor(ly / dy);
    double x, y;
    std::vector<std::pair<double, double>> points (Nx*Ny);
    int k = 0;
    for (int i{0}; i < Nx; i++){
        for (int j{0}; j < Ny; j++){
            x = double(i)*dx + double(j%2)*dx/2.0;
            y = double(j)*dy;
            points[k] = {x, y};
            k += 1;
        }
    }
    return points;
}

std::vector<std::pair<double, double>> set_area_fraction(auto& points, double area_fraction){
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<> distr(0.0, 1.0);
    std::vector<std::pair<double, double>> new_points;
    size_t L = points.size();
    for (size_t i = 0; i<L; i++){
        double rnd = distr(eng);
        if (rnd < area_fraction) {
            new_points.push_back(points[i]);
        }
    }
    return new_points;
}

void setup_experiment(double area_fraction) {
    std::ofstream fout("initial.data");
    SystemProperties props {
        1e-6,
        0.3,
        0.3,
        0.1,
        1e-3,
        4.8e-3,
        2e-4,
        1e-3
    };

    dump_preamble(fout, props);
    /// Fill the system
    std::vector<std::pair<double, double>> xy_points = create_lattice(props.lx, props.ly, 4e-3);

    double number_fraction = area_fraction * xy_points.size()*ballProps.radius*ballProps.radius*3.141/(props.lx*props.ly);

    std::vector<std::pair<double, double>> new_xy_points = set_area_fraction(xy_points, number_fraction);

    double z = 4e-3;
    for (std::pair<double, double> p: new_xy_points){
        dump(fout, p.first, p.second, z, ballProps);
    }

}


#endif //INC_3DMOLECULARDYNAMICS_SETUP_H
