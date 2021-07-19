//
// Created by ppxjd3 on 14/07/2021.
//

#ifndef INC_3DMOLECULARDYNAMICS_ENGINE_H
#define INC_3DMOLECULARDYNAMICS_ENGINE_H

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <vector>
#include "Particle.h"
#include "BasePlate.h"
#include <random>
#include "nanoflann.h"
#include "KDTreeVectorOfVectorsAdaptor.h"

typedef  std::vector<std::vector<double>> my_vector_of_vectors_t;
typedef KDTreeVectorOfVectorsAdaptor<my_vector_of_vectors_t, double> my_kd_tree_t;

namespace fs = std::filesystem;

const double SQRT3 = sqrt(3);

struct ProgramOptions {
    fs::path savepath;
    int save_interval;
};

class Engine {

public:
    /**
     * Initialises the Engine class
     *
     * \param fname Filename of initialisation data
     * \param options Struct containing various options for the program
     */
     Engine(const char* fname, ProgramOptions& options);

     ///Iterates the simulation by one timestep
     void step();

     void set_baseplate(double A, double T){basePlate.set_A(A); basePlate.set_T(T);}

private:
    /// Setup the system

    void read_system_params(const char* fname);

    void init_system();

    void add_particles();

    void add_base_particles();

    void make_tree();

    /// Calculate the collisional forces between all particles
    void make_forces();

    /// Calculate the forces with the base
    void make_plate_forces();

    /// Calculate forces and updates positions/velocities
    void integrate();

    /////////////////////////////////////////////////////////////////////////////
    /// Lattice
    /////////////////////////////////////////////////////////////////////////////

    /// 2D vector representing the lattice cells
    /// Value contains -1 if empty or the index of the particle.
    std::vector<std::vector<int>> pindex;

    /// Vector containing a list of neighbours for each particle.
    std::vector<std::vector<int>> partners;

    /// Updates pindex and partners.
    void make_ilist();

    /// Checks if pindex will change.
    bool ilist_needs_update();
    void clear_pindex();
    void init_lattice_algorithm();
    double rmin{0}, rmax{0}, gk{0};
    int gm{0}, Nx{0}, Ny{0};

    void dump();
    void check_dump();
    int save{1};

    ///////////////////////////////////////////////////////////
    /// Particle data
    //////////////////////////////////////////////////////////
    size_t no_of_particles{0};
    std::vector<Particle> particles;
    std::vector<Particle> base_particles;
    std::unique_ptr<my_kd_tree_t> tree;
    ProgramOptions _options;
    my_vector_of_vectors_t base_particles_for_tree;
    size_t no_of_base_particles{0};

    /////////////////////////////////////////////////////////
    /// File Writing
    /////////////////////////////////////////////////////////
    std::FILE* f1;

    double timestep{0};
    double lx{0};
    double ly{0};
    double lz{0};
    double x_0{0.0};
    double y_0{0.0};
    double z_0{0.0};
    double dimple_rad{0};
    double dimple_spacing{0};
    double dimple_depth{0};
    double base_height{0};
    double ball_rad{0};
    double base_rad{0};
    double ball_height{0};
    double ball_youngs{0};
    double ball_poisson{0};
    double base_youngs{0};
    double base_poisson{0};
    double ball_base_normal_constant{0};
    double area_fraction{0};


    double _base_amplitude{0};
    double _base_period{0};
    double Time{0};
    Vector G{0, 0, -9.81};
    BasePlate basePlate{0, 0, 0};

};


#endif //INC_3DMOLECULARDYNAMICS_ENGINE_H
