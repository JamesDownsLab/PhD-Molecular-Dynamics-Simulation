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
#include <chrono>
#include "Options.h"
#include "nanoflann.h"
#include "KDTreeVectorOfVectorsAdaptor.h"

typedef  std::vector<std::vector<double>> my_vector_of_vectors_t;
typedef KDTreeVectorOfVectorsAdaptor<my_vector_of_vectors_t, double> my_kd_tree_t;

namespace fs = std::filesystem;

const double SQRT3 = sqrt(3);

class Engine {

public:
    /**
     * Initialises the Engine class
     *
     * \param fname Filename of initialisation data
     * \param options Struct containing various options for the program
     */
     Engine(Options& options);

     ///Iterates the simulation by one timestep
     void step();

     void set_baseplate(double A, double T){basePlate.set_A(A); basePlate.set_T(T);}

private:
    /// Setup the system

//    void read_system_params(const char* fname);

    void init_system();

    void add_particles();

    void add_base_particles();

    void create_dimples();

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
    std::vector<std::vector<int>> partners_base;

    /// Updates pindex and partners.
    void make_ilist();

    /// Checks if pindex will change.
    bool ilist_needs_update();
    void clear_pindex();
    void init_lattice_algorithm();
    void init_lattice_algorithm_for_base_particles();
    void clear_pindex_base();
    double rmin{0}, rmax{0}, gk{0};
    double r_base{0}, gk_base{0};
    int gm{0}, Nx{0}, Ny{0};
    int gm_base{0}, nx_base{0}, ny_base{0};

    std::vector<std::vector<int>> pindex_base;

    void dump();
    void check_dump();
    int save{1};

    ///////////////////////////////////////////////////////////
    /// Particle data
    //////////////////////////////////////////////////////////
    size_t no_of_particles{0};
    std::vector<Particle> particles;
    std::vector<Particle> base_particles;
//    std::unique_ptr<my_kd_tree_t> tree;
    Options _options;
    my_vector_of_vectors_t base_particles_for_tree;
    size_t no_of_base_particles{0};

    /////////////////////////////////////////////////////////
    /// File Writing
    /////////////////////////////////////////////////////////
    std::FILE* f1;

    double timestep{0};

    double Time{0};
    Vector G{0, 0, -9.81};
    BasePlate basePlate{0, 0, 0};
    double lx;
    double ly;
    double lz;

    std::chrono::steady_clock::time_point begin;

};


#endif //INC_3DMOLECULARDYNAMICS_ENGINE_H
