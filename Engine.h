//
// Created by ppxjd3 on 14/07/2021.
//

#ifndef INC_3DMOLECULARDYNAMICS_ENGINE_H
#define INC_3DMOLECULARDYNAMICS_ENGINE_H

#include <cmath>
#include <fstream>
#include <filesystem>
#include <vector>
#include "Particle.h"

namespace fs = std::filesystem;

const double SQRT3 = sqrt(3);

struct ProgramOptions {
    fs::path savepath;
};

class Engine {

public:
    /**
     * Initialises the Engine class
     *
     * \param fname Filename of initialisation data
     * \param options Struct containing various options for the program
     */
     Engine(const char* fname, ProgramOptions options);

     ///Iterates the simulation by one timestep
     void step();

private:
    /// Setup the system
    void init_system(const char* fname);

    /// Calculate the collisional forces between all particles
    void make_forces();

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
    double rmin, rmax, gk;
    int gm, Nx, Ny;

    void dump();

    ///////////////////////////////////////////////////////////
    /// Particle data
    //////////////////////////////////////////////////////////
    size_t no_of_particles;
    std::vector<Particle> particles;
    ProgramOptions _options;

    /////////////////////////////////////////////////////////
    /// File Writing
    /////////////////////////////////////////////////////////
    std::FILE* f1;

    double timestep;
    double lx;
    double ly;
    double lz;
    double x_0{0.0};
    double y_0{0.0};
    double z_0{0.0};
    double dimple_rad;
    double dimple_spacing;
    double dimple_depth;
    double base_height;
    double Time{0};
};


#endif //INC_3DMOLECULARDYNAMICS_ENGINE_H
