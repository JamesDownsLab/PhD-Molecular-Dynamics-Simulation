//
// Created by ppxjd3 on 14/07/2021.
//

#include "Engine.h"

Engine::Engine(const char *fname, ProgramOptions& options)
        : _options{options} {
    f1 = fopen(_options.savepath.string().c_str(), "w");
    init_system(fname);
    init_lattice_algorithm();
    basePlate.set_zi(base_height);
}

void Engine::init_system(const char *fname) {
    std::ifstream fparticle{fname};

    // Read the system properties
    while (fparticle.peek() == '#'){
        std::string type;
        fparticle >> type;

        if (type == "#timestep:"){
            fparticle >> timestep;
        }
        else if (type == "#lx:"){
            fparticle >> lx;
        }
        else if (type == "#ly:"){
            fparticle >> ly;
        }
        else if (type == "#lz:"){
            fparticle >> lz;
        }
        else if (type == "#dimple_rad:"){
            fparticle >> dimple_rad;
        }
        else if (type == "#dimple_spacing:"){
            fparticle >> dimple_spacing;
        }
        else if (type == "#dimple_depth:"){
            fparticle >> dimple_depth;
        }
        else if (type == "#base_height:"){
            fparticle >> base_height;
        }
        else {
            std::cout << "Unknown type: " << type << std::endl;
        }
        fparticle.ignore(100, '\n');
    }

    // Read the particles
    while (fparticle){
        Particle pp;
        fparticle >> pp;
        if (fparticle) {
            particles.push_back(pp);
        }
    }

    no_of_particles = particles.size();
    std::cout << no_of_particles << " particles read" << std::endl;

    dump();
}

void Engine::dump() {
    std::fprintf(f1, "ITEM: TIMESTEP\n%d\n", int(Time / timestep));
    std::fprintf(f1, "ITEM: TIME\n%.8f\n", Time);
    std::fprintf(f1, "ITEM: BOX BOUNDS pp pp f\n%.4f %.4f\n%.4f %.4f\n%.4f %.4f\n", 0.0, lx, 0.0, ly, 0.0, lz);
    std::fprintf(f1, "ITEM: NUMBER OF ATOMS\n%d\n", int(no_of_particles));
    std::fprintf(f1, "ITEM: ATOMS x y z vx vy vz radius type\n");
    for (Particle& p : particles) {
        std::fprintf(f1, "%.9f %.9f %.9f %.9f %.9f %.9f %.9f %d\n", p.x(), p.y(), p.z(), p.vx(), p.vy(), p.vz(), p.r(), 0);
    }
}

void Engine::step() {
     // Check whether the optimiser needs updating
     if (ilist_needs_update()) {make_ilist();}

     basePlate.update(Time);

     integrate();

     dump();
}

void Engine::integrate() {
    // Set forces to zero
    std::for_each(particles.begin(), particles.end(),
                  [&](Particle& p){
        p.set_force_to_zero();
    });

    // Calculate all the forces between particles
    make_forces();

    // Calculate all the forces between the particles and the plate
    make_plate_forces();

    // Update  the positions of all the particles
    std::for_each(particles.begin(), particles.end(),
                  [&](Particle& p){
        p.velocity_verlet(timestep, G);
    });

    // Apply periodic boundary conditions
    std::for_each(particles.begin(), particles.end(),
                  [&](Particle& p) {p.periodic_bc(x_0, y_0, lx, ly);});

    Time += timestep;
}




///////////////////////////////////////////////////////////////////////////////
/// Lattice Method
///////////////////////////////////////////////////////////////////////////////

void Engine::init_lattice_algorithm() {
    rmin = particles[0].r();
    rmax = particles[0].r();

    // Calculate gk, the size of the lattice sites
    gk = sqrt(2) * rmin;

    // Calcualte gm, the number of lattice sites that the largest particle covers
    gm = int(2*rmax/gk) + 1;

    // calculate the numebr of lattice sites in each dimension
    Nx = int(lx / gk) + 1;
    Ny = int(ly / gk + 1);
    partners.resize(no_of_particles);
    pindex.resize(Nx);
    for (auto& p : pindex){
        p.resize(Ny);
    }
    clear_pindex();
    make_ilist();
}

void Engine::make_ilist() {
    // For each particle, add it to the pindex lattice site
    // and clear its partners
    for (unsigned int i{0}; i<no_of_particles; i++){
        double x = particles[i].x();
        double y = particles[i].y();
        int ix = int((x-x_0) / gk);
        int iy = int((y-y_0) / gk);
        pindex[ix][iy] = i;
        partners[i].clear();
    }

    // Generate the partners list for each particle
    for (unsigned int i{ 0 }; i < no_of_particles; i++) {
        double x = particles[i].x();
        double y = particles[i].y();
        if ((x >= x_0) && (x < x_0 + lx) && (y >= y_0) && (y < y_0 + ly)) {
            int ix = int((x - x_0) / gk);
            int iy = int((y - y_0) / gk);
            // Check the adjacent gm lattice sites for particles
            for (int dx = -gm; dx <= gm; dx++) {
                for (int dy = -gm; dy <= gm; dy++) {
                    int iix = (ix + dx + Nx) % Nx;
                    int iiy = (iy + dy + Ny) % Ny;
                    int k = pindex[iix][iiy];
                    // Only record the particle once
                    if (k > (int)i) {
                        partners[i].push_back(k);
                    }
                }
            }
        }
    }
}

bool Engine::ilist_needs_update() {
    // If the pindex is the same as last time then it doesn't need updating
    for (unsigned int i{ 0 }; i < no_of_particles; i++) {
        double x = particles[i].x();
        double y = particles[i].y();
        int ix = int((x - x_0) / gk);
        int iy = int((y - y_0) / gk);
        if (pindex[ix][iy] != i) {
            clear_pindex();
            return true;
        }
    }
    return false;
}

void Engine::clear_pindex()
{
    for (auto& p : pindex) {
        for (auto& q : p) {
            q = -1;
        }
    }
}

void Engine::make_forces() {
    // Loop over the partners list for each particle
    for (unsigned int i{ 0 }; i < no_of_particles; i++) {
        for (unsigned int k{ 0 }; k < partners[i].size(); k++) {
            int pk = partners[i][k];
            force(particles[i], particles[pk], lx, ly, lz);
        }
    }
}

void Engine::make_plate_forces() {
    for (unsigned  int i{0}; i < no_of_particles; i++){
        force(particles[i], basePlate);
    }
}


