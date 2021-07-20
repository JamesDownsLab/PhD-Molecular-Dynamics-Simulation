//
// Created by ppxjd3 on 14/07/2021.
//

#include "Engine.h"

#include <memory>

Engine::Engine(const char *fname, ProgramOptions& options)
        : _options{options} {
    begin = std::chrono::steady_clock::now();
    f1 = fopen(_options.savepath.string().c_str(), "w");
    read_system_params("params.txt");
    init_system();

    basePlate.set_zi(base_height);
    ball_base_normal_constant =
            4*sqrt(particles[0].r())/3*
                    (1/((1-ball_poisson*ball_poisson)/ball_youngs + (1-base_poisson*base_poisson)/base_youngs));

    ball_normal_constant = 2*ball_youngs * sqrt(ball_rad)/(3*(1-ball_poisson*ball_poisson));
    ball_damping_constant = 0.01;
    ball_base_damping_constant = 0.01;
}

void Engine::init_system() {
    add_particles();
    add_base_particles();
    dump();
    init_lattice_algorithm();
    init_lattice_algorithm_for_base_particles();
}

void Engine::dump() {
    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    std::cout << "DUMP Simulation Time : " << Time << " s\t"
        << "Timestep : " << Time/timestep << "\t"
        << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(now-begin).count() << "s" << std::endl;

    std::fprintf(f1, "ITEM: TIMESTEP\n%d\n", int(Time / timestep));
    std::fprintf(f1, "ITEM: TIME\n%.8f\n", Time);
    std::fprintf(f1, "ITEM: BOX BOUNDS pp pp f\n%.4f %.4f\n%.4f %.4f\n%.4f %.4f\n", 0.0, lx, 0.0, ly, 0.0, lz);
    std::fprintf(f1, "ITEM: NUMBER OF ATOMS\n%d\n", int(no_of_particles)+int(no_of_base_particles));
    std::fprintf(f1, "ITEM: ATOMS x y z vx vy vz radius type\n");
    for (Particle& p : particles) {
        std::fprintf(f1, "%.9f %.9f %.9f %.9f %.9f %.9f %.9f %d\n", p.x(), p.y(), p.z(), p.vx(), p.vy(), p.vz(), p.r(), 0);
    }
    for (Particle& p: base_particles) {
        std::fprintf(f1, "%.9f %.9f %.9f %.9f %.9f %.9f %.9f %d\n", p.x(), p.y(), basePlate.z(), p.vx(), p.vy(), basePlate.vz(), p.r(), 1);
    }
}

void Engine::step() {
     // Check whether the optimiser needs updating
     if (ilist_needs_update()) {make_ilist();}

     basePlate.update(Time);

     integrate();

     check_dump();
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
        p.velocity_verlet(timestep, G, ball_mass);
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
    rmin = particles.at(0).r();
    rmax = particles.at(0).r();

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
        pindex[ix][iy] = (int)i;
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
        if (pindex.at(ix).at(iy) != i) {
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
            force(particles[i], particles[pk], ball_normal_constant, ball_damping_constant, lx, ly, lz);
        }
    }
}

void Engine::make_plate_forces() {
    for (auto& p: particles){
        double z = p.z();
        if (z-basePlate.z() < ball_rad+base_rad) {
            double x = p.x();
            double y = p.y();
            int ix = int(x / gk_base);
            int iy = int(y / gk_base);
            for (int dx = -gm_base; dx <= gm_base; dx++) {
                for (int dy = -gm_base; dy <= gm_base; dy++) {
                    size_t k = pindex_base[(ix + dx + nx_base) % nx_base][(iy + dy + ny_base) % ny_base];
                    force(p, base_particles.at(k), basePlate, ball_base_normal_constant, ball_base_damping_constant);
                }
            }
        }
    }
}

void Engine::check_dump() {
    if (save != _options.save_interval) {
        save++;
    }
    else {
        save = 1;
        dump();
    }
}

void Engine::read_system_params(const char *fname) {

    std::ifstream fparticle{fname};
    // Read the system properties
    while (fparticle.peek() == '#'){
        std::string type;
        fparticle >> type;

        if (type == "#timestep:"){
            fparticle >> timestep;
            std::cout << "Timestep = " << timestep << std::endl;
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
        else if (type == "#ball_youngs:"){
            fparticle >> ball_youngs;
        }
        else if (type == "#ball_poisson:"){
            fparticle >> ball_poisson;
        }
        else if (type == "#base_youngs:"){
            fparticle >> base_youngs;
        }
        else if (type == "#base_poisson:"){
            fparticle >> base_poisson;
        }
        else if(type == "#ball_rad:"){
            fparticle >> ball_rad;
        }
        else if (type == "#base_rad:"){
            fparticle >> base_rad;
        }
        else if (type == "#area_fraction:"){
            fparticle >> area_fraction;
        }
        else if (type == "#ball_height:"){
            fparticle >> ball_height;
        }
        else if (type == "#ball_mass:"){
            fparticle >> ball_mass;
        }
        else {
            std::cout << "Unknown type: " << type << std::endl;
        }
        fparticle.ignore(100, '\n');
    }

}

void Engine::add_particles() {
    double dx = 2*ball_rad;
    double dy = 2*ball_rad*sqrt(3.0)/2.0;
    int nx = floor(lx / dx);
    int ny = floor(ly / dy);

    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<> distr(0.0, 1.0);

    for (int i{0}; i < nx; i++){
        for (int j{0}; j < ny; j++){
            double x = double(i)*dx + double(j%2)*dx/2.0;
            double y = double(j)*dy;
            Particle pp(x, y, ball_height, ball_rad);
            if (distr(eng) < area_fraction) {
                particles.push_back(pp);
            }
        }
    }
    no_of_particles = particles.size();
}

void Engine::add_base_particles() {
    double dx = 2*base_rad;
    double dy = 2*base_rad*sqrt(3.0)/2.0;
    int nx = floor(lx / dx);
    int ny = floor(ly / dy);
    for (int i{0}; i < nx; i++){
        for (int j{0}; j < ny; j++){
            double x = double(i)*dx + double(j%2)*dx/2.0;
            double y = double(j)*dy;
            Particle pp(x, y, base_height, base_rad);
            base_particles.push_back(pp);
        }
    }
    no_of_base_particles = base_particles.size();
}

void Engine::init_lattice_algorithm_for_base_particles() {
    r_base = base_particles.at(0).r();
    gk_base = sqrt(2)*r_base;
    gm_base = int(2*particles.at(0).r()/gk+1);
    nx_base = int(lx/gk_base)+1;
    ny_base = int(ly/gk_base)+1;

    partners_base.resize(no_of_particles);
    pindex_base.resize(nx_base);
    for (auto & px : pindex_base){
        px.resize(ny_base);
    }
    clear_pindex_base();

    for (size_t i=0; i<base_particles.size(); i++) {
        double x = base_particles[i].x();
        double y = base_particles[i].y();
        size_t ix = int(x/gk_base);
        size_t iy = int(y/gk_base);
        pindex_base.at(ix).at(iy)=i;
    }
}

void Engine::clear_pindex_base() {
    for (auto &px: pindex) {
        for (auto &py: px) {
            py = -1;
        }
    }
}
