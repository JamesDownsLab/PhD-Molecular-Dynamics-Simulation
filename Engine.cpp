//
// Created by ppxjd3 on 14/07/2021.
//

#include "Engine.h"

Engine::Engine(const char *fname, ProgramOptions options) : _options{options} {
    f1 = fopen(_options.savepath.string().c_str(), "w");
    init_system(fname);

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
