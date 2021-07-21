//
// Created by ppxjd3 on 21/07/2021.
//
#include "Options.h"

Options read_input_file(const char* fname){
    std::ifstream stream{fname};
    ProgramOptions programOptions = ProgramOptions();
    SystemProps systemProps = SystemProps();
    ParticleProps ballProps = ParticleProps();
    ParticleProps baseProps = ParticleProps();
    while(stream.peek() == '#'){
        std::string type;
        stream >> type;

        if (type == "#savepath:"){
            std::string savepath;
            stream >> savepath;
            programOptions.savepath = savepath;
//            stream >> programOptions.savepath;
        }
        else if (type == "#steps:"){
            stream >> programOptions.steps;
        }
        else if (type == "#save_interval:"){
            stream >> programOptions.save_interval;
        }
        else if (type == "#timestep:"){
            stream >> programOptions.timestep;
        }
        else if (type == "#experiment:"){
            stream >> programOptions.experiment;
        }
        else if (type == "#lx:"){
            stream >> systemProps.lx;
        }
        else if (type == "#ly:"){
            stream >> systemProps.ly;
        }
        else if (type == "#lz:"){
            stream >> systemProps.lz;
        }
        else if (type == "#area_fraction:"){
            stream >> systemProps.area_fraction;
        }
        else if (type == "#base_height:"){
            stream >> systemProps.base_height;
        }
        else if (type == "#ball_height:"){
            stream >> systemProps.ball_height;
        }
        else if (type == "#dimple_spacing:"){
            stream >> systemProps.dimple_spacing;
        }
        else if (type == "#dimple_radius:"){
            stream >> systemProps.dimple_radius;
        }
        else if (type == "#dimple_depth:"){
            stream >> systemProps.dimple_depth;
        }
        else if (type == "#ball_radius:"){
            stream >> ballProps.radius;
        }
        else if (type == "#ball_mass:"){
            stream >> ballProps.mass;
        }
        else if (type == "#ball_youngs:"){
            stream >> ballProps.youngs_modulus;
        }
        else if (type == "#ball_poisson:"){
            stream >> ballProps.poisson;
        }
        else if (type == "#ball_damping:"){
            stream >> ballProps.damping_factor;
        }
        else if (type == "#base_radius:"){
            stream >> baseProps.radius;
        }
        else if (type == "#base_mass:"){
            stream >> baseProps.mass;
        }
        else if (type == "#base_youngs:"){
            stream >> baseProps.youngs_modulus;
        }
        else if (type == "#base_poisson:"){
            stream >> baseProps.poisson;
        }
        else if (type == "#base_damping:"){
            stream >> baseProps.damping_factor;
        }
        else {
            std::cout << "Unknown type: " << type << std::endl;
        }
        stream.ignore(100, '\n');
    }
    Options options;
    options.programOptions = programOptions;
    options.systemProps = systemProps;
    options.ballProps = ballProps;
    options.baseProps = baseProps;
    return options;
}