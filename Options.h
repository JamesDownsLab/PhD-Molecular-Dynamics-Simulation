//
// Created by ppxjd3 on 21/07/2021.
//

#pragma once

#include <filesystem>
#include <string>
#include <fstream>
#include <iostream>

struct ProgramOptions {
    std::filesystem::path savepath{""};
    int steps{0};
    int save_interval{0};
    double timestep{0};
    std::string experiment{""};
};

struct SystemProps {
    double lx;
    double ly;
    double lz;
    double area_fraction;
    double base_height;
    double ball_height;
    double dimple_spacing;
    double dimple_radius;
    double dimple_depth;
};

struct ParticleProps {
    double radius;
    double mass;
    double youngs_modulus;
    double poisson;
    double damping_factor; // A in Poschel's book

};

struct Options {
    ProgramOptions programOptions;
    SystemProps systemProps;
    ParticleProps ballProps;
    ParticleProps baseProps;
};

Options read_input_file(const char* fname);
