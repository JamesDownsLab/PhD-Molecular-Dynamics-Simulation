//
// Created by ppxjd3 on 14/07/2021.
//

#include "Engine.h"

Engine::Engine(const char *fname, ProgramOptions options) : _options{options} {
    f1 = fopen(_options.savepath.string().c_str(), "w");
    init_system(fname);

}

void Engine::init_system(const char *fname) {
    std::cout << fname << std::endl;
}
