#include <iostream>
#include "Engine.h"
#include <string>
#include "Options.h"


int main(int argc, char** argv){
    const char* fname;
    for (int i = 0; i<argc; i++){
        std::cout << argv[i] << "\n";
        std::string command = argv[i];
        if (command == "--in"){
            fname = argv[i+1];
        }
    }

    Options options = read_input_file(fname);
    Engine engine(options);


    if (options.programOptions.experiment == "stable") {
        engine.set_baseplate(options.programOptions.amplitude, 0.02);
        for (int s{0}; s <= options.programOptions.steps; s++) {
            engine.step();
        }
    }

    else if (options.programOptions.experiment == "startstop"){
        engine.set_baseplate(options.programOptions.amplitude, 0.02);
        for (int s{0}; s <= options.programOptions.steps/2; s++) {
            engine.step();
        }
        engine.set_baseplate(0.0, 0.02);
        for (int s{0}; s<=options.programOptions.steps/2; s++){
            engine.step();
        }
    }

    else {
        std::cout << "Experiment not specified" << std::endl;
    }
}