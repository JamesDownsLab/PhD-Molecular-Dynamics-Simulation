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

    else if (options.programOptions.experiment == "ramp"){
        std::cout << "Stating ramp experiment" << std::endl;
        engine.set_baseplate(options.programOptions.amplitude_start, 0.02);
        double start = options.programOptions.amplitude_start;
        double amp = start;
        double end = options.programOptions.amplitude_end;
        double rate = options.programOptions.ramp_rate;
        double time = std::abs((start - end)/rate);
        double steps = round(time / options.programOptions.timestep);
        double interval = (end - start)/steps;
        for (int s{0}; s<=steps; s++){
            amp += interval;
            engine.set_baseplate(amp, 0.02);
            engine.step();
        }
    }

    else {
        std::cout << "Experiment not specified" << std::endl;
    }
}