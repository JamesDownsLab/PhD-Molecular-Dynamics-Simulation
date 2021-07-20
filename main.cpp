#include <iostream>
#include "Engine.h"
#include <sstream>
#include <string>

int main(int argc, char** argv){
    int steps{1000};
    int save_interval{100};
    for (int i = 0; i<argc; i++){
        std::cout << argv[i] << "\n";
        std::string command = argv[i];
        if (command == "--steps"){
            steps = atoi(argv[i+1]);
            i++;
        }
        if (command == "--save_interval"){
            save_interval = atoi(argv[i+1]);
            i++;
        }
    }

    ProgramOptions options{
            "data.dump",
            save_interval
    };
//    setup_experiment(0.8);
    Engine engine("initial.data", options);
    engine.set_baseplate(2e-4, 0.02);
    for (int s{0}; s <= steps; s++) {
        engine.step();
    }
}