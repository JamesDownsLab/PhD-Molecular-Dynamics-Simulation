#include <iostream>
#include "Engine.h"
#include "setup.h"

int main() {
    ProgramOptions options{
        "data.dump"
    };
    setup_experiment(0.8);
    Engine engine("initial.data", options);
}
