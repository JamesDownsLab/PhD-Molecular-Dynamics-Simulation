#include <iostream>
#include "Engine.h"


int main() {
    ProgramOptions options{
        "data.dump"
    };
    Engine engine("initial.data", options);

}
