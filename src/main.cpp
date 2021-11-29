#include <iostream>
#include "Environment.h"

/*
 * First: shifting properties of pg and acidThresh, and damping
 *
 * Second: at chosen properties,  vary prolProbShift
 *
 * Start angiogenesis once tumor reaches X size -> different strengths of treatment
 * Use rate of prolProbShift to model different metastatic potentials
 */

int main(int argc, char **argv) {
    std::cout << "Reading arguments\n";
    std::string folder = argv[1];
    std::string set = argv[2];

    std::cout << "Making save folder\n";
    std::string str = "mkdir -p ./"+folder+"/set_" + set + "/timeStates";
    const char *command = str.c_str();
    std::system(command);

    std::cout << "Generating parameters\n";
    str = "python3 genParams.py ./"+folder+"/set_"+set+" "+set;
    command = str.c_str();
    std::system(command);

    std::cout << "Starting simulation\n";
    double start = omp_get_wtime();
    Environment model(folder, set, 2000);
    model.simulate(0.05);
    double stop = omp_get_wtime();
    std::cout << "Duration: " << (stop-start)/(60*60) << std::endl;

    return 0;
}
