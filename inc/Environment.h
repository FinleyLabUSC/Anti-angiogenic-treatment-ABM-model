#ifndef CBMODEL_ENVIRONMENT_H
#define CBMODEL_ENVIRONMENT_H

#include <vector>
#include "CancerCell.h"
#include "Diffusibles.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>
#include "Vessel.h"

/*
 * Main simulation object. Uses vectors to store cells. Loads all the parameters.
 * Runs cell functions.
 */

class Environment{
public:
    /*
     * FUNCTIONS that are called by main.cpp
     */
    Environment(std::string folder, std::string set, double size);
    void simulate(double tstep);

private:
    /*
     * FUNCTIONS
     */
    // RunCells
    void runCancer(double tstep, Diffusibles &diff);
    void calculateForces(double tstep, Diffusibles &diff);

    // Misc
    void diffusion(double tstep, Diffusibles &diff);
    void initializeCells();
    void angiogenesis(Diffusibles &diff);
    void updateCellGrids(Diffusibles &diff);
    double tumorDiameter();
    void treatmentOnOff(double tstep);

    // LoadSave
    void save(Diffusibles &diff, double tstep);
    void loadParams();

    // Info
    void printStep(Diffusibles &diff, double time);
    void updateTimeSeries(Diffusibles &diff);

    /*
     * VARIABLES
     */
    // cell lists
    std::vector<Cancer> cc_list;
    std::vector<Vessel> vessel_list;

    // time courses
    std::vector<int> cancerTS;
    std::vector<int> vesselTS;
    std::vector<int> hypoxicCancerTS;
    std::vector<int> necroticCancerTS;
    std::vector<int> prolCancerTS;
    std::vector<int> quiCancerTS;
    std::vector<double> phTS;
    std::vector<double> o2TS;
    std::vector<double> gluTS;
    std::vector<double> diameterTS;
    std::vector<double> vegfTS;
    std::vector<int> aliveCancerTS;

    // parameter lists
    std::vector<double> cancerParams;
    std::vector<double> cellDiameters;
    std::vector<double> envParams;
    std::vector<double> diffParams;
    std::vector<double> treatParams;

    // cell diameters
    double standardDiameter;
    double cancerDiameter;

    // treatment
    bool angioTreatmentOn;
    bool chemotherapyOn;
    double angioTreatment;
    double startDiameterAngio;
    double startDiameterChemo;
    double endDiameter;
    double tumorD;
    double chemoStartTime;
    double angioOnDuration;
    double angioOnTime;
    double angioOffDuration;
    double angioOffTime;
    double chemoOnDuration;
    double chemoOnTime;
    double chemoOffDuration;
    double chemoOffTime;
    bool angioTreatmentInit;
    bool chemotherapyInit;

    // other parameters
    std::string saveDir;
    int steps;
    double dx;
    double dt;
    double envSize;
    double vesRec;
    double vesAge;
    double vesMu;
    double vesDec;
    double vesDens;

    // visualization
    float renderScale;
    float windowSize;

    std::mt19937 mt;
};

#endif //CBMODEL_ENVIRONMENT_H
