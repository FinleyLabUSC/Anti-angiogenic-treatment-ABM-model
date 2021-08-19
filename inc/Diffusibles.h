#include <vector>
#include <omp.h>
#include <cmath>
#include <iostream>

#ifndef CBMODEL_DIFFUSIBLES_H
#define CBMODEL_DIFFUSIBLES_H

#include "Vessel.h"

/*
 * holds grids of diffusible factors along with diffusion functions
 */

class Diffusibles{
public:
    /*
     * FUNCTIONS
     */
    // initialize diffusibles
    Diffusibles(double size, double spatialStep, std::vector<double> &diffParams);
    void setBounds();

    // set current state and run diffusion
    void runDiffusion(double time, bool chemoOn);
    void vesselLocations(std::vector<Vessel> vessel_list);
    void setMaxPG(double cellPG);

    /*
     * GRIDS
     */
    // metabolism related
    std::vector<std::vector<double>> O2;
    std::vector<std::vector<double>> glu;
    std::vector<std::vector<double>> H;

    // cytokines
    std::vector<std::vector<double>> VEGF;

    // chemotherapy
    std::vector<std::vector<double>> C;

    // cell locations (not using some of these)
    std::vector<std::vector<double>> cancerCells;
    std::vector<std::vector<double>> vessels;
    std::vector<std::vector<double>> hypoxicCells;
    std::vector<std::vector<double>> normalCells;

    // parameter values at each location
    std::vector<std::vector<double>> Vo;
    std::vector<std::vector<double>> kO;
    std::vector<std::vector<double>> kG;
    std::vector<std::vector<double>> Ao;
    std::vector<std::vector<double>> pg;
    std::vector<std::vector<double>> kH;
    std::vector<std::vector<double>> dc;

    // spatial step
    double dx;

    // diffusion coefficients
    double Do;
    double Dg;
    double Dh;
    double DVEGF;
    double Dc;

    // cytokine secretion rates
    double kVEGF;

    // chemotherapy uptake and decay rate
    double kupC;
    double kdegC;

    // metabolism of normal tissue
    double baseVo;
    double basekO;
    double basekG;
    double baseAo;
    double basepg;
    double basekH;
    double maxPG;

private:
    // numerical solvers
    void finDiff(std::vector<std::vector<double>> &g, std::vector<std::vector<double>> &g_next, int i, int j, double dt, double D);

    // number of grid points
    int N;

    // diffusible factor values at time = t+1
    std::vector<std::vector<double>> O2_next;
    std::vector<std::vector<double>> glu_next;
    std::vector<std::vector<double>> H_next;
    std::vector<std::vector<double>> VEGF_next;
    std::vector<std::vector<double>> C_next;


    // maximum % difference in diffusible factor concentration between time t+1 and t
    std::vector<std::vector<double>> c;
};

#endif //CBMODEL_DIFFUSIBLES_H
