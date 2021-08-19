#ifndef CBMODEL_CANCERCELL_H
#define CBMODEL_CANCERCELL_H

#include <array>
#include <vector>
#include <random>
#include <cmath>
#include <random>
#include <string>
#include "Diffusibles.h"

class Cancer{
public:
    /*
     * FUNCTIONS
     */
    // initialize cell
    Cancer(std::array<double, 2> loc, double CD, double SD, int idx, std::vector<double> &cancerParams, bool initial, int N);

    // proliferation and passing on traits
    void prolifState();
    void inherit(double motherPG, double motherAcidThreshold, double motherProlProb, double motherChemoDamage, double motherChemoTolerance, double motherChemoAccumulated);
    void startCellCycle();
    double prolThresh();

    // interactions with environment
    void envEffects(double dt, Diffusibles &diff);
    void progressNecrotic(double dt);
    void migrate(double dt, Diffusibles &diff);
    void chemotherapy(double tstep, Diffusibles &diff);

    // force functions
    std::array<double, 2> attractiveForce(std::array<double, 2> dx, double otherRadius, double otherMu);
    std::array<double, 2> repulsiveForce(std::array<double, 2> dx, double otherRadius, double otherMu);
    void calculateForces(std::array<double, 2> otherX, double otherRadius, double otherMu, double otherProlProb, int otherType);
    void resolveForces(double dt, Diffusibles &diff);
    void resetForces();

    // overlap functions
    void calculateOverlap(std::array<double, 2> otherX, double otherRadius);
    void resetOverlap();
    void isCompressed();

    // cell behavior functions
    std::array<double, 2> proliferate();
    void age(double dt);

    // atp and nutrient functions
    void atpProduction(Diffusibles &diff);

    // other functions
    double calcDistance(std::array<double, 2> otherX);
    void generalLocation(double dx);
    void updateID(int idx);
    void neighboringCells(std::array<double, 2> otherX, int otherID, int otherType);

    /*
     * VARIABLES
     */
    // change in size when necrotic
    double necRadiusDecrease;
    double diameter;

    // impacts of environment
    double pgShift;
    double maxPG;
    double acidThresholdShift;
    double chemoAcidRes;
    double acRes;
    double acRate;
    double acThresh;

    // chemotherapy
    double chemoDamage;
    double chemoTolerance;
    double chemoTolRate;
    double chemoAccumulated;
    double chemoAccThresh;
    double chemoTime;
    double chemoTimeThresh;
    double chemoUptake;
    double chemoRepair;

    // proliferation and migration
    double prolProb;
    double prolProbShift;
    double migRate;
    double migBias;
    bool migratory;
    double prolProbBound;

    // location
    std::array<double, 2> x;
    std::array<int, 2> generalX;

    // physical properties
    double radius;
    bool compressed;
    double currentOverlap;
    std::vector<int> neighborsC;

    // age, division, and lifespan
    double totalAge;
    double divAge;
    double div;
    double lifespan;
    bool canProlif;
    int numDivs;

    // force properties
    double standardDiameter;
    double mu;
    double kc;
    double damping;
    double dampingFactor;
    double maxOverlap;
    double rmax;
    std::array<double, 2> currentForces;

    // atp and nutrients
    double Vo;
    double kO;
    double kG;
    double Ao;
    double pg;
    double kH;
    std::string metabState;
    bool hypoxic;
    double acidThreshold;
    double baseAcidThreshold;
    double minAcidThreshold;

    // identification
    int id;
    int type;
    std::string state;
    int maxN;

private:
    std::mt19937 mt;
};

#endif //CBMODEL_CANCERCELL_H
