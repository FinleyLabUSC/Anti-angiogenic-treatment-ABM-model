#ifndef CBMODEL_VESSEL_H
#define CBMODEL_VESSEL_H

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <string>

class Vessel{
public:
    Vessel(std::array<double, 2> loc, double CD, double mol, double hm, double hd, int recruited);
    void treatment(bool treatmentOn);

    // overlap functions
    void calculateOverlap(std::array<double, 2> otherX, double otherRadius);
    void resetOverlap();
    void isCompressed();
    double calcDistance(std::array<double, 2> otherX);
    void generalLocation(double dx);

    double radius;
    std::array<double, 2> x;
    std::array<int, 2> generalX;
    double maxOverlap;
    double currentOverlap;
    std::string state;
    double age;
    double mu;

    double healthyMu;
    double health;
};

#endif //CBMODEL_VESSEL_H
