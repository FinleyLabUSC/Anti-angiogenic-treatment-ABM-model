#include "Vessel.h"
#include <iostream>

Vessel::Vessel(std::array<double, 2> loc, double CD, double mol, double hm, double hd, int recruited) {
    /*
     * treat like simple cells that are fixed in place
     */
    x = loc;
    radius = 0.5*CD;
    maxOverlap = mol*CD;
    state = "alive";
    currentOverlap = 0;
    age = 0;

    healthyMu = hm;

    if(recruited == 1){
        health = hd;
    } else{
        health = 1.0;
    }

    mu = healthyMu*health;
}

void Vessel::calculateOverlap(std::array<double, 2> otherX, double otherRadius) {
    double distance = calcDistance(otherX);
    if(distance < radius + otherRadius){
        currentOverlap += radius + otherRadius - distance;
    }
}

void Vessel::resetOverlap() {
    currentOverlap = 0;
}

void Vessel::isCompressed() {
    /*
     * remove if compressed by cancer cell
     */
    if(currentOverlap > maxOverlap){
        state = "dead";
    }
}

double Vessel::calcDistance(std::array<double, 2> otherX) {
    double d0 = (otherX[0] - x[0]);
    double d1 = (otherX[1] - x[1]);
    return sqrt(d0*d0 + d1*d1);
}

void Vessel::generalLocation(double dx) {
    double i = x[0] + dx/2;
    i -= fmod(i, dx);
    i /= dx;
    double j = x[1] + dx/2;
    j -= fmod(j, dx);
    j /= dx;

    generalX = {static_cast<int>(i), static_cast<int>(j)};
}

void Vessel::treatment(bool treatmentOn) {
    if(treatmentOn){
        mu = healthyMu;
        health = 1.0;
    }
}