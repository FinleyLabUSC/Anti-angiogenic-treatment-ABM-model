#include "CancerCell.h"
#include <cmath>
#include <iostream>

/*
 * REVISIONS
 * - necrotic cells
 *  - shrink to a certain size
 *  - no attractive force
 *  - reduce mu
 */

Cancer::Cancer(std::array<double, 2> loc, double CD, double SD, int idx, std::vector<double> &cancerParams, bool initial, int N): mt((std::random_device())()) {
    /*
     * x -> location of cell at inception
     * id -> index in cell storage vector
     * state -> current state. enters as "alive"
     * canProlif -> is able to proliferate. enters as true
     * standardDiameter -> used for calculating cell forces. set equal to cancer cell diameter
     * radius -> half of the cell diameter
     * diameter -> cell diameter in um
     * type -> cell type. 0 is cancer cell
     * maxN -> size of the PDE grid. keeps cells within simulation bounds
     */
    x = loc;
    id = idx;
    state = "alive";
    canProlif = true;
    standardDiameter = SD;
    radius = 0.5*CD;
    diameter = CD;
    type = 0;
    maxN = N-1;

    /*
     * FORCE PARAMETERS
     * same for all cells
     */
    mu = cancerParams[0];
    kc = cancerParams[1];
    damping = cancerParams[2];
    dampingFactor = cancerParams[23];
    rmax = 1.5*CD; // maximum interaction radius
    currentForces = {0,0}; // initialize current forces as 0

    /*
     * AGE AND DIVISION
     */
    totalAge = 0; // how old the cell is
    divAge = 0; // current spot in the cell cycle. cell divides at divAge=div and this resets
    // randomize cell cycle duration over a small range to prevent synchronous division
    std::uniform_real_distribution<double> dis(0.95*cancerParams[10],1.05*cancerParams[10]);
    div = dis(mt);
    // convert lifespan to hours. cell dies when totalAge=lifespan
    lifespan = cancerParams[11]*24;
    // if present at start of simulation, randomize cell age and cell cycle spot
    if(initial){
        std::uniform_real_distribution<double> divTime(0.0, div);
        std::uniform_real_distribution<double> lifeTime(0.0, lifespan);
        divAge = divTime(mt);
        totalAge = lifeTime(mt);
    }
    numDivs = cancerParams[12]; // maximum number of divisions

    /*
     * OVERLAP
     */
    maxOverlap = cancerParams[3]*CD; // fraction of cell diameter
    compressed = false; // assume not overlapped. doesn't really matter as this gets immediately updated

    /*
     * ATP AND NUTRIENTS
     */
    Vo = cancerParams[4];
    kO = cancerParams[5];
    kG = cancerParams[6];
    Ao = cancerParams[7];
    pg = cancerParams[8];
    kH = cancerParams[9];
    metabState = "proliferative"; // assume proliferative. updates immediately. if atp production is too low, no longer proliferative
    hypoxic = false; // assume not hypoxic. updates immediately
    necRadiusDecrease = cancerParams[14]; // rate of size decrease when necrotic
    acidThreshold = cancerParams[15]; // pH at which cell dies
    baseAcidThreshold =  cancerParams[15];

    /*
     * ENVIRONMENTAL EFFECTS
     * if hypoxic, increase glucose uptake -> upregulation of GLUT
     * if pH is close to acidThreshold, shift threshold lower
     */
    pgShift = cancerParams[16];
    maxPG = cancerParams[17];
    acidThresholdShift = cancerParams[18];
    minAcidThreshold = cancerParams[19];
    chemoAcidRes = cancerParams[25];
    acRes = 0.0;
    acRate = cancerParams[26];
    acThresh = cancerParams[27];

    /*
     * PROLIFERATION
     */
    prolProb = 1;
    prolProbShift = cancerParams[20];
    migratory = false;
    migRate = cancerParams[21];
    migBias = cancerParams[22];
    prolProbBound = cancerParams[24];

    /*
     * CHEMOTHERAPY
     */
    chemoDamage = 0;
    chemoTolerance = cancerParams[29];
    chemoTolRate = cancerParams[26];
    chemoAccumulated = 0;
    chemoAccThresh = cancerParams[28];
    chemoTime = 0;
    chemoTimeThresh = cancerParams[27];
    chemoUptake = cancerParams[30];
    chemoRepair = cancerParams[31];
}

void Cancer::prolifState() {
    /*
     * determine if the cell can proliferate
     */
    canProlif = !(state == "dead" || state == "necrotic" || metabState != "proliferative" || compressed || numDivs <= 0
            || migratory);
}

void Cancer::migrate(double dt, Diffusibles &diff, std::vector<Vessel> &vessels) {
    /*
     * if cell can migrate, randomly selects direction based on chemoattractant concentration
     * uses cell's general location on the PDE grid
     */
    if(state == "dead" || state == "necrotic" || !migratory){return;}

    for(auto &v : vessels){
        if(calcDistance(v.x) <= rmax + v.radius){
            return;
        }
    }

    std::array<double, 8> dx = {-1, -1, -1, 0, 0, 1, 1, 1};
    std::array<double, 8> dy = {-1, 0, 1, -1, 1, -1, 0, 1};

    std::vector<double> probs;
    std::vector<double> normProbs;
    double sum = 0;
    for (int i = 0; i < dx.size(); ++i) {
        probs.push_back(diff.O2[generalX[0] + dx[i]][generalX[1] + dy[i]]);
        normProbs.push_back(0);
        sum += probs[i];
    }

    int maxId = 0;
    double maxProb = probs[0];
    for(int i=1; i<probs.size(); ++i){
        if(probs[i] > maxProb){
            maxProb = probs[i];
            maxId = i;
        }
    }

    probs[maxId] *= migBias; // influence probability of choosing the direction with largest chemokine concentration
    sum = 0;
    for(int i=0; i<probs.size(); ++i){
        sum += probs[i];
    }

    for(int i=0; i<probs.size(); ++i){
        normProbs[i] = probs[i]/sum;
    }

    for(int i=1; i<probs.size(); ++i){
        normProbs[i] += normProbs[i-1];
    }

    std::uniform_real_distribution<double> dis(0.0,1.0);
    double p = dis(mt);

    int choice = 0;
    for(double norm_prob : normProbs){
        if(p > norm_prob){
            choice++;
        }
    }

    std::array<double, 2> targetX = {generalX[0] + dx[choice], generalX[1] + dy[choice]};
    // add a little randomness so that cell isn't moving at 90degree or 45degree angles
    double angle = atan2(targetX[1] - generalX[1], targetX[0] - generalX[0]);// + (dis(mt) - 0.5);
    double speed = migRate + static_cast<double>(hypoxic)*2*migRate; // Norton, JTB, 2018
    x[0] += dt * speed * cos(angle);
    x[1] += dt * speed * sin(angle);
}

void Cancer::progressNecrotic(double dt) {
    /*
     * if necrotic, shrink radius until removing cell
     */
    if(state != "necrotic"){return;}

    radius -= dt*necRadiusDecrease;
    /*if(radius <= 0.1*diameter){
        state = "dead";
    }*/
    radius = std::max(radius, 0.25*diameter);
}

void Cancer::envEffects(double dt, Diffusibles &diff) {
    /*
     * if pH is close to threshold, increase acid resistance
     * if hypoxic, upregulate GLUT to increase glucose uptake
     */

    if(hypoxic){
        pg += dt*pgShift;
        pg = std::min(pg, maxPG);
        prolProb -= dt*prolProbShift;
        prolProb = std::max(prolProb, prolProbBound);
    }

    double pH = -log10(diff.H[generalX[0]][generalX[1]]/(1e3));

    if(pH <= acidThreshold+0.1) {
        acidThreshold -= dt * acidThresholdShift;
        acidThreshold = std::max(acidThreshold, minAcidThreshold);
    }
}

void Cancer::chemotherapy(double tstep, Diffusibles &diff) {
    if(state == "dead" || state == "necrotic"){return;}

    double dt = 0.01;
    double steps = 60*60*tstep/dt;
    for(int i=0; i<steps; ++i) {
	    double drugUptake = chemoUptake*diff.C[generalX[0]][generalX[1]];
        chemoAccumulated += dt * drugUptake;
        chemoDamage += dt * (drugUptake - chemoRepair * chemoDamage);
        if (chemoAccumulated > chemoAccThresh) {
            chemoTime += dt;
        }
        if (chemoTime > chemoTimeThresh) {
            chemoTolerance += dt * chemoTolRate;
        }
    }
    if (chemoDamage > chemoTolerance) {
        state = "necrotic";
	mu *= 0.25;
    }
}

void Cancer::inherit(double motherPG, double motherAcidThreshold, double motherProlProb, double motherChemoDamage, double motherChemoTolerance, double motherChemoAccumulated) {
    /*
     * set glucose uptake and acidThreshold equal to the mother cell following division
     */
    pg = motherPG;
    acidThreshold = motherAcidThreshold;
    prolProb = motherProlProb;
    chemoDamage = motherChemoDamage;
    chemoTolerance = motherChemoTolerance;
    chemoAccumulated = motherChemoAccumulated;
}

void Cancer::startCellCycle() {
    if(state == "proliferative" || state == "dead" || state == "necrotic"){return;}

    std::uniform_real_distribution<double> dis(0.0,1.0);

    double thresh = prolProb;//prolThresh();

    if(dis(mt) < thresh){
        migratory = false;
        state = "proliferative";
    } else{
        migratory = true;
        state = "quiescent";
    }
}

double Cancer::prolThresh() {
    return prolProbBound + (1.0 - prolProbBound)/(1 + exp(-prolProb));
}

// force functions
std::array<double, 2> Cancer::attractiveForce(std::array<double, 2> dx, double otherRadius, double otherMu) {
    /*
     * attraction between neighboring, non-touching cells
     */
    double dxNorm = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
    std::array<double, 2> dxUnit = {dx[0]/dxNorm, dx[1]/dxNorm};
    double sij = radius/standardDiameter + otherRadius;

    double scaleFactor = otherMu*(dxNorm - sij)*exp(-kc*(dxNorm - sij)/sij);
    double F0 = dxUnit[0]*scaleFactor;
    double F1 = dxUnit[1]*scaleFactor;

    return {F0, F1};
}

std::array<double, 2> Cancer::repulsiveForce(std::array<double, 2> dx, double otherRadius, double otherMu) {
    /*
     * repulsion between overlapping cells
     */
    double dxNorm = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
    std::array<double, 2> dxUnit = {dx[0]/dxNorm, dx[1]/dxNorm};
    double sij = radius/standardDiameter + otherRadius;

    double scaleFactor = otherMu*sij*log10(1 + (dxNorm - sij)/sij);
    double F0 = dxUnit[0]*scaleFactor;
    double F1 = dxUnit[1]*scaleFactor;

    return {F0, F1};
}

void Cancer::calculateForces(std::array<double, 2> otherX, double otherRadius, double otherMu, double otherProlProb, int otherType) {
    /*
     * determine x and y components of forces between two cells
     * add to the current forces on the cells
     */
    // calculate actual distance between cell centers
    double distance = calcDistance(otherX);
    if(distance < rmax){
        // determine x and y distances between cell centers
        std::array<double, 2> dx = {(otherX[0]-x[0])/standardDiameter, (otherX[1]-x[1])/standardDiameter};
        if(distance < (radius + otherRadius)){
            std::array<double, 2> force = repulsiveForce(dx, otherRadius/standardDiameter, otherMu);
            currentForces[0] += force[0];
            currentForces[1] += force[1];
        } else if((type == 0 && otherType == 0) && state != "necrotic"){ // attraction if both are cancer cells
            double migratoryState = std::min(otherProlProb, prolThresh());
            std::array<double, 2> force = attractiveForce(dx, otherRadius/standardDiameter, otherMu);
            currentForces[0] += force[0]*migratoryState;
            currentForces[1] += force[1]*migratoryState;
        }
    }
}

void Cancer::resolveForces(double dt, Diffusibles &diff) {
    /*
     * numerically solve change in cell location based on forces
     */
    bool outer = (diff.normalCells[generalX[0]-1][generalX[1]-1]+diff.normalCells[generalX[0]][generalX[1]-1]+diff.normalCells[generalX[0]+1][generalX[1]-1]
                  +diff.normalCells[generalX[0]-1][generalX[1]]+diff.normalCells[generalX[0]][generalX[1]]+diff.normalCells[generalX[0]+1][generalX[1]]
                  +diff.normalCells[generalX[0]-1][generalX[1]+1]+diff.normalCells[generalX[0]][generalX[1]+1]+diff.normalCells[generalX[0]+1][generalX[1]+1]) >= 1;

    double d = damping;
    if(outer){d*=dampingFactor;}

    x[0] += (dt/d)*currentForces[0]*standardDiameter;
    x[1] += (dt/d)*currentForces[1]*standardDiameter;

    currentForces = {0,0};
}

void Cancer::resetForces() {
    /*
     * add random component to forces
     */
    std::random_device rd;
    double D = 1;
    std::uniform_real_distribution<double> dis(-D, D);
    currentForces = {dis(rd),dis(rd)};
}

void Cancer::neighboringCells(std::array<double, 2> otherX, int otherID, int otherType){
    /*
     * determine which cells are neighboring the current cell
     * used to speed up calculation of cell forces
     */
    double dis = calcDistance(otherX);
    if(dis <= 2*rmax){
        if(otherType == 0) {
            neighborsC.push_back(otherID);
        }
    }
}

// overlap functions
void Cancer::calculateOverlap(std::array<double, 2> otherX, double otherRadius, std::string otherType) {
    /*
     * determine if  two  cells are overlapping
     */
    double distance = calcDistance(otherX);
    if(distance < radius + otherRadius && otherType != "necrotic"){
        currentOverlap += radius + otherRadius - distance;
    }
}

void Cancer::resetOverlap() {
    // reset overlap
    currentOverlap = 0;
}

void Cancer::isCompressed() {
    /*
     * determine if a cell is compressed
     */
    compressed = currentOverlap > maxOverlap;

}

// cell behavior functions
std::array<double, 2> Cancer::proliferate() {
    /*
     * if able to proliferate, do so in direction away from current forces on the cell
     * place overlapping the mother cell and then the cells repel each other
     */
    if(!canProlif){return {0,0};}

    // cost associated with increased acid resistance -> minAcidThreshold = double cell cycle time
    double atScale = (baseAcidThreshold - acidThreshold)/(baseAcidThreshold - minAcidThreshold);

    if(divAge >= (div + atScale*div)){
        double angle = 0;
        if(currentForces[0]!=0 && currentForces[1]!=0){
            angle = atan2(currentForces[1],currentForces[0]);
        } else{
            std::random_device rd;
            std::uniform_real_distribution<double> dis(0.0, 2*3.1415);
            angle = dis(rd);
        }
        divAge = 0; // reset division  age
        numDivs -= 1; // subtract from max  number of divisions
        state = "quiescent";
        chemoAccumulated /= 2.0;
        return{0.5*radius*cos(angle)+x[0], 0.5*radius*sin(angle)+x[1]};
    } else{
        return {0,0};
    }
}

void Cancer::age(double dt) {
    /*
     * increase cell age
     */
    if(state == "necrotic"){return;}

    if(state == "proliferative") {
        divAge += dt;
    }

    totalAge += dt;
    if(totalAge >= lifespan){
        state = "dead";
    }
}

// atp and nutrient functions
void Cancer::atpProduction(Diffusibles &diff) {
    /*
     * calculate ATP production rate based on nutrient concentrations
     * then determine cell state
     */
    if(state == "necrotic"){return;}

    double wO = -Vo * diff.O2[generalX[0]][generalX[1]] / (kO + diff.O2[generalX[0]][generalX[1]]);
    double wG = -(pg * Ao / 2 + 27 * wO / 10) * diff.glu[generalX[0]][generalX[1]] / (kG + diff.glu[generalX[0]][generalX[1]]);
    double wA = -(2*wG + 27*wO/5);
    double A = wA/Ao; // A is ratio between ATP production rate and ATP requirement

    if(A >= 0.8){
        metabState = "proliferative";
    } else if(A >=0.3){
        metabState = "quiescent";
    } else{
        if(type == 0){
            // cancer cells become necrotic and shrink before dying
            state = "necrotic";
            mu *= 0.25;
            metabState = "non-proliferative";
        } else {
            // other cells just die
            metabState = "non-proliferative";
            state = "dead";
        }
    }

    // if pH is too low, cell dies
    if(-log10(diff.H[generalX[0]][generalX[1]]*1e-3) < acidThreshold){
        if(type == 0){
            state = "necrotic";
            mu *= 0.25;
            metabState = "non-proliferative";
        } else {
            metabState = "non-proliferative";
            state = "dead";
        }
    }

    // if O2 concentration is too low, cell becomes hypoxic
    hypoxic = (diff.O2[generalX[0]][generalX[1]] < 0.002) && state!="dead" && state!="necrotic";
}

// other functions
double Cancer::calcDistance(std::array<double, 2> otherX) {
    /*
     * distance between cell centers
     */
    double d0 = (otherX[0] - x[0]);
    double d1 = (otherX[1] - x[1]);
    return sqrt(d0*d0 + d1*d1);
}

void Cancer::generalLocation(double dx) {
    /*
     * discretize location to points on the PDE grid
     */
    double i = x[0] + dx/2;
    i -= fmod(i, dx);
    i /= dx;
    double j = x[1] + dx/2;
    j -= fmod(j, dx);
    j /= dx;

    generalX = {static_cast<int>(i), static_cast<int>(j)};
    // if out of bounds, kill cell
    if(generalX[0]<=1 || generalX[1]<=1 || generalX[0]>=maxN-1 || generalX[1]>=maxN-1){
        state = "dead";
    }
}

void Cancer::updateID(int idx) {
    /*
     * update index in cell vector
     */
    id = idx;
}
