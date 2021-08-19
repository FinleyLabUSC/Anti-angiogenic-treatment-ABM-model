#include "Environment.h"

void Environment::calculateForces(double tstep, Diffusibles &diff) {
    /*
     * calculate forces between cells
     * first, determine neighboring cells
     * then for steps of duration dt in tstep, calculate and resolve forces between neighbors
     */

    int Nsteps = static_cast<int>(tstep/dt);

    // determine cell neighbors
#pragma omp parallel for
    for(int i=0; i<cc_list.size(); ++i){
        cc_list[i].neighborsC.clear();
        cc_list[i].resetForces();
        for(auto &cell : cc_list){
            if(cc_list[i].id != cell.id) {
                cc_list[i].neighboringCells(cell.x,cell.id, 0);
            }
        }
    }

    // iterate through steps: calculate and resolve forces
    for(int i=0; i<Nsteps; ++i){
#pragma omp parallel for
        for(int i=0; i<cc_list.size(); ++i){
            cc_list[i].resetForces();
            for(auto &j : cc_list[i].neighborsC){
                cc_list[i].calculateForces(cc_list[j].x, cc_list[j].radius, cc_list[j].mu, cc_list[j].prolThresh(), cc_list[j].type);
            }
            for(auto& ves : vessel_list){
                cc_list[i].calculateForces(ves.x, ves.radius, ves.mu, 1.0, 3);
            }
        }

#pragma omp parallel for
        for(int i=0; i<cc_list.size(); ++i){
            cc_list[i].resolveForces(dt, diff);
        }
    }

    // calculate overlap between cells
#pragma omp parallel for
    for(int i=0; i<cc_list.size(); ++i){
        cc_list[i].resetOverlap();
        for(auto &cell2 : cc_list){
            if(cc_list[i].id != cell2.id) {
                cc_list[i].calculateOverlap(cell2.x, cell2.radius);
            }
        }
        cc_list[i].isCompressed();
        cc_list[i].prolifState();
        cc_list[i].generalLocation(dx);
    }
    // see if vessels collapse from cancer cell pressure
#pragma omp parallel for
    for(int i=0; i<vessel_list.size(); ++i){
        vessel_list[i].resetOverlap();
        vessel_list[i].age += tstep;
        for(auto &cell2 : cc_list){
            vessel_list[i].calculateOverlap(cell2.x, cell2.radius);
        }
        if(vessel_list[i].age  >= vesAge) {
            vessel_list[i].isCompressed();
        }
        vessel_list[i].treatment(angioTreatmentOn);
    }

    std::vector<Vessel> new_vessels;
    for(auto &vessel : vessel_list){
        if(vessel.state != "dead"){
            new_vessels.push_back(vessel);
        }
    }
    vessel_list = new_vessels;
}

void Environment::runCancer(double tstep, Diffusibles &diff) {
    /*
     * cancer cell functions
     */

#pragma omp parallel for
    for(int i=0; i<cc_list.size(); ++i){
        cc_list[i].generalLocation(dx);
        cc_list[i].migrate(tstep, diff);
        cc_list[i].resetForces();
        for(auto &cell2 : cc_list){
            if(cc_list[i].id != cell2.id){
                cc_list[i].calculateForces(cell2.x, cell2.radius, cell2.mu, cell2.prolThresh(), cell2.type);
            }
        }
        for(auto& ves : vessel_list){
            cc_list[i].calculateForces(ves.x, ves.radius, ves.mu, 1.0, 3);
        }
        cc_list[i].isCompressed();

        cc_list[i].generalLocation(dx);
        cc_list[i].atpProduction(diff);
        cc_list[i].envEffects(tstep, diff);
        cc_list[i].chemotherapy(tstep, diff);
        cc_list[i].progressNecrotic(tstep);
        cc_list[i].startCellCycle();
        cc_list[i].prolifState();
        cc_list[i].age(tstep);
    }

    int numCells = cc_list.size();
    for(int i=0; i<numCells; ++i){
        std::array<double, 2> newLoc = cc_list[i].proliferate();
        if(newLoc[0]!=0 && newLoc[1]!=0){
            cc_list.push_back(Cancer(newLoc, cancerDiameter, standardDiameter, cc_list.size(), cancerParams, false,
                                     static_cast<int>(envSize/dx)));
            cc_list[cc_list.size() - 1].inherit(cc_list[i].pg,
                                                cc_list[i].acidThreshold,
                                                cc_list[i].prolProb,
                                                cc_list[i].chemoDamage,
                                                cc_list[i].chemoTolerance,
                                                cc_list[i].chemoAccumulated);
        }
    }

    std::vector<Cancer> new_cc;
    for(auto &cell : cc_list){
        if(cell.state != "dead"){new_cc.push_back(cell);}
    }
    cc_list = new_cc;
    for(int i=0; i<cc_list.size(); ++i){
        cc_list[i].updateID(i);
    }
}