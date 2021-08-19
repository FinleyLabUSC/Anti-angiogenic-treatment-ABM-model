#include "Environment.h"

void Environment::diffusion(double tstep, Diffusibles &diff) {
    /*
     * set cell locations on PDE grids
     * then diffuse
     */
#pragma omp parallel for default(none)
    for(int i=0; i<cc_list.size(); ++i){
        cc_list[i].generalLocation(dx);
    }

    updateCellGrids(diff);
    diff.vesselLocations(vessel_list);
    diff.runDiffusion(tstep, chemotherapyOn);
}

void Environment::initializeCells() {
    /*
     * randomize initial cells and vasculature
     * only starting with cancer cells right now
     */
    int N = static_cast<int>(envSize/dx);

    // cancer cells
    int z = 0;
    for(int i=(N/2)-2; i<(N/2)+3; ++i){
        for(int j=(N/2)-2; j<(N/2)+3; ++j){
            cc_list.push_back(Cancer({i*dx, j*dx}, cancerDiameter, standardDiameter, z, cancerParams, true, static_cast<int>(envSize/dx)));
            z++;
        }
    }

    // vessels
    // from Robertson-Tessi 2015
    int Nves = envSize*envSize/(vesDens*vesDens);
    std::uniform_real_distribution<double> spacing(vesDens-50,vesDens+50);
    std::uniform_real_distribution<double> vesDis(50.0,envSize-50.0);
    for(int q=0; q<Nves; ++q){
        double i = vesDis(mt);
        double j = vesDis(mt);
        double space = spacing(mt);
        bool badSpace = true;
        int tries = 0;
        while(badSpace){
            tries++;
            if(tries > 10000){
                // if can't find a suitable location, just choose a random one
                badSpace = false;
            }
            double minDis = 250;
            for(auto & vessel : vessel_list){
                double di = i - vessel.x[0];
                double dj = j - vessel.x[1];
                double distance = sqrt(di*di + dj*dj);
                minDis = std::min(distance, minDis);
            }
            if(minDis < space){
                i = vesDis(mt);
                j = vesDis(mt);
            } else{
                badSpace = false;
            }
        }
        vessel_list.push_back(Vessel({i, j}, cellDiameters[1], envParams[2], vesMu, vesDec, 0));
        vessel_list[vessel_list.size()-1].generalLocation(dx);
    }
}

void Environment::updateCellGrids(Diffusibles &diff) {
    /*
     * set  locations of cells on PDE grids and set secretion rates
     */
    for(int i=0; i<diff.O2.size(); ++i){
        for(int j=0; j<diff.O2.size(); ++j){
            diff.cancerCells[i][j]  = 0;
            // set normal tissue uptake
            diff.Vo[i][j] = diff.baseVo*diff.normalCells[i][j];
            diff.kO[i][j] = diff.basekO*diff.normalCells[i][j];
            diff.kG[i][j] = diff.basekG*diff.normalCells[i][j];
            diff.pg[i][j] = diff.basepg*diff.normalCells[i][j];
            diff.kH[i][j] = diff.basekH*diff.normalCells[i][j];
            diff.Ao[i][j] = diff.baseAo*diff.normalCells[i][j];
            diff.hypoxicCells[i][j] = 0;
        }
    }

    for(auto &cell : cc_list){
        if(cell.state == "necrotic"){continue;}

        diff.cancerCells[cell.generalX[0]][cell.generalX[1]] = 1;
        diff.normalCells[cell.generalX[0]][cell.generalX[1]] = 0; // permanently remove normal tissue where there's a cancer cell
        diff.Vo[cell.generalX[0]][cell.generalX[1]] = cell.Vo;
        diff.kO[cell.generalX[0]][cell.generalX[1]] = cell.kO;
        diff.Ao[cell.generalX[0]][cell.generalX[1]] = cell.Ao;
        diff.kG[cell.generalX[0]][cell.generalX[1]] = cell.kG;
        diff.pg[cell.generalX[0]][cell.generalX[1]] = cell.pg;
        diff.kH[cell.generalX[0]][cell.generalX[1]] = cell.kH;
        if(cell.hypoxic){
            diff.hypoxicCells[cell.generalX[0]][cell.generalX[1]] = 1;
        }
        diff.setMaxPG(cell.pg);
    }

}

void Environment::angiogenesis(Diffusibles &diff) {
    /*
     * recruit vessels randomly to a spot based on VEGF concentration
     */
    std::uniform_real_distribution<double> dis(0.0,1.0);

    for(int i=5; i<diff.VEGF.size()-5; ++i){
        for(int j=5; j<diff.VEGF.size()-5; ++j){
            if(dis(mt) < diff.VEGF[i][j]*vesRec && (diff.normalCells[i][j] + diff.cancerCells[i][j] >= 1)){
                vessel_list.push_back((Vessel({i*dx,j*dx}, cellDiameters[1], envParams[2], vesMu, vesDec, 1)));
                vessel_list[vessel_list.size()-1].generalLocation(dx);
            }
        }
    }
}

double Environment::tumorDiameter() {
    /*
     * calculates the rough diameter of the tumor as the distance between the farthest appart tumor cells
     */

    double maxDist = 0;
    for(auto& cell : cc_list){
        for (auto &cell2 : cc_list) {
            if (cell.id != cell2.id) {
                double dx0 = cell.x[0] - cell2.x[0];
                double dx1 = cell.x[1] - cell2.x[1];
                double dist = sqrt(dx0 * dx0 + dx1 * dx1);
                maxDist = std::max(maxDist, dist);
            }
        }
    }

    return maxDist;
}

void Environment::treatmentOnOff(double tstep) {
    if(angioTreatmentInit){
        // if treatment has been started, cycle treatment
        if(!angioTreatmentOn){
            // if treatment is off, progress the "off" clock
            angioOffTime += tstep;
            if(angioOffTime >= angioOffDuration){
                // if "off" clock reaches its duration, turn treatment on
                angioTreatmentOn = true;
                // reset "off" clock
                angioOffTime = 0.0;
                // adjust vessel recruitment rate to treatment value
                vesRec /= angioTreatment;
            }
        } else{
            // if treatment is on, progress "on" clock
            angioOnTime += tstep;
            if(angioOnTime >= angioOnDuration){
                // if "on" clock reaches its duration, turn treatment off
                angioTreatmentOn = false;
                // reset "on" clock
                angioOnTime = 0.0;
                // reset vessel recruitment to base value
                vesRec = envParams[1];
            }
        }
    }

    if(chemotherapyInit){
        // if treatment has been started, cycle treatment
        if(!chemotherapyOn){
            // if treatment is off, progress the "off" clock
            chemoOffTime += tstep;
            if(chemoOffTime >= chemoOffDuration){
                // if "off" clock reaches its duration, turn treatment on
                chemotherapyOn = true;
                // reset "off" clock
                chemoOffTime = 0.0;
            }
        } else{
            // if treatment is on, progress "on" clock
            chemoOnTime += tstep;
            if(chemoOnTime >= chemoOnDuration){
                // if "on" clock reaches its duration, turn treatment off
                chemotherapyOn = false;
                // reset "on" clock
                chemoOnTime = 0.0;
            }
        }
    }
}