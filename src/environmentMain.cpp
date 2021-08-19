#include "Environment.h"

Environment::Environment(std::string folder, std::string set, double size): mt((std::random_device())()) {
    /*
     * creates a simulation
     * set save directory and load parameters
     */
    saveDir = "./"+folder+"/set_"+set;
    loadParams();

    cancerDiameter = cellDiameters[0];
    standardDiameter = cancerDiameter;

    angioTreatment = treatParams[0];
    angioTreatmentOn = false;
    chemotherapyOn = false;
    vesRec = envParams[1];
    vesAge = envParams[3];
    startDiameterAngio = treatParams[1];
    startDiameterChemo = treatParams[4];
    endDiameter = treatParams[7];
    tumorD = 0;
    vesMu = envParams[4];
    vesDec = envParams[5];
    vesDens = envParams[0];

    angioOnDuration = treatParams[2];
    angioOnTime = 0;
    angioOffDuration = treatParams[3];
    angioOffTime = 0;
    chemoOnDuration = treatParams[5];
    chemoOnTime = 0;
    chemoOffDuration = treatParams[6];
    chemoOffTime = 0;

    angioTreatmentInit = false;
    chemotherapyInit = false;

    chemoStartTime = 0;

    // set PDE spatial step to slightly smaller than cancer cell diameter
    dx = cancerDiameter*0.75;
    envSize = size;
    steps = 0;
    // windowSize and renderScale are for when visualization is added
    windowSize = 2000.f;
    renderScale = windowSize/size;

    // time step for calculating forces
    dt = 0.005;
}

void Environment::simulate(double tstep) {
    // initialize environment and set diffusion bounds
    Diffusibles diff(envSize, dx, diffParams);
    initializeCells();
    for(auto & cell : cc_list){
        cell.generalLocation(dx);
    }
    updateCellGrids(diff);
    diff.vesselLocations(vessel_list);
    diff.setBounds();

    while(true){
        diffusion(tstep, diff);
        calculateForces(tstep, diff);
        runCancer(tstep, diff);
        angiogenesis(diff);

        steps += 1;
        //printStep(diff, steps*tstep);
        if(fmod(steps*tstep, 1) == 0) {
            updateTimeSeries(diff);
        }
        if(fmod(steps*tstep, 24) == 0) {
            save(diff,tstep);
            treatmentOnOff(tstep);
            if(tumorD >= startDiameterAngio && !angioTreatmentInit){
                // first time reaching start diameter,
                // turn on treatment
                angioTreatmentInit = true;
                angioTreatmentOn = true;
                vesRec /= angioTreatment;
            }
            if(tumorD >= startDiameterChemo && !chemotherapyInit){
                // first time reaching start diameter,
                // turn on treatment
                chemotherapyInit = true;
                chemotherapyOn = true;
                chemoStartTime = static_cast<double>(steps)*tstep/24;
            }
            if(tumorD >= endDiameter){break;}
	    if(steps*tstep/24 >= 100){break;}
	    //if(steps*tstep/24 >= endDiameter){break;}
        }

        if(cc_list.empty()){break;}
    }
    save(diff, tstep);
}
