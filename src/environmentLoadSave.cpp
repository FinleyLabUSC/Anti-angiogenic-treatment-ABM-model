#include "Environment.h"

void Environment::loadParams() {
    /*
     * loads parameter files
     */
    std::ifstream dataCP(saveDir+"/params/cancerParams.csv");
    std::string line;
    while(std::getline(dataCP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        cancerParams.push_back(parsedRow[0]);
    }
    dataCP.close();

    std::ifstream dataCD(saveDir+"/params/cellDiameters.csv");
    while(std::getline(dataCD, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        cellDiameters.push_back(parsedRow[0]);
    }
    dataCD.close();

    std::ifstream dataEP(saveDir+"/params/envParams.csv");
    while(std::getline(dataEP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        envParams.push_back(parsedRow[0]);
    }
    dataEP.close();

    std::ifstream dataDP(saveDir+"/params/diffParams.csv");
    while(std::getline(dataDP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        diffParams.push_back(parsedRow[0]);
    }
    dataDP.close();

    std::ifstream dataDT(saveDir+"/params/treatParams.csv");
    while(std::getline(dataDT, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        treatParams.push_back(parsedRow[0]);
    }
    dataDT.close();
}

void Environment::save(Diffusibles &diff, double tstep) {
    /*
     * saves timecourse and spatial data
     */

    std::ofstream myfile;

    int numCancer = 0;
    for(auto & cell : cc_list){
        if(cell.state != "dead"){
            numCancer++;
        }
    }

    double time = steps*tstep/24;
    int numVessels = vessel_list.size();
    double td = tumorDiameter();

    myfile.open(saveDir+"/outputs.csv");
    myfile << time << "," << numCancer << "," << td << "," << chemoStartTime << "," << numVessels << std::endl;
    myfile.close();

    myfile.open(saveDir+"/cancerCells.csv");
    for(auto &cell : cc_list){
        int state = 0;
        int mState = 0;
        if(cell.state == "alive"){
            state = 0;
        } else if(cell.state == "necrotic"){
            state = 1;
        }
        if(cell.metabState == "proliferative"){
            mState = 0;
        } else if(cell.metabState == "quiescent"){
            mState = 1;
        } else if(cell.metabState == "non-proliferative"){
            mState = 2;
        }
        double thresh = cell.prolProb;//cell.prolThresh();
        myfile << cell.x[0] << "," << cell.x[1] << "," << cell.radius << "," << state << "," << mState << ","
               << static_cast<int>(cell.hypoxic) << "," << cell.pg << "," << cell.acidThreshold << ","
               << thresh << "," << cell.chemoDamage << "," << cell.chemoTolerance << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/vessels.csv");
    for(auto &ves : vessel_list){
        myfile << ves.x[0] << "," << ves.x[1] << "," << ves.radius << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/O2.csv");
    for(int i=0; i<diff.O2.size(); ++i){
        myfile << diff.O2[i][0];
        for(int j=1; j<diff.O2.size(); ++j){
            myfile << "," << diff.O2[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/glu.csv");
    for(int i=0; i<diff.glu.size(); ++i){
        myfile << diff.glu[i][0];
        for(int j=1; j<diff.glu.size(); ++j){
            myfile << "," << diff.glu[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/H.csv");
    for(int i=0; i<diff.H.size(); ++i){
        myfile << diff.H[i][0];
        for(int j=1; j<diff.H.size(); ++j){
            myfile << "," << diff.H[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/VEGF.csv");
    for(int i=0; i<diff.VEGF.size(); ++i){
        myfile << diff.VEGF[i][0];
        for(int j=1; j<diff.VEGF.size(); ++j){
            myfile << "," << diff.VEGF[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/chemo.csv");
    for(int i=0; i<diff.C.size(); ++i){
        myfile << diff.C[i][0];
        for(int j=1; j<diff.C.size(); ++j){
            myfile << "," << diff.C[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/cancerTS.csv");
    myfile << cancerTS[0];
    for(int i=1; i<cancerTS.size(); ++i){
        myfile << "," << cancerTS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/diameterTS.csv");
    myfile << diameterTS[0];
    for(int i=1; i<diameterTS.size(); ++i){
        myfile << "," << diameterTS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/o2TS.csv");
    myfile << o2TS[0];
    for(int i=1; i<o2TS.size(); ++i){
	myfile << "," << o2TS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/gluTS.csv");
    myfile << gluTS[0];
    for(int i=1; i<gluTS.size(); ++i){
	myfile << "," << gluTS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/phTS.csv");
    myfile << phTS[0];
    for(int i=1; i<phTS.size(); ++i){
	myfile << "," << phTS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/vesselTS.csv");
    myfile << vesselTS[0];
    for(int i=1; i<vesselTS.size(); ++i){
        myfile << "," << vesselTS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/hypoxicCancerTS.csv");
    myfile << hypoxicCancerTS[0];
    for(int i=1; i<hypoxicCancerTS.size(); ++i){
        myfile << "," << hypoxicCancerTS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/necroticCancerTS.csv");
    myfile << necroticCancerTS[0];
    for(int i=1; i<necroticCancerTS.size(); ++i){
        myfile << "," << necroticCancerTS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/prolCancerTS.csv");
    myfile << prolCancerTS[0];
    for(int i=1; i<prolCancerTS.size(); ++i){
        myfile << "," << prolCancerTS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/quiCancerTS.csv");
    myfile << quiCancerTS[0];
    for(int i=1; i<quiCancerTS.size(); ++i){
        myfile << "," << quiCancerTS[i];
    }
    myfile << std::endl;
    myfile.close();
}
