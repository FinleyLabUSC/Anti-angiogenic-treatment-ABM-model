#include "Environment.h"

void Environment::printStep(Diffusibles &diff, double time) {
    /*
     * prints current state of the model
     */
    double maxH = 0;
    double maxVEGF = 0;
    double avgO2 = 0;
    double avgGlu = 0;
    double maxO2 = 0;
    double maxGlu = 0;
    double minH = 1e8;
    for(int i=0; i<diff.O2.size(); ++i){
        for(int j=0; j<diff.O2.size(); ++j){
            maxVEGF = std::max(maxVEGF, diff.VEGF[i][j]);
            maxH = std::max(maxH, diff.H[i][j]);
            avgO2 += diff.O2[i][j];
            avgGlu += diff.glu[i][j];
            maxO2 = std::max(maxO2, diff.O2[i][j]);
            maxGlu = std::max(maxGlu, diff.glu[i][j]);
            minH = std::min(minH, diff.H[i][j]);
        }
    }
    avgO2 /= static_cast<double>(diff.O2.size()*diff.O2.size());
    avgGlu /= static_cast<double>(diff.O2.size()*diff.O2.size());

    int hc = 0;
    int necC = 0;
    int prolC = 0;
    double avgPG = 0;
    double avgAT = 0;
    double avgProlProb = 0;
    for(auto &cell : cc_list){
        if(cell.hypoxic){
            hc++;
        }
        if(cell.state  == "necrotic"){
            necC++;
        }
        if(cell.state == "proliferative"){
            prolC++;
        }
        avgPG += cell.pg;
        avgAT += cell.acidThreshold;
        avgProlProb += cell.prolProb;//cell.prolThresh();
    }
    avgPG /= cc_list.size();
    avgAT /= cc_list.size();
    avgProlProb /= cc_list.size();

    std::cout << "************************************\n"
              << "Time (d): " << time/24 << std::endl
              << "Tumor Diameter: " << tumorD << std::endl
              << "Cancer Cells: " << cc_list.size() << std::endl
              << "Proliferative cells: " << prolC << std::endl
              << "Hypoxic cancer cells: " << hc << std::endl
              << "Necrotic cancer cells: " << necC << std::endl
              << "Vessels: " << vessel_list.size() << std::endl
              << "Max VEGF: " << maxVEGF << std::endl
              << "Min pH: " << -log10(maxH/(1e3)) << std::endl
              << "Max pH: " << -log10(minH/1e3) << std::endl
              << "Avg O2: " << avgO2 << std::endl
              <<  "Max O2: " << maxO2 << std::endl
              << "Avg Glu: " << avgGlu << std::endl
              << "Max Glu: " << maxGlu << std::endl
              << "Avg pG: " << avgPG << std::endl
              << "Avg acidThreshold: " << avgAT << std::endl
              << "Avg prolProb: " << avgProlProb << std::endl;
}

void Environment::updateTimeSeries(Diffusibles &diff) {
    /*
     * saves timecourse data
     */
    tumorD = tumorDiameter();
    diameterTS.push_back(tumorD);
    cancerTS.push_back(cc_list.size());
    vesselTS.push_back(vessel_list.size());

    int necC = 0;
    int hc = 0;
    int prolC = 0;
    int quiC = 0;
    for(auto &cell : cc_list){
        if(cell.hypoxic){hc++;}
        if(cell.metabState == "proliferative"){prolC++;}
        if(cell.metabState == "quiescent"){quiC++;}
        if(cell.state == "necrotic"){necC++;}
    }

    hypoxicCancerTS.push_back(hc);
    necroticCancerTS.push_back(necC);
    prolCancerTS.push_back(prolC);
    quiCancerTS.push_back(quiC);

    double minO2 = 0.056;
    double minGlu = 5.0;
    double minpH = 7.4;
    for(int i=0; i<diff.O2.size(); ++i){
	for(int j=0; j<diff.O2.size(); ++j){
	    minO2 = std::min(minO2, diff.O2[i][j]);
	    minGlu = std::min(minGlu, diff.glu[i][j]);
	    minpH = std::min(minpH, -log10(diff.H[i][j]*1e-3));
	}
    }
    o2TS.push_back(minO2);
    gluTS.push_back(minGlu);
    phTS.push_back(minpH);
}
