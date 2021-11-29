#include "Diffusibles.h"

Diffusibles::Diffusibles(double size, double spatialStep, std::vector<double> &diffParams){
    dx = spatialStep;
    maxPG = 1;

    N = static_cast<int>(size/dx); // size of PDE grids

    // initialize all the grids to 0
    for(int i=0;i<N;++i){
        std::vector<double> row;
        // initialize PDE grids to zero
        for(int j=0;j<N;++j){
            row.push_back(0.0);
        }
        O2.push_back(row);
        glu.push_back(row);
        H.push_back(row);
        O2_next.push_back(row);
        glu_next.push_back(row);
        H_next.push_back(row);
        cancerCells.push_back(row);
        c.push_back(row);
        VEGF.push_back(row);
        VEGF_next.push_back(row);
        vessels.push_back(row);
        Vo.push_back(row);
        kO.push_back(row);
        Ao.push_back(row);
        kG.push_back(row);
        pg.push_back(row);
        kH.push_back(row);
        hypoxicCells.push_back(row);
        normalCells.push_back(row);
        dc.push_back(row);
        C.push_back(row);
        C_next.push_back(row);
    }
    for(int i=0; i<N; ++i){
        for(int j=0; j<N; ++j){
            // set nutrient concentrations to very low to prevent errors when initializing
            // may not be needed anymore
            O2[i][j] = 1e-10;//0;//0.056;
            glu[i][j] = 1e-10;//0;//5.0;
            H[i][j] = 0;//pow(10,-7.4)*1e3;
            // assume environment is filled with normal tissue
            normalCells[i][j] = 1;
            dc[i][j] = 1;
        }
    }

    // metabolism for normal tissue
    baseVo = diffParams[0];
    basekO = diffParams[1];
    basekG = diffParams[2];
    baseAo = diffParams[3];
    basepg = diffParams[4];
    basekH = diffParams[5];

    // diffusion coefficients um^2 / sec
    Do = 1820.0;
    Dg = 500.0;
    Dh = 1080.0;
    Dc = diffParams[8];

    // Malbacher 2018
    // from ^ paper, they scale to oxygen diffusion
    // Curtis 2020
    DVEGF = 0.01880*Do/3.7606;

    // secretion rates for cytokines. non-dimensional
    kVEGF = 1;

    // chemotherapy uptake
    kupC = diffParams[6];
    kdegC = diffParams[7];
}

void Diffusibles::runDiffusion(double time, bool chemoOn) {
    /*
     * use finite difference to calculate diffusion
     * if max % change between steps is < cutoff, assume equilibrium and end early
     */

    double cutoff = 1e-4;//1e-6;

    // set time step based on Do (the largest diffusion coeff) and pg (glucose uptake)
    double dt = 0.2*dx*dx/(Do*maxPG);
    double steps = std::max(1.0, (time*60*60)/dt);

    // diffuse  nutrients
    for(int q=0; q<steps; ++q) {
#pragma omp parallel for collapse(2)
        for (int i = 1; i < O2.size() - 1; ++i) {
            for (int j = 1; j < O2.size() - 1; ++j) {
                double maxDiff = 0;

                // METABOLSIM
                finDiff(O2, O2_next, i, j, dt, Do);
                //double fO = (-Vo[i][j] * O2[i][j] / (kO[i][j] + O2[i][j]));
                double fO = (-Vo[i][j] * O2_next[i][j] / (kO[i][j] + O2_next[i][j]));
                O2_next[i][j] += dt * fO;

                finDiff(glu, glu_next, i, j, dt, Dg);
                //double fG = (-(pg[i][j] * Ao[i][j] / 2 + 27 * fO / 10) * glu[i][j] / (kG[i][j] + glu[i][j]));
                double fG = (-(pg[i][j] * Ao[i][j] / 2 + 27 * fO / 10) * glu_next[i][j] / (kG[i][j] + glu_next[i][j]));
                glu_next[i][j] += dt * fG;

                finDiff(H, H_next, i, j, dt, Dh);
                double fH = (-kH[i][j] * 2 * fG);
                H_next[i][j] += dt * fH;

                // for healthy vessels, keep nutrient concentration constant at blood levels
                // for unhealthy vessels, slow the rate of nutrient replacement
                //  this should reach a steady-state with a lower concentration than blood
                if(vessels[i][j] == 1){
                    O2_next[i][j] += dc[i][j]*(0.056 - O2_next[i][j]);
                    glu_next[i][j] += dc[i][j]*(5.0 - glu_next[i][j]);
                    H_next[i][j] += dc[i][j]*((pow(10,-7.4)*1e3) - H_next[i][j]);
                    /*O2_next[i][j] += dc[i][j]*(0.056 - O2[i][j]);
                    glu_next[i][j] += dc[i][j]*(5.0 - glu[i][j]);
                    H_next[i][j] += dc[i][j]*((pow(10,-7.4)*1e3) - H[i][j]);*/
                }

                // FIND MAX % DIFFERENCE AT EACH SITE
                if(O2[i][j] > 0){
                    maxDiff = std::max(maxDiff, fabs(O2_next[i][j] - O2[i][j])/O2[i][j]);
                }
                if(glu[i][j] > 0){
                    maxDiff = std::max(maxDiff, fabs(glu_next[i][j] - glu[i][j])/glu[i][j]);
                }
                if(H[i][j] > 0){
                    maxDiff = std::max(maxDiff, fabs(H_next[i][j] - H[i][j])/H[i][j]);
                }
                c[i][j] = maxDiff;
            }
        }

        // find total max difference
        double maxDiff = 0;
        for(int i=0; i<N; ++i){
            for(int j=0; j<N; ++j){
                maxDiff = std::max(maxDiff, c[i][j]);
                if(O2_next[i][j] < 0){
                    std::cout << "O2 negative\n";
                    throw std::runtime_error(" ");
                }
                if(glu_next[i][j] < 0){
                    std::cout << "glu negative\n";
                    throw std::runtime_error(" ");
                }
                if(H_next[i][j] < 0){
                    std::cout << "H negative\n";
                    throw std::runtime_error(" ");
                }
            }
        }

#pragma omp parallel for collapse(2)
        for (int i = 1; i < O2.size() - 1; ++i) {
            for (int j = 1; j < O2.size() - 1; ++j) {
                // update grids
                O2[i][j] = O2_next[i][j];
                glu[i][j] = glu_next[i][j];
                H[i][j] = H_next[i][j];
            }
        }
        if(maxDiff < cutoff && q > 10){
            //std::cout << "nutrient diff: " << maxDiff << std::endl; 
            break;}
    }

    // diffuse cytokines
    dt = 0.2*dx*dx/DVEGF;
    dt = std::min(1.0, dt);
    steps = std::max(1.0, (time*60*60)/dt);
    double stepsTaken = 0;
    for(int q=0; q<steps; ++q) {
#pragma omp parallel for collapse(2)
        for (int i = 1; i < O2.size() - 1; ++i) {
            for (int j = 1; j < O2.size() - 1; ++j) {
                double maxDiff = 0;
                // CYTOKINES
                // scale values from 0 - 1. degredation in environment. removal via vessels
                finDiff(VEGF, VEGF_next, i, j, dt, DVEGF);
                VEGF_next[i][j] += dt*(kVEGF*hypoxicCells[i][j]*(1 - VEGF_next[i][j]) - 0.001*VEGF[i][j] - 0.006*VEGF[i][j]*vessels[i][j]);

                // FIND MAX % DIFFERENCE AT EACH SITE
                if(VEGF[i][j] > 0){
                    maxDiff = std::max(maxDiff, fabs(VEGF_next[i][j] - VEGF[i][j])/VEGF[i][j]);
                }
                c[i][j] = maxDiff;
            }
        }

        // find total max difference
        double maxDiff = 0;
        for(int i=0; i<N; ++i){
            for(int j=0; j<N; ++j){
                maxDiff = std::max(maxDiff, c[i][j]);
                if(VEGF_next[i][j] < 0){
                    std::cout << "VEGF negative\n";
                    throw std::runtime_error(" ");
                }
            }
        }

#pragma omp parallel for collapse(2)
        for (int i = 1; i < O2.size() - 1; ++i) {
            for (int j = 1; j < O2.size() - 1; ++j) {
                // update grids
                VEGF[i][j] = VEGF_next[i][j];
            }
        }

        stepsTaken++;
        if(maxDiff < cutoff && q > 10){
            //std::cout << "vegf diff: " << maxDiff << std::endl;
            break;}
    }


    // diffuse chemotherapy
    dt = 0.2*dx*dx/Dc;
    steps = std::max(1.0, (time*60*60)/dt);
    for(int q=0; q<steps; ++q) {
#pragma omp parallel for collapse(2)
        for (int i = 1; i < O2.size() - 1; ++i) {
            for (int j = 1; j < O2.size() - 1; ++j) {
                double maxDiff = 0;
                finDiff(C, C_next, i, j, dt, Dc);
                C_next[i][j] += dt*(-kupC*C_next[i][j]*cancerCells[i][j] - kdegC*C_next[i][j]);

                if(vessels[i][j] == 1 && chemoOn){
                    C_next[i][j] += dc[i][j]*(1.0 - C_next[i][j]);
                }

                // FIND MAX % DIFFERENCE AT EACH SITE
                if(C[i][j] > 0){
                    maxDiff = std::max(maxDiff, fabs(C_next[i][j] - C[i][j])/C[i][j]);
                }
                c[i][j] = maxDiff;
            }
        }

        // find total max difference
        double maxDiff = 0;
        for(int i=0; i<N; ++i){
            for(int j=0; j<N; ++j){
                maxDiff = std::max(maxDiff, c[i][j]);
                if(C_next[i][j] < 0){
                    C_next[i][j] = 0;
                }
            }
        }

#pragma omp parallel for collapse(2)
        for (int i = 1; i < O2.size() - 1; ++i) {
            for (int j = 1; j < O2.size() - 1; ++j) {
                // update grids
                C[i][j] = C_next[i][j];
            }
        }
        if(maxDiff < cutoff && q > 10){break;}
    }
}

void Diffusibles::finDiff(std::vector<std::vector<double>> &g,
                                                      std::vector<std::vector<double>> &g_next, int i, int j, double dt,
                                                      double D) {
    g_next[i][j] = g[i][j] + dt*D/(dx*dx)*(g[i-1][j] + g[i+1][j] + g[i][j-1] + g[i][j+1] - 4*g[i][j]);
}

void Diffusibles::vesselLocations(std::vector<Vessel> vessel_list) {
    /*
     * set spots that contain vessels
     */
    for(int i=0; i<N; ++i){
        for(int j=0; j<N; ++j){
            vessels[i][j] = 0;
            dc[i][j] = 1;
        }
    }
    for(auto &vessel : vessel_list){
        int i = vessel.x[0] + dx/2;
        i -= fmod(i, dx);
        i /= dx;
        int j = vessel.x[1] + dx/2;
        j -= fmod(j, dx);
        j /= dx;

        // if vessels are larger than dx, they take up multiple spots
        int n = vessel.radius/dx;
        for(int s=0; s<n+1; ++s){
            vessels[i+s][j+s] = 1;
            normalCells[i+s][j+s] = 0;
            dc[i+s][j+s] = vessel.health;
        }
    }
}

void Diffusibles::setBounds() {
    /*
     * run at start of simulation
     * sets the bounds to the average values for each nutrient
     */
    double dt = 0.2*dx*dx/Do;

    bool go = true;
    while(go) {
#pragma omp parallel for
        for (int i = 1; i < O2.size() - 1; ++i) {
            for (int j = 1; j < O2.size() - 1; ++j) {
                double maxDiff = 0;

                finDiff(O2, O2_next, i, j, dt, Do);
                double fO = (-Vo[i][j] * O2[i][j] / (kO[i][j] + O2[i][j]));
                O2_next[i][j] += dt * fO;

                finDiff(glu, glu_next, i, j, dt, Dg);
                double fG = (-(pg[i][j] * Ao[i][j] / 2 + 27 * fO / 10) * glu[i][j] / (kG[i][j] + glu[i][j]));
                glu_next[i][j] += dt * fG;

                finDiff(H, H_next, i, j, dt, Dh);
                double fH = (-kH[i][j] * 2 * fG);
                H_next[i][j] += dt * fH;

                if(vessels[i][j] == 1){
                    O2_next[i][j] = 0.056;
                    glu_next[i][j] = 5.0;
                    H_next[i][j] = pow(10,-7.4)*1e3;
                }

                if(O2[i][j] > 0){
                    maxDiff = std::max(maxDiff, fabs(O2_next[i][j] - O2[i][j])/O2[i][j]);
                }
                if(glu[i][j] > 0){
                    maxDiff = std::max(maxDiff, fabs(glu_next[i][j] - glu[i][j])/glu[i][j]);
                }
                if(H[i][j] > 0){
                    maxDiff = std::max(maxDiff, fabs(H_next[i][j] - H[i][j])/H[i][j]);
                }
                c[i][j] = maxDiff;
            }
        }

        double maxDiff = 0;
        double avgO2 = 0;
        double avgGlu = 0;
        double avgH = 0;
        for(int i=0; i<N; ++i){
            for(int j=0; j<N; ++j){
                maxDiff = std::max(maxDiff, c[i][j]);
                if(O2_next[i][j] < 0){
                    std::cout << "O2 negative\n";
                    throw std::runtime_error(" ");
                }
                if(glu_next[i][j] < 0){
                    std::cout << "glu negative\n";
                    throw std::runtime_error(" ");
                }
                if(H_next[i][j] < 0){
                    std::cout << "H negative\n";
                    throw std::runtime_error(" ");
                }
                avgO2 += O2_next[i][j];
                avgGlu += glu_next[i][j];
                avgH += H_next[i][j];
            }
        }
        avgO2 /= static_cast<double>(O2.size()*O2.size());
        avgGlu /= static_cast<double>(O2.size()*O2.size());
        avgH /= static_cast<double>(O2.size()*O2.size());

        for(int i=0; i<O2.size(); ++i){
            O2_next[0][i] = avgO2;
            O2_next[i][0] = avgO2;
            O2_next[N-1][i] = avgO2;
            O2_next[i][N-1] = avgO2;

            glu_next[0][i] = avgGlu;
            glu_next[i][0] = avgGlu;
            glu_next[N-1][i] = avgGlu;
            glu_next[i][N-1] = avgGlu;

            H_next[0][i] = avgH;
            H_next[i][0] = avgH;
            H_next[N-1][i] = avgH;
            H_next[i][N-1] = avgH;
        }

#pragma omp parallel for
        for (int i = 0; i < O2.size(); ++i) {
            for (int j = 0; j < O2.size(); ++j) {
                O2[i][j] = O2_next[i][j];
                glu[i][j] = glu_next[i][j];
                H[i][j] = H_next[i][j];
            }
        }
        if(maxDiff < 1e-4){go = false;}
    }
}

void Diffusibles::setMaxPG(double cellPG) {
    /*
     * determine the maximum shift in glucose uptake
     * makes the time step small enough to prevent negative glucose concentrations
     */
    maxPG = std::max(maxPG, cellPG);
}
