import numpy as np
import sys
import os

cancerParams = np.zeros((32,1))
## force parameters
cancerParams[0] = 50 # mu
cancerParams[1] = 10 # kc
cancerParams[2] = 10 # damping
cancerParams[3] = 0.2 # overlap
## metabolism parameters
cancerParams[4] = 0.012 # Vo
cancerParams[5] = 0.005 # kO
cancerParams[6] = 0.04 # kG
cancerParams[7] = 29*cancerParams[4]/5 # Ao
cancerParams[8] = 1 # pg
cancerParams[9] = 2.5e-4 # kH
## age and div
cancerParams[10] = 30 # div time (hrs)
cancerParams[11] = 1e3 # lifespan (days)
cancerParams[12] = 1e3 # maximum number of divisions
## other
cancerParams[13] = 6 # engagement time (hrs)
cancerParams[14] = 0.5 # necrotic radius decrease
cancerParams[15] = 6.5 # acidThreshold
## hypoxic effects
cancerParams[16] = 0.005 # pgShift
cancerParams[17] = 5.0 # maxPG
cancerParams[18] = 0.01 # acidThresholdShift
cancerParams[19] = 6.0 # minAcidThreshold
cancerParams[20] = 0.003 # prolProbShift
## migration
cancerParams[21] = 8.3 # um/hr -> Norton, Modeling triple-negative breast cancer heterogeneity: Effects of stromal macrophages, fibroblasts and tumor vasculature, JTB, 2018. MB231 CCR5+
cancerParams[22] = 50 # migration bias
##
cancerParams[23] = 10 # dampingFactor
cancerParams[24] = 0.001 # prolProbBound
## chemotherapy resistance
cancerParams[25] = 1 # acid-inactivation of chemotherapeutic NOT USED IN THIS STUDY. IGNORE
cancerParams[26] = 10**(-4.5) # chemoTolRate -> rate at which the cells gains tolerance
cancerParams[27] = 0.02*cancerParams[10] # chemoTimeThresh -> fraction of cell cycle cell needs to be exposed to drug for to induce tolerance
cancerParams[28] = 0.01 # chemoAccThresh -> concentration of drug needed to induce tolerance
cancerParams[29] = 90 # chemoTolerance -> base amount of drug damage to kill a cell
cancerParams[30] = 0.007 # chemoUptake -> rate of drug uptake
cancerParams[31] = 0.00001 # chemoRepair -> rate of damage repair

cellDiameters = np.zeros((2,1))
cellDiameters[0] = 15 # cancer
cellDiameters[1] = 15 # vessels

envParams = np.zeros((6,1))
envParams[0] = 158 # initial vessel spacing
envParams[1] = 2e-5 # base vessel recruitment rate. in study: [1e-5, 2e-5]
envParams[2] = 0.5 # vessel max overlap
envParams[3] = 24 # vesAge -> how long a vessel has to be alive before it can be removed
envParams[4] = 100 # healthy vessel mu
envParams[5] = 0.25 # decrease factor for unstable vessel mu and nutrient conc. in study: [0.25, 0.5]

# start at 1240
treatParams = np.zeros((8,1))
treatParams[0] = 1 # angioTreatment (1 = no treatment). Set very high (1e6) when turning treatment on
treatParams[1] = 20000 # tumor diameter at start of angioTreatment -> if no treatment, set beyond treatment duration or else vessel health will be affected
treatParams[2] = 1000 # days anio treatment is on
treatParams[3] = 0 # days angio treatment is off
treatParams[4] = 20000 # tumor diameter at start of chemotherapy
treatParams[5] = 5 # days chemotherapy is on
treatParams[6] = 10 # days chemotherapy is off
treatParams[7] = 1560 # tumor diameter at end of simulation

diffParams = np.zeros(((9,1)))
diffParams[0] = 0.012 # Vo
diffParams[1] = 0.005 # kO
diffParams[2] = 0.04 # kG
diffParams[3] = 29*diffParams[0]/5 # Ao
diffParams[4] = 1 # pg
diffParams[5] = 0.5e-4 # kH
diffParams[6] = cancerParams[30] # chemo uptake
diffParams[7] = 0.005 # chemoDecay
diffParams[8] = 50 # chemoDiff coeff

os.system('mkdir -p '+sys.argv[1]+'/params')

np.savetxt(sys.argv[1]+'/params/cancerParams.csv', cancerParams, delimiter=',')
np.savetxt(sys.argv[1]+'/params/cellDiameters.csv', cellDiameters, delimiter=',')
np.savetxt(sys.argv[1]+'/params/envParams.csv', envParams, delimiter=',')
np.savetxt(sys.argv[1]+'/params/diffParams.csv', diffParams, delimiter=',')
np.savetxt(sys.argv[1]+'/params/treatParams.csv', treatParams, delimiter=',')
