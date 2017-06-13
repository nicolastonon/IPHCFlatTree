# CMSSW Release
RELEASE=8_0_25

# Setup release
cmsrel CMSSW_$RELEASE
cd CMSSW_8_0_25/src
cmsenv
git cms-init

# Clone this repo
git clone https://github.com/IPHC/IPHCFlatTree.git

# Egamma
git cms-merge-topic shervin86:Moriond2017_JEC_energyScales
cd EgammaAnalysis/ElectronTools/data; git clone git@github.com:ECALELFS/ScalesSmearings.git; cd -
git cms-merge-topic ikrav:egm_id_80X_v2

# DeepCSV
cms-merge-topic mverzett:DeepFlavour-from-CMSSW_8_0_21 
mkdir RecoBTag/DeepFlavour/data/ 
cd RecoBTag/DeepFlavour/data/
wget http://home.fnal.gov/~verzetti//DeepFlavour/training/DeepFlavourNoSL.json
cd - 

# Add MET filters
git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter

# Tools needed for AK10 jet collection
git clone https://github.com/cms-jet/JetToolbox JMEAnalysis/JetToolbox 

# Compile the monster
scram b
