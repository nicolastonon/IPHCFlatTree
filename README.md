# FlatTreeProducer installation and setup

README for the IPHCFllatTree -- tHq branch, describing the basic steps to run the FlatTree production.

*Do not forget to source :*
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

*To create proxy :*
```
voms-proxy-init -voms cms -hours 192
```

## FlatTreeProducer

Installing and running the IPHCFlatTree code to produce Flat Trees. Using branch "tHq" (based on tag "Walrus-patch-2").

### Installation

```
cd /home-pbs/username/
mkdir MyAnalysis
cd MyAnalysis

# CMSSW Release
RELEASE=8_0_25

# Setup release
cmsrel CMSSW_$RELEASE
cd CMSSW_X_Y_Z/src
cmsenv
git cms-init

# Clone this repo (tHq branch!)
git clone -b tHq https://github.com/IPHC/IPHCFlatTree.git

# Egamma
git cms-merge-topic shervin86:Moriond2017_JEC_energyScales
cd EgammaAnalysis/ElectronTools/data; git clone https://github.com/ECALELFS/ScalesSmearings; cd -
git cms-merge-topic ikrav:egm_id_80X_v2

# Add MET filters
git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter

# Tools needed for AK10 jet collection
git clone https://github.com/cms-jet/JetToolbox JMEAnalysis/JetToolbox 

# Add DeepCSV tagger
git cms-merge-topic -u mverzett:DeepFlavour-from-CMSSW_8_0_21
mkdir RecoBTag/DeepFlavour/data/; cd RecoBTag/DeepFlavour/data/; wget http://home.fnal.gov/~verzetti//DeepFlavour/training/DeepFlavourNoSL.json; cd -

# Compile the monster (use -jN for multicore)
scram b -j5
```

(( Instructions taken from [IPHCFlatTree's README](https://github.com/IPHC/IPHCFlatTree/tree/Walrus-patch2) ))


### Set-up


```
cd XXX/IPHCFlatTree/FlatTreeProducer/test/PROD
```
* **list.txt** - create it and add the names of the datasets/samples you want to process, e.g. : 
```
...
/THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
...
```
(NB : if you add several samples at once, they will each yield a separate task. The merging of files, e.g. of several extensions of one sample, has to be done at the NTupleProducer level)



* **submit.zsh** - modify the following :
```
...
ver="XXX" //Version name, e.g. "tHqProd"
prodv="/store/user/YOUR_USERNAME/FlatTree/${ver}/" //Will store output files on dpm/store
...
```

* **crabConfigTemplate.py** - modify the following :
```
...
isData=0 #Or 1 for data
...
config.Data.unitsPerJob = 2 #For MC, ~2 MC files per job (else too large)
config.Data.unitsPerJob = 10 #For Data, ~10 lumisections per job
...
config.Data.splitting = 'FileBased' #For MC
#config.Data.splitting = 'LumiBased' #For data
...
#config.Data.lumiMask = 'GRL/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt' #Need to comment for MC prod
...
```


### Interactive test

Can use a modified cfg file in directory test/ (else, errors from some relative paths), to check that you can access the data. Some lines needs to be changed (sample name, etc.)

```
cmsRun crabConfig_test.py
```


### Launch the jobs

```
./submit.zsh
```
