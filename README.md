IPHCFlatTree
============

IPHC analysis framework based on FlatTree

Install
-------

```
# CMSSW Release
RELEASE=8_0_12

# Setup release
cmsrel CMSSW_$RELEASE
cd CMSSW_X_Y_Z/src
cmsenv
git cms-init

# Clone this repo
git clone https://github.com/IPHC/IPHCFlatTree.git

----------------------------------------------------------------

# For EGMSmearer

based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer

git remote add -f -t ecal_smear_fix_80X emanueledimarco https://github.com/emanueledimarco/cmssw.git
git cms-addpkg EgammaAnalysis/ElectronTools
git checkout -b from-277de3c 277de3c

cd EgammaAnalysis/ElectronTools/data
git clone -b ICHEP2016_approval_7p65fb https://github.com/emanueledimarco/ScalesSmearings.git
----------------------------------------------------------------


# Add MET filters
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate

# Tools needed for AK10 jet collection
git clone https://github.com/cms-jet/JetToolbox JMEAnalysis/JetToolbox 

# Compile the monster
scram b
```


Last modif
----------
- add a new section in conf.xml in order to apply a preselection (filtering events)
- to activate it, one needs to run on the boolean 'activate' in the preselection section
- preselection based on n_muons, n_electrons (or their sum), n_jets, MET passing pt, eta cuts
- current limitations:
  - can only use cut_min for pt and cut_max for abs(eta)

To do
---------
- preselection does not include taus, b-jets, trigger, ...
- both pt and eta cuts for objects need to be float - one should pay attention in the config.xml file
