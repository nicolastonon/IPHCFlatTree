IPHCFlatTree
============

IPHC analysis framework based on FlatTree

Install
-------

```
# CMSSW Release
RELEASE=8_0_4

# Setup release
cmsrel CMSSW_$RELEASE
cd CMSSW_X_Y_Z/src
cmsenv
git cms-init

# Clone this repo
git clone https://github.com/IPHC/IPHCFlatTree.git

# Tools needed for AK10 jet collection
git clone https://github.com/cms-jet/JetToolbox JMEAnalysis/JetToolbox 

# Compile the monster
scram b
```


Last modif
----------
- add a new section in conf.xml in order to apply a preselection (filtering events)
- preselection based on n_muons, n_electrons (or their sum), n_jets, MET passing pt, eta cuts

To do
---------
- preselection does not include taus, b-jets, trigger, ...
- both pt and eta cuts for objects need to be float - one should pay attention in the config.xml file
