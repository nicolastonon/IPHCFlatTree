IPHCFlatTree
============

IPHC analysis framework based on FlatTree

Install
-------

```
# CMSSW Release
RELEASE=7_6_4

# Setup release
cmsrel CMSSW_$RELEASE
cd CMSSW_X_Y_Z/src
cmsenv
git cms-init

# Clone this repo
git clone https://github.com/IPHC/IPHCFlatTree.git

# Tools needed for AK10 jet collection
git clone https://github.com/cms-jet/JetToolbox JMEAnalysis/JetToolbox 

# Switch to particular release (if you don't want to use the HEAD version)
cd IPHCFlatTree
git checkout Akoula-patch3
cd ../

# Compile the monster
scram b
```
