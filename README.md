IPHCFlatTree
============

IPHC analysis framework based on FlatTree

Install
-------

```
# CMSSW Release
RELEASE=7_4_2

# Setup release
cmsrel CMSSW_$RELEASE
cd CMSSW_X_Y_Z/src
cmsenv
git cms-init

# Add the dependencies for POG CMSSW packages
git cms-addpkg EgammaAnalysis/ElectronTools
cd EgammaAnalysis/ElectronTools/data
cat download.url | xargs wget
cd -

# Clone this repo
git clone https://github.com/IPHC/IPHCFlatTree.git

# Switch to particular release, otherwise ignore these lines to use the HEAD
# IGNORE THIS STEP FOR 74X
cd IPHCFlatTree
git checkout MantaRay-patch4
cd ../

# Compile the monster
scram b
```
