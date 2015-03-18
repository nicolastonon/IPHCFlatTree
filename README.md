IPHCFlatTree
============

IPHC analysis framework based on FlatTree

Install
-------

```
# Release to use (to be changed according to your need and CMSSW evolution)
RELEASE=7_3_0

# Setup release
cmsrel CMSSW_$RELEASE
cd CMSSW_X_Y_Z/src
cmsrel
git cms-init

# Clone this repo
git clone https://github.com/IPHC/IPHCFlatTree.git

# Switch to release
cd IPHCFlatTree
git checkout v20150314_patch1
cd ../

# Add the dependencies
git cms-addpkg EgammaAnalysis/ElectronTool
cd EgammaAnalysis/ElectronTool/data
for FILE in `cat download.url`; do wget $FILE; done;
cd ../../..
```
