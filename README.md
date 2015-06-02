IPHCFlatTree
============

IPHC analysis framework based on FlatTree

Install
-------

```
# CMSSW Release
RELEASE=7_2_3

# Setup release
cmsrel CMSSW_$RELEASE
cd CMSSW_X_Y_Z/src
cmsenv
git cms-init

# Clone this repo
git clone https://github.com/IPHC/IPHCFlatTree.git

# Switch to particular release, otherwise ignore these lines to use the HEAD
cd IPHCFlatTree
git checkout MantaRay-patch1
cd ../

# Add the dependencies
git cms-merge-topic HuguesBrun:trigElecIdInCommonIsoSelection720
```
