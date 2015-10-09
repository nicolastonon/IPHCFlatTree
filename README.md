IPHCFlatTree
============

IPHC analysis framework based on FlatTree

Install
-------

```
# CMSSW Release
RELEASE=7_4_12_patch4

# Setup release
cmsrel CMSSW_$RELEASE
cd CMSSW_X_Y_Z/src
cmsenv
git cms-init

# Clone this repo
git clone https://github.com/IPHC/IPHCFlatTree.git

# Egamma ID
git cms-merge-topic ikrav:egm_id_7.4.12_v1

# Switch to particular release, otherwise ignore these lines to use the HEAD
cd IPHCFlatTree
git checkout MantaRay-patch8
cd ../

# Compile the monster
scram b
```
