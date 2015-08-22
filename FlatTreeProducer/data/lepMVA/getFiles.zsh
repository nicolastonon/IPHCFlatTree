#!/bin/env zsh

# use blob in the link
wpath="https://github.com/CERN-PH-CMG/cmg-cmssw/blob/CMGTools-from-CMSSW_7_4_7/CMGTools/TTHAnalysis/data/leptonMVA/tth"

wget "${wpath}/el_pteta_high_cb_BDTG.weights.xml?raw=true"
mv el_pteta_high_cb_BDTG.weights.xml\?raw\=true el_pteta_high_cb_BDTG.weights.xml

wget "${wpath}/el_pteta_high_ec_BDTG.weights.xml?raw=true"
mv el_pteta_high_ec_BDTG.weights.xml\?raw\=true el_pteta_high_ec_BDTG.weights.xml

wget "${wpath}/el_pteta_high_fb_BDTG.weights.xml?raw=true"
mv el_pteta_high_fb_BDTG.weights.xml\?raw\=true el_pteta_high_fb_BDTG.weights.xml

wget "${wpath}/el_pteta_low_BDTG.weights.xml?raw=true"
mv el_pteta_low_BDTG.weights.xml\?raw\=true el_pteta_low_BDTG.weights.xml

wget "${wpath}/el_pteta_medium_cb_BDTG.weights.xml?raw=true"
mv el_pteta_medium_cb_BDTG.weights.xml\?raw\=true el_pteta_medium_cb_BDTG.weights.xml

wget "${wpath}/el_pteta_medium_ec_BDTG.weights.xml?raw=true"
mv el_pteta_medium_ec_BDTG.weights.xml\?raw\=true el_pteta_medium_ec_BDTG.weights.xml

wget "${wpath}/el_pteta_medium_fb_BDTG.weights.xml?raw=true"
mv el_pteta_medium_fb_BDTG.weights.xml\?raw\=true el_pteta_medium_fb_BDTG.weights.xml

wget "${wpath}/mu_pteta_high_b_BDTG.weights.xml?raw=true"
mv mu_pteta_high_b_BDTG.weights.xml\?raw\=true mu_pteta_high_b_BDTG.weights.xml

wget "${wpath}/mu_pteta_high_e_BDTG.weights.xml?raw=true"
mv mu_pteta_high_e_BDTG.weights.xml\?raw\=true mu_pteta_high_e_BDTG.weights.xml

wget "${wpath}/mu_pteta_low_BDTG.weights.xml?raw=true"
mv mu_pteta_low_BDTG.weights.xml\?raw\=true mu_pteta_low_BDTG.weights.xml

wget "${wpath}/mu_pteta_medium_b_BDTG.weights.xml?raw=true"
mv mu_pteta_medium_b_BDTG.weights.xml\?raw\=true mu_pteta_medium_b_BDTG.weights.xml

wget "${wpath}/mu_pteta_medium_e_BDTG.weights.xml?raw=true"
mv mu_pteta_medium_e_BDTG.weights.xml\?raw\=true mu_pteta_medium_e_BDTG.weights.xml

