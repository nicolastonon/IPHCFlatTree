#!/bin/env zsh

wpath="https://github.com/CERN-PH-CMG/cmg-cmssw/blob/CMG_MiniAOD_Lite_V6_0_from-CMSSW_7_0_6/CMGTools/TTHAnalysis/data/leptonMVA"

wget "${wpath}/el_pteta_high_cb_BDTG.weights.xml?raw=true"
mv el_pteta_high_cb_BDTG.weights.xml\?raw\=true el_pteta_high_cb_BDTG.weights.xml

wget "${wpath}/el_pteta_high_ec_BDTG.weights.xml?raw=true"
mv el_pteta_high_ec_BDTG.weights.xml\?raw\=true el_pteta_high_ec_BDTG.weights.xml

wget "${wpath}/el_pteta_high_fb_BDTG.weights.xml?raw=true"
mv el_pteta_high_fb_BDTG.weights.xml\?raw\=true el_pteta_high_fb_BDTG.weights.xml

wget "${wpath}/el_pteta_low_cb_BDTG.weights.xml?raw=true"
mv el_pteta_low_cb_BDTG.weights.xml\?raw\=true el_pteta_low_cb_BDTG.weights.xml

wget "${wpath}/el_pteta_low_ec_BDTG.weights.xml?raw=true"
mv el_pteta_low_ec_BDTG.weights.xml\?raw\=true el_pteta_low_ec_BDTG.weights.xml

wget "${wpath}/el_pteta_low_fb_BDTG.weights.xml?raw=true"
mv el_pteta_low_fb_BDTG.weights.xml\?raw\=true el_pteta_low_fb_BDTG.weights.xml

wget "${wpath}/mu_pteta_high_b_BDTG.weights.xml?raw=true"
mv mu_pteta_high_b_BDTG.weights.xml\?raw\=true mu_pteta_high_b_BDTG.weights.xml

wget "${wpath}/mu_pteta_high_e_BDTG.weights.xml?raw=true"
mv mu_pteta_high_e_BDTG.weights.xml\?raw\=true mu_pteta_high_e_BDTG.weights.xml

wget "${wpath}/mu_pteta_low_b_BDTG.weights.xml?raw=true"
mv mu_pteta_low_b_BDTG.weights.xml\?raw\=true mu_pteta_low_b_BDTG.weights.xml

wget "${wpath}/mu_pteta_low_e_BDTG.weights.xml?raw=true"
mv mu_pteta_low_e_BDTG.weights.xml\?raw\=true mu_pteta_low_e_BDTG.weights.xml
