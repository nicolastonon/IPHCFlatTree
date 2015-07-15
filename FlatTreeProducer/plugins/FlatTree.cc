#include "IPHCFlatTree/FlatTreeProducer/interface/FlatTree.hh"

void FlatTree::Init()
{
   n_presel_jets = 0;
   n_presel_btag = 0;
   n_presel_electron = 0;
   n_presel_muon = 0;
   n_presel_tau = 0;

   ev_run = DEFVAL;
   ev_id = DEFVAL;
   ev_lumi = DEFVAL;
   ev_rho = DEFVAL;

   met_pt = DEFVAL;
   met_phi = DEFVAL;
   met_sumet = DEFVAL;

   metPuppi_pt = DEFVAL;
   metPuppi_phi = DEFVAL;
   metPuppi_sumet = DEFVAL;
   
   pv_x = DEFVAL;
   pv_y = DEFVAL;
   pv_z = DEFVAL;

   pv_ndof = DEFVAL;
   pv_rho = DEFVAL;
   pv_isFake = DEFVAL;
   
   mc_weight = DEFVAL;
   mc_id = DEFVAL;
   mc_f1 = DEFVAL;
   mc_f2 = DEFVAL;
   mc_x1 = DEFVAL;
   mc_x2 = DEFVAL;
   mc_scale = DEFVAL;
   mc_ptHat = DEFVAL;

   mc_pu_intime_NumInt = DEFVAL;
   mc_pu_trueNumInt = DEFVAL;
   mc_pu_before_npu = DEFVAL;
   mc_pu_after_npu = DEFVAL;

   mc_pu_Npvi = DEFVAL;
   mc_pu_Nzpositions.clear();
   mc_pu_BunchCrossing.clear();
   for(unsigned int i=0;i<mc_pu_zpositions.size();i++) mc_pu_zpositions[i].clear();
   mc_pu_zpositions.clear();
   for(unsigned int i=0;i<mc_pu_sumpT_lowpT.size();i++) mc_pu_sumpT_lowpT[i].clear();
   mc_pu_sumpT_lowpT.clear();
   for(unsigned int i=0;i<mc_pu_sumpT_highpT.size();i++) mc_pu_sumpT_highpT[i].clear();
   mc_pu_sumpT_highpT.clear();
   for(unsigned int i=0;i<mc_pu_ntrks_lowpT.size();i++) mc_pu_ntrks_lowpT[i].clear();
   mc_pu_ntrks_lowpT.clear();
   for(unsigned int i=0;i<mc_pu_ntrks_highpT.size();i++) mc_pu_ntrks_highpT[i].clear();
   mc_pu_ntrks_highpT.clear();

   trigger_n = 0;
   trigger.clear();
   trigger_name.clear();
   trigger_pass.clear();
   trigger_prescale.clear();

   triggerobject_n = 0;
   triggerobject_pt.clear();
   triggerobject_eta.clear();
   triggerobject_phi.clear();

   triggerobject_collection.clear();

   triggerobject_filterIds_n.clear();
   triggerobject_filterIds.clear();

   triggerobject_isTriggerL1Mu.clear();
   triggerobject_isTriggerL1NoIsoEG.clear();
   triggerobject_isTriggerL1IsoEG.clear();
   triggerobject_isTriggerL1CenJet.clear();
   triggerobject_isTriggerL1ForJet.clear();
   triggerobject_isTriggerL1TauJet.clear();
   triggerobject_isTriggerL1ETM.clear();
   triggerobject_isTriggerL1ETT.clear();
   triggerobject_isTriggerL1HTT.clear();
   triggerobject_isTriggerL1HTM.clear();
   triggerobject_isTriggerL1JetCounts.clear();
   triggerobject_isTriggerL1HfBitCounts.clear();
   triggerobject_isTriggerL1HfRingEtSums.clear();
   triggerobject_isTriggerL1TechTrig.clear();
   triggerobject_isTriggerL1Castor.clear();
   triggerobject_isTriggerL1BPTX.clear();
   triggerobject_isTriggerL1GtExternal.clear();

   triggerobject_isHLT_TriggerPhoton.clear();
   triggerobject_isHLT_TriggerElectron.clear();
   triggerobject_isHLT_TriggerMuon.clear();
   triggerobject_isHLT_TriggerTau.clear();
   triggerobject_isHLT_TriggerJet.clear();
   triggerobject_isHLT_TriggerBJet.clear();
   triggerobject_isHLT_TriggerMET.clear();
   triggerobject_isHLT_TriggerTET.clear();
   triggerobject_isHLT_TriggerTHT.clear();
   triggerobject_isHLT_TriggerMHT.clear();
   triggerobject_isHLT_TriggerTrack.clear();
   triggerobject_isHLT_TriggerCluster.clear();
   triggerobject_isHLT_TriggerMETSig.clear();
   triggerobject_isHLT_TriggerELongit.clear();
   triggerobject_isHLT_TriggerMHTSig.clear();
   triggerobject_isHLT_TriggerHLongit.clear();

   triggerobject_filterLabels_n.clear();
   triggerobject_filterLabels.clear();

   triggerobject_pathNamesAll_n.clear();
   triggerobject_pathNamesAll.clear();
   triggerobject_pathNamesAll_isBoth.clear();
   triggerobject_pathNamesAll_isL3.clear();
   triggerobject_pathNamesAll_isLF.clear();
   triggerobject_pathNamesAll_isNone.clear();

   el_n = 0;
   el_pt.clear();
   el_eta.clear();
   el_phi.clear();
   el_m.clear();
   el_E.clear();
   el_id.clear();
   el_charge.clear();

   el_scleta.clear();
   el_isGsfCtfScPixChargeConsistent.clear();
   el_passConversionVeto.clear();
   el_dB3D.clear();
   el_edB3D.clear();

   el_neutralHadronIso.clear();
   el_chargedHadronIso.clear();
   el_puChargedHadronIso.clear();
   el_ecalIso.clear();
   el_hcalIso.clear();
   el_particleIso.clear();
   el_photonIso.clear();
   el_trackIso.clear();

   el_pfIso_sumChargedHadronPt.clear();
   el_pfIso_sumNeutralHadronEt.clear();
   el_pfIso_sumPhotonEt.clear();
   el_pfIso_sumPUPt.clear();

   el_miniIso.clear();

   el_isLoose.clear();
   el_isTight.clear();
   el_isRobustLoose.clear();
   el_isRobustTight.clear();
   el_isRobustHighEnergy.clear();

   el_vx.clear();
   el_vy.clear();
   el_vz.clear();

   el_isGsf.clear();
   el_dxy.clear();
   el_dz.clear();
   el_dxyError.clear();
   el_dzError.clear();

   el_numberOfHits.clear();

   el_sigmaIetaIeta.clear();
   el_sigmaIphiIphi.clear();
   el_hadronicOverEm.clear();
   el_dr03TkSumPt.clear();
   el_dr03EcalRecHitSumEt.clear();
   el_dr03HcalTowerSumEt.clear();
   el_numberOfLostHits.clear();

   el_fbrem.clear();
   el_kf_normalizedChi2.clear();
   el_trackerLayersWithMeasurement.clear();
   el_gsf_normalizedChi2.clear();
   el_deltaEtaSuperClusterTrackAtVtx.clear();
   el_deltaPhiSuperClusterTrackAtVtx.clear();
   el_deltaEtaSeedClusterTrackAtCalo.clear();
   el_see.clear();
   el_spp.clear();
   el_superClusterEtaWidth.clear();
   el_superClusterPhiWidth.clear();
   el_full5x5_OneMinusE1x5E5x5.clear();
   el_OneMinusE1x5E5x5.clear();
   el_full5x5_r9.clear();
   el_r9.clear();
   el_eSuperClusterOverP.clear();
   el_IoEmIoP.clear();
   el_ooEmooP.clear();
   el_eleEoPout.clear();
   el_PreShowerOverRaw.clear();
   el_ecalEnergy.clear();

   el_mvaNonTrigV0.clear();

   el_lepMVA.clear();

   el_lepMVA_neuRelIso.clear();
   el_lepMVA_chRelIso.clear();
   el_lepMVA_jetDR.clear();
   el_lepMVA_jetPtRatio.clear();
   el_lepMVA_jetBTagCSV.clear();
   el_lepMVA_sip3d.clear();
   el_lepMVA_dxy.clear();
   el_lepMVA_dz.clear();
   el_lepMVA_mvaId.clear();

   el_hasMCMatch.clear();
   el_gen_pt.clear();
   el_gen_eta.clear();
   el_gen_phi.clear();
   el_gen_m.clear();
   el_gen_status.clear();
   el_gen_id.clear();
   el_gen_charge.clear();
   el_gen_dr.clear();

   el_hasMatchedConversion.clear();
   el_expectedMissingInnerHits.clear();

   mu_n = 0;
   mu_pt.clear();
   mu_eta.clear();
   mu_phi.clear();
   mu_m.clear();
   mu_E.clear();
   mu_id.clear();
   mu_charge.clear();

   mu_dB3D.clear();
   mu_edB3D.clear();

   mu_dB.clear();
   mu_edB.clear();
   
   mu_muonBest_dz.clear();
   
   mu_neutralHadronIso.clear();
   mu_chargedHadronIso.clear();
   mu_puChargedHadronIso.clear();
   mu_ecalIso.clear();
   mu_hcalIso.clear();
   mu_photonIso.clear();
   mu_trackIso.clear();

   mu_pfIso03_sumChargedHadronPt.clear();
   mu_pfIso03_sumNeutralHadronEt.clear();
   mu_pfIso03_sumPhotonEt.clear();
   mu_pfIso03_sumPUPt.clear();

   mu_pfIso04_sumChargedHadronPt.clear();
   mu_pfIso04_sumNeutralHadronEt.clear();
   mu_pfIso04_sumPhotonEt.clear();
   mu_pfIso04_sumPUPt.clear();
   
   mu_miniIso.clear();

   mu_isGlobalMuon.clear();
   mu_isTrackerMuon.clear();
   mu_isStandAloneMuon.clear();
   mu_isCaloMuon.clear();
   mu_isPFMuon.clear();

   mu_vx.clear();
   mu_vy.clear();
   mu_vz.clear();
   
   mu_segmentCompatibility.clear();

   mu_isTightMuon.clear();

   mu_hasTrack.clear();
   mu_track_trackerLayersWithMeasurement.clear();
   
   mu_hasGlobalTrack.clear();
   mu_globalTrack_dxy.clear();
   mu_globalTrack_dz.clear();
   mu_globalTrack_dxyError.clear();
   mu_globalTrack_dzError.clear();  
   mu_globalTrack_normalizedChi2.clear();
   
   mu_combinedQuality_chi2LocalPosition.clear();
   mu_combinedQuality_trkKink.clear();

   mu_hasInnerTrack.clear();
   mu_innerTrack_dxy.clear();
   mu_innerTrack_dz.clear();
   mu_innerTrack_dxyError.clear();
   mu_innerTrack_dzError.clear();
   mu_innerTrack_normalizedChi2.clear();
   
   mu_innerTrack_validFraction.clear();

   mu_bestTrackType.clear();
   mu_bestTrack_dxy.clear();
   mu_bestTrack_dz.clear();
   mu_bestTrack_dxyError.clear();
   mu_bestTrack_dzError.clear();
   mu_bestTrack_normalizedChi2.clear();

   mu_innerTrack_pt.clear();
   mu_innerTrack_ptError.clear();
   mu_innerTrack_numberOfValidPixelHits.clear();

   mu_numberOfMatches.clear();
   mu_numberOfMatchedStations.clear();
   mu_numberOfValidMuonHits.clear();

   mu_lepMVA.clear();

   mu_lepMVA_neuRelIso.clear();
   mu_lepMVA_chRelIso.clear();
   mu_lepMVA_jetDR.clear();
   mu_lepMVA_jetPtRatio.clear();
   mu_lepMVA_jetBTagCSV.clear();
   mu_lepMVA_sip3d.clear();
   mu_lepMVA_dxy.clear();
   mu_lepMVA_dz.clear();
   mu_lepMVA_mvaId.clear();

   mu_hasMCMatch.clear();
   mu_gen_pt.clear();
   mu_gen_eta.clear();
   mu_gen_phi.clear();
   mu_gen_m.clear();
   mu_gen_status.clear();
   mu_gen_id.clear();
   mu_gen_charge.clear();
   mu_gen_dr.clear();

   tau_n = 0;
   tau_pt.clear();
   tau_eta.clear();
   tau_phi.clear();
   tau_m.clear();
   tau_E.clear();
   tau_id.clear();
   tau_charge.clear();
   
   tau_hasLeadChargedHadrCand.clear();
   tau_leadingTrackPt.clear();
   tau_leadingTrackCharge.clear();
   
   tau_decayMode.clear();
//   tau_decayModeFindingOldDMs.clear();
   tau_decayModeFindingNewDMs.clear();
   
   tau_puCorrPtSum.clear();
   tau_neutralIsoPtSum.clear();
   tau_chargedIsoPtSum.clear();
   tau_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();
   
   tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
   tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
   tau_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();
   
   tau_againstMuonLoose3.clear();
   tau_againstMuonTight3.clear();
   
   tau_againstElectronVLooseMVA5.clear();
   tau_againstElectronLooseMVA5.clear();
   tau_againstElectronMediumMVA5.clear();
   
   tau_pfEssential_jet_pt.clear();
   tau_pfEssential_jet_eta.clear();
   tau_pfEssential_jet_phi.clear();
   tau_pfEssential_jet_m.clear();

   tau_pfEssential_jetCorr_pt.clear();
   tau_pfEssential_jetCorr_eta.clear();
   tau_pfEssential_jetCorr_phi.clear();
   tau_pfEssential_jetCorr_m.clear();
   
   tau_pfEssential_hasSV.clear();
   tau_pfEssential_sv_x.clear();
   tau_pfEssential_sv_y.clear();
   tau_pfEssential_sv_z.clear();
   
   tau_pfEssential_flightLengthSig.clear();
   tau_pfEssential_dxy.clear();
   tau_pfEssential_dxy_error.clear();
   tau_pfEssential_dxy_Sig.clear();
   
   jet_n = 0;
   jet_pt.clear();
   jet_eta.clear();
   jet_phi.clear();
   jet_m.clear();
   jet_E.clear();

   jet_ntrk.clear();

   jet_JBP.clear();
   jet_JP.clear();
   jet_TCHP.clear();
   jet_TCHE.clear();
   jet_SSVHE.clear();
   jet_SSVHP.clear();
   jet_CMVA.clear();
   jet_CSV.clear();
   jet_CSVv2.clear();
   jet_partonFlavour.clear();
   jet_hadronFlavour.clear();

   jet_neutralHadronEnergy.clear();
   jet_neutralEmEnergy.clear();
   jet_chargedHadronEnergy.clear();
   jet_chargedEmEnergy.clear();
   jet_electronEnergy.clear();
   jet_muonEnergy.clear();
   jet_photonEnergy.clear();

   jet_chargedMultiplicity.clear();
   jet_neutralMultiplicity.clear();
   jet_chargedHadronMultiplicity.clear();
   
   jet_jecFactorUncorrected.clear();
   jet_jecFactorL1FastJet.clear();
   jet_jecFactorL2Relative.clear();
   jet_jecFactorL3Absolute.clear();
   
   jet_Unc.clear();
   
   jet_pileupJetId.clear();

   jet_hasGenJet.clear();   
   jet_genJet_pt.clear();
   jet_genJet_eta.clear();
   jet_genJet_phi.clear();
   jet_genJet_m.clear();
   jet_genJet_E.clear();
   jet_genJet_status.clear();
   jet_genJet_id.clear();

   jet_hasGenParton.clear();   
   jet_genParton_pt.clear();
   jet_genParton_eta.clear();
   jet_genParton_phi.clear();
   jet_genParton_m.clear();
   jet_genParton_E.clear();
   jet_genParton_status.clear();
   jet_genParton_id.clear();
   
   jetPuppi_n = 0;
   jetPuppi_pt.clear();
   jetPuppi_eta.clear();
   jetPuppi_phi.clear();
   jetPuppi_m.clear();
   jetPuppi_E.clear();

   jetPuppi_ntrk.clear();

   jetPuppi_JBP.clear();
   jetPuppi_JP.clear();
   jetPuppi_TCHP.clear();
   jetPuppi_TCHE.clear();
   jetPuppi_SSVHE.clear();
   jetPuppi_SSVHP.clear();
   jetPuppi_CMVA.clear();
   jetPuppi_CSV.clear();
   jetPuppi_CSVv2.clear();
   jetPuppi_partonFlavour.clear();
   jetPuppi_hadronFlavour.clear();

   jetPuppi_neutralHadronEnergy.clear();
   jetPuppi_neutralEmEnergy.clear();
   jetPuppi_chargedHadronEnergy.clear();
   jetPuppi_chargedEmEnergy.clear();
   jetPuppi_electronEnergy.clear();
   jetPuppi_muonEnergy.clear();
   jetPuppi_photonEnergy.clear();

   jetPuppi_chargedMultiplicity.clear();
   jetPuppi_neutralMultiplicity.clear();
   jetPuppi_chargedHadronMultiplicity.clear();
   
   jetPuppi_jecFactorUncorrected.clear();
   jetPuppi_jecFactorL1FastJet.clear();
   jetPuppi_jecFactorL2Relative.clear();
   jetPuppi_jecFactorL3Absolute.clear();
   
   jetPuppi_pileupJetId.clear();

   jetPuppi_hasGenJet.clear();   
   jetPuppi_genJet_pt.clear();
   jetPuppi_genJet_eta.clear();
   jetPuppi_genJet_phi.clear();
   jetPuppi_genJet_m.clear();
   jetPuppi_genJet_E.clear();
   jetPuppi_genJet_status.clear();
   jetPuppi_genJet_id.clear();

   jetPuppi_hasGenParton.clear();   
   jetPuppi_genParton_pt.clear();
   jetPuppi_genParton_eta.clear();
   jetPuppi_genParton_phi.clear();
   jetPuppi_genParton_m.clear();
   jetPuppi_genParton_E.clear();
   jetPuppi_genParton_status.clear();
   jetPuppi_genParton_id.clear();
   
   genJet_n = 0;
   genJet_pt.clear();
   genJet_eta.clear();
   genJet_phi.clear();
   genJet_m.clear();
   genJet_E.clear();
   genJet_emEnergy.clear();
   genJet_hadEnergy.clear();
   genJet_invisibleEnergy.clear();
   genJet_auxiliaryEnergy.clear();
   genJet_flavour.clear();
}

void FlatTree::CreateBranches(int buffersize = 32000)
{
   if( doWrite("ev_run") ) tree->Branch("ev_run", &ev_run, "ev_run/I", buffersize);
   if( doWrite("ev_id") ) tree->Branch("ev_id", &ev_id, "ev_id/I", buffersize);
   if( doWrite("ev_lumi") ) tree->Branch("ev_lumi", &ev_lumi, "ev_lumi/I", buffersize);
   if( doWrite("ev_rho") ) tree->Branch("ev_rho", &ev_rho, "ev_rho/F", buffersize);
   
   if( doWrite("met_pt") ) tree->Branch("met_pt", &met_pt, "met_pt/F", buffersize);
   if( doWrite("met_phi") ) tree->Branch("met_phi", &met_phi, "met_phi/F", buffersize);
   if( doWrite("met_sumet") ) tree->Branch("met_sumet", &met_sumet, "met_sumet/F", buffersize);

   if( doWrite("metPuppi_pt") ) tree->Branch("metPuppi_pt", &metPuppi_pt, "metPuppi_pt/F", buffersize);
   if( doWrite("metPuppi_phi") ) tree->Branch("metPuppi_phi", &metPuppi_phi, "metPuppi_phi/F", buffersize);
   if( doWrite("metPuppi_sumet") ) tree->Branch("metPuppi_sumet", &metPuppi_sumet, "metPuppi_sumet/F", buffersize);
   
   if( doWrite("pv_x") ) tree->Branch("pv_x", &pv_x, "pv_x/F", buffersize);
   if( doWrite("pv_y") ) tree->Branch("pv_y", &pv_y, "pv_y/F", buffersize);
   if( doWrite("pv_z") ) tree->Branch("pv_z", &pv_z, "pv_z/F", buffersize);
   
   if( doWrite("pv_ndof") ) tree->Branch("pv_ndof", &pv_ndof, "pv_ndof/I", buffersize);
   if( doWrite("pv_rho") ) tree->Branch("pv_rho", &pv_rho, "pv_rho/F", buffersize);
   if( doWrite("pv_isFake") ) tree->Branch("pv_isFake", &pv_isFake, "pv_isFake/I", buffersize);

   if( doWrite("mc_weight") ) tree->Branch("mc_weight", &mc_weight, "mc_weight/F", buffersize);
   if( doWrite("mc_id") ) tree->Branch("mc_id", &mc_id, "mc_id/I", buffersize);
   if( doWrite("mc_f1") ) tree->Branch("mc_f1", &mc_f1, "mc_f1/I", buffersize);
   if( doWrite("mc_f2") ) tree->Branch("mc_f2", &mc_f2, "mc_f2/I", buffersize);
   if( doWrite("mc_x1") ) tree->Branch("mc_x1", &mc_x1, "mc_x1/F", buffersize);
   if( doWrite("mc_x2") ) tree->Branch("mc_x2", &mc_x2, "mc_x2/F", buffersize);
   if( doWrite("mc_scale") ) tree->Branch("mc_scale", &mc_scale, "mc_scale/F", buffersize);
   if( doWrite("mc_ptHat") ) tree->Branch("mc_ptHat", &mc_ptHat, "mc_ptHat/F", buffersize);
   
   if( doWrite("mc_pu_intime_NumInt") ) tree->Branch("mc_pu_intime_NumInt", &mc_pu_intime_NumInt, "mc_pu_intime_NumInt/I", buffersize);
   if( doWrite("mc_pu_trueNumInt") ) tree->Branch("mc_pu_trueNumInt", &mc_pu_trueNumInt, "mc_pu_trueNumInt/I", buffersize);
   if( doWrite("mc_pu_before_npu") ) tree->Branch("mc_pu_before_npu", &mc_pu_before_npu, "mc_pu_before_npu/I", buffersize);
   if( doWrite("mc_pu_after_npu") ) tree->Branch("mc_pu_after_npu", &mc_pu_after_npu, "mc_pu_after_npu/I", buffersize);
   
   if( doWrite("mc_pu_Npvi") ) tree->Branch("mc_pu_Npvi", &mc_pu_Npvi, "mc_pu_Npvi/I", buffersize);
   if( doWrite("mc_pu_Nzpositions") ) tree->Branch("mc_pu_Nzpositions", "std::vector<int>", &mc_pu_Nzpositions, buffersize);
   if( doWrite("mc_pu_BunchCrossing") ) tree->Branch("mc_pu_BunchCrossing", "std::vector<int>", &mc_pu_BunchCrossing, buffersize );
   if( doWrite("mc_pu_zpositions") ) tree->Branch("mc_pu_zpositions", "std::vector<std::vector<float> >", &mc_pu_zpositions, buffersize);
   if( doWrite("mc_pu_sumpT_lowpT") ) tree->Branch("mc_pu_sumpT_lowpT", "std::vector<std::vector<float> >", &mc_pu_sumpT_lowpT, buffersize);
   if( doWrite("mc_pu_sumpT_highpT") ) tree->Branch("mc_pu_sumpT_highpT", "std::vector<std::vector<float> >", &mc_pu_sumpT_highpT, buffersize);
   if( doWrite("mc_pu_ntrks_lowpT") ) tree->Branch("mc_pu_ntrks_lowpT", "std::vector<std::vector<int> >", &mc_pu_ntrks_lowpT, buffersize);
   if( doWrite("mc_pu_ntrks_highpT") ) tree->Branch("mc_pu_ntrks_highpT", "std::vector<std::vector<int> >", &mc_pu_ntrks_highpT, buffersize);
 
   if( doWrite("trigger_n") ) tree->Branch("trigger_n", &trigger_n, "trigger_n/I", buffersize); 
   if( doWrite("trigger") ) tree->Branch("trigger", "std::vector<int>", &trigger, buffersize);
   if( doWrite("trigger_name") ) tree->Branch("trigger_name", "std::vector<string>", &trigger_name, buffersize);
   if( doWrite("trigger_pass") ) tree->Branch("trigger_pass", "std::vector<bool>", &trigger_pass, buffersize);
   if( doWrite("trigger_prescale") ) tree->Branch("trigger_prescale", "std::vector<int>", &trigger_prescale, buffersize);
   
   if( doWrite("triggerobject_n") ) tree->Branch("triggerobject_n", &triggerobject_n, "triggerobject_n/I",buffersize);
   if( doWrite("triggerobject_pt") ) tree->Branch("triggerobject_pt", "std::vector<float>", &triggerobject_pt, buffersize);
   if( doWrite("triggerobject_eta") ) tree->Branch("triggerobject_eta", "std::vector<float>", &triggerobject_eta, buffersize);
   if( doWrite("triggerobject_phi") ) tree->Branch("triggerobject_phi", "std::vector<float>", &triggerobject_phi, buffersize);

   if( doWrite("triggerobject_collection") ) tree->Branch("triggerobject_collection", "std::vector<std::string>", &triggerobject_collection, buffersize);

   if( doWrite("triggerobject_filterIds_n") ) tree->Branch("triggerobject_filterIds_n", "std::vector<int>", &triggerobject_filterIds_n, buffersize);
   if( doWrite("triggerobject_filterIds") ) tree->Branch("triggerobject_filterIds", "std::vector<int>", &triggerobject_filterIds, buffersize);

   if( doWrite("triggerobject_isTriggerL1Mu") ) tree->Branch("triggerobject_isTriggerL1Mu", "std::vector<bool>", &triggerobject_isTriggerL1Mu, buffersize);
   if( doWrite("triggerobject_isTriggerL1NoIsoEG") ) tree->Branch("triggerobject_isTriggerL1NoIsoEG", "std::vector<bool>", &triggerobject_isTriggerL1NoIsoEG, buffersize);
   if( doWrite("triggerobject_isTriggerL1IsoEG") ) tree->Branch("triggerobject_isTriggerL1IsoEG", "std::vector<bool>", &triggerobject_isTriggerL1IsoEG, buffersize);
   if( doWrite("triggerobject_isTriggerL1CenJet") ) tree->Branch("triggerobject_isTriggerL1CenJet", "std::vector<bool>", &triggerobject_isTriggerL1CenJet, buffersize);
   if( doWrite("triggerobject_isTriggerL1ForJet") ) tree->Branch("triggerobject_isTriggerL1ForJet", "std::vector<bool>", &triggerobject_isTriggerL1ForJet, buffersize);
   if( doWrite("triggerobject_isTriggerL1TauJet") ) tree->Branch("triggerobject_isTriggerL1TauJet", "std::vector<bool>", &triggerobject_isTriggerL1TauJet, buffersize);
   if( doWrite("triggerobject_isTriggerL1ETM") ) tree->Branch("triggerobject_isTriggerL1ETM", "std::vector<bool>", &triggerobject_isTriggerL1ETM, buffersize);
   if( doWrite("triggerobject_isTriggerL1ETT") ) tree->Branch("triggerobject_isTriggerL1ETT", "std::vector<bool>", &triggerobject_isTriggerL1ETT, buffersize);
   if( doWrite("triggerobject_isTriggerL1HTT") ) tree->Branch("triggerobject_isTriggerL1HTT", "std::vector<bool>", &triggerobject_isTriggerL1HTT, buffersize);
   if( doWrite("triggerobject_isTriggerL1HTM") ) tree->Branch("triggerobject_isTriggerL1HTM", "std::vector<bool>", &triggerobject_isTriggerL1HTM, buffersize);
   if( doWrite("triggerobject_isTriggerL1JetCounts") ) tree->Branch("triggerobject_isTriggerL1JetCounts", "std::vector<bool>", &triggerobject_isTriggerL1JetCounts, buffersize);
   if( doWrite("triggerobject_isTriggerL1HfBitCounts") ) tree->Branch("triggerobject_isTriggerL1HfBitCounts", "std::vector<bool>", &triggerobject_isTriggerL1HfBitCounts, buffersize);
   if( doWrite("triggerobject_isTriggerL1HfRingEtSums") ) tree->Branch("triggerobject_isTriggerL1HfRingEtSums", "std::vector<bool>", &triggerobject_isTriggerL1HfRingEtSums, buffersize);
   if( doWrite("triggerobject_isTriggerL1TechTrig") ) tree->Branch("triggerobject_isTriggerL1TechTrig", "std::vector<bool>", &triggerobject_isTriggerL1TechTrig, buffersize);
   if( doWrite("triggerobject_isTriggerL1Castor") ) tree->Branch("triggerobject_isTriggerL1Castor", "std::vector<bool>", &triggerobject_isTriggerL1Castor, buffersize);
   if( doWrite("triggerobject_isTriggerL1BPTX") ) tree->Branch("triggerobject_isTriggerL1BPTX", "std::vector<bool>", &triggerobject_isTriggerL1BPTX, buffersize);
   if( doWrite("triggerobject_isTriggerL1GtExternal") ) tree->Branch("triggerobject_isTriggerL1GtExternal", "std::vector<bool>", &triggerobject_isTriggerL1GtExternal, buffersize);

   if( doWrite("triggerobject_isHLT_TriggerPhoton") ) tree->Branch("triggerobject_isHLT_TriggerPhoton", "std::vector<bool>", &triggerobject_isHLT_TriggerPhoton, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerElectron") ) tree->Branch("triggerobject_isHLT_TriggerElectron", "std::vector<bool>", &triggerobject_isHLT_TriggerElectron, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerMuon") ) tree->Branch("triggerobject_isHLT_TriggerMuon", "std::vector<bool>", &triggerobject_isHLT_TriggerMuon, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerTau") ) tree->Branch("triggerobject_isHLT_TriggerTau", "std::vector<bool>", &triggerobject_isHLT_TriggerTau, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerJet") ) tree->Branch("triggerobject_isHLT_TriggerJet", "std::vector<bool>", &triggerobject_isHLT_TriggerJet, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerBJet") ) tree->Branch("triggerobject_isHLT_TriggerBJet", "std::vector<bool>", &triggerobject_isHLT_TriggerBJet, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerMET") ) tree->Branch("triggerobject_isHLT_TriggerMET", "std::vector<bool>", &triggerobject_isHLT_TriggerMET, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerTET") ) tree->Branch("triggerobject_isHLT_TriggerTET", "std::vector<bool>", &triggerobject_isHLT_TriggerTET, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerTHT") ) tree->Branch("triggerobject_isHLT_TriggerTHT", "std::vector<bool>", &triggerobject_isHLT_TriggerTHT, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerMHT") ) tree->Branch("triggerobject_isHLT_TriggerMHT", "std::vector<bool>", &triggerobject_isHLT_TriggerMHT, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerTrack") ) tree->Branch("triggerobject_isHLT_TriggerTrack", "std::vector<bool>", &triggerobject_isHLT_TriggerTrack, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerCluster") ) tree->Branch("triggerobject_isHLT_TriggerCluster", "std::vector<bool>", &triggerobject_isHLT_TriggerCluster, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerMETSig") ) tree->Branch("triggerobject_isHLT_TriggerMETSig", "std::vector<bool>", &triggerobject_isHLT_TriggerMETSig, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerELongit") ) tree->Branch("triggerobject_isHLT_TriggerELongit", "std::vector<bool>", &triggerobject_isHLT_TriggerELongit, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerMHTSig") ) tree->Branch("triggerobject_isHLT_TriggerMHTSig", "std::vector<bool>", &triggerobject_isHLT_TriggerMHTSig, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerHLongit") ) tree->Branch("triggerobject_isHLT_TriggerHLongit", "std::vector<bool>", &triggerobject_isHLT_TriggerHLongit, buffersize);

   if( doWrite("triggerobject_filterLabels_n") ) tree->Branch("triggerobject_filterLabels_n", "std::vector<int>", &triggerobject_filterLabels_n, buffersize);
   if( doWrite("triggerobject_filterLabels") ) tree->Branch("triggerobject_filterLabels", "std::vector<string>", &triggerobject_filterLabels, buffersize);

   if( doWrite("triggerobject_pathNamesAll_n") ) tree->Branch("triggerobject_pathNamesAll_n", "std::vector<int>", &triggerobject_pathNamesAll_n, buffersize);
   if( doWrite("triggerobject_pathNamesAll") ) tree->Branch("triggerobject_pathNamesAll", "std::vector<string>", &triggerobject_pathNamesAll, buffersize);
   if( doWrite("triggerobject_pathNamesAll_isBoth") ) tree->Branch("triggerobject_pathNamesAll_isBoth", "std::vector<bool>", &triggerobject_pathNamesAll_isBoth, buffersize);
   if( doWrite("triggerobject_pathNamesAll_isL3") ) tree->Branch("triggerobject_pathNamesAll_isL3", "std::vector<bool>", &triggerobject_pathNamesAll_isL3, buffersize);
   if( doWrite("triggerobject_pathNamesAll_isLF") ) tree->Branch("triggerobject_pathNamesAll_isLF", "std::vector<bool>", &triggerobject_pathNamesAll_isLF, buffersize);
   if( doWrite("triggerobject_pathNamesAll_isNone") ) tree->Branch("triggerobject_pathNamesAll_isNone", "std::vector<bool>", &triggerobject_pathNamesAll_isNone, buffersize);

   if( doWrite("el_n") ) tree->Branch("el_n", &el_n, "el_n/I", buffersize);
   if( doWrite("el_pt") ) tree->Branch("el_pt", "std::vector<float>", &el_pt, buffersize);
   if( doWrite("el_eta") ) tree->Branch("el_eta", "std::vector<float>", &el_eta, buffersize);
   if( doWrite("el_phi") ) tree->Branch("el_phi", "std::vector<float>", &el_phi, buffersize);
   if( doWrite("el_m") ) tree->Branch("el_m", "std::vector<float>", &el_m, buffersize);
   if( doWrite("el_E") ) tree->Branch("el_E", "std::vector<float>", &el_E, buffersize);
   if( doWrite("el_id") ) tree->Branch("el_id", "std::vector<int>", &el_id, buffersize);
   if( doWrite("el_charge") ) tree->Branch("el_charge", "std::vector<int>", &el_charge, buffersize);
   
   if( doWrite("el_scleta") ) tree->Branch("el_scleta", "std::vector<float>", &el_scleta, buffersize);
   if( doWrite("el_passConversionVeto") ) tree->Branch("el_passConversionVeto", "std::vector<int>", &el_passConversionVeto, buffersize);   

   if( doWrite("el_isGsfCtfScPixChargeConsistent") ) tree->Branch("el_isGsfCtfScPixChargeConsistent", "std::vector<int>", &el_isGsfCtfScPixChargeConsistent, buffersize);
   if( doWrite("el_dB3D") ) tree->Branch("el_dB3D", "std::vector<float>", &el_dB3D, buffersize);
   if( doWrite("el_edB3D") ) tree->Branch("el_edB3D", "std::vector<float>", &el_edB3D, buffersize);
   
   if( doWrite("el_neutralHadronIso") ) tree->Branch("el_neutralHadronIso", "std::vector<float>", &el_neutralHadronIso, buffersize);
   if( doWrite("el_chargedHadronIso") ) tree->Branch("el_chargedHadronIso", "std::vector<float>", &el_chargedHadronIso, buffersize);
   if( doWrite("el_puChargedHadronIso") ) tree->Branch("el_puChargedHadronIso", "std::vector<float>", &el_puChargedHadronIso, buffersize);
   if( doWrite("el_ecalIso") ) tree->Branch("el_ecalIso", "std::vector<float>", &el_ecalIso, buffersize);
   if( doWrite("el_hcalIso") ) tree->Branch("el_hcalIso", "std::vector<float>", &el_hcalIso, buffersize);
   if( doWrite("el_particleIso") ) tree->Branch("el_particleIso", "std::vector<float>", &el_particleIso, buffersize);
   if( doWrite("el_photonIso") ) tree->Branch("el_photonIso", "std::vector<float>", &el_photonIso, buffersize);
   if( doWrite("el_trackIso") ) tree->Branch("el_trackIso", "std::vector<float>", &el_trackIso, buffersize);
   
   if( doWrite("el_pfIso_sumChargedHadronPt") ) tree->Branch("el_pfIso_sumChargedHadronPt", "std::vector<float>", &el_pfIso_sumChargedHadronPt, buffersize);
   if( doWrite("el_pfIso_sumNeutralHadronEt") ) tree->Branch("el_pfIso_sumNeutralHadronEt", "std::vector<float>", &el_pfIso_sumNeutralHadronEt, buffersize);
   if( doWrite("el_pfIso_sumPhotonEt") ) tree->Branch("el_pfIso_sumPhotonEt", "std::vector<float>", &el_pfIso_sumPhotonEt, buffersize);
   if( doWrite("el_pfIso_sumPUPt") ) tree->Branch("el_pfIso_sumPUPt", "std::vector<float>", &el_pfIso_sumPUPt, buffersize);
   
   if( doWrite("el_miniIso") ) tree->Branch("el_miniIso", "std::vector<float>", &el_miniIso, buffersize);
   
   if( doWrite("el_isLoose") ) tree->Branch("el_isLoose", "std::vector<int>", &el_isLoose, buffersize);
   if( doWrite("el_isTight") ) tree->Branch("el_isTight", "std::vector<int>", &el_isTight, buffersize);
   if( doWrite("el_isRobustLoose") ) tree->Branch("el_isRobustLoose", "std::vector<int>", &el_isRobustLoose, buffersize);
   if( doWrite("el_isRobustTight") ) tree->Branch("el_isRobustTight", "std::vector<int>", &el_isRobustTight, buffersize);
   if( doWrite("el_isRobustHighEnergy") ) tree->Branch("el_isRobustHighEnergy", "std::vector<int>", &el_isRobustHighEnergy, buffersize);
   
   if( doWrite("el_vx") ) tree->Branch("el_vx", "std::vector<float>", &el_vx, buffersize);
   if( doWrite("el_vy") ) tree->Branch("el_vy", "std::vector<float>", &el_vy, buffersize);
   if( doWrite("el_vz") ) tree->Branch("el_vz", "std::vector<float>", &el_vz, buffersize);
   
   if( doWrite("el_isGsf") ) tree->Branch("el_isGsf", "std::vector<bool>", &el_isGsf, buffersize);
   if( doWrite("el_dxy") ) tree->Branch("el_dxy", "std::vector<float>", &el_dxy, buffersize);
   if( doWrite("el_dz") ) tree->Branch("el_dz", "std::vector<float>", &el_dz, buffersize);
   if( doWrite("el_dxyError") ) tree->Branch("el_dxyError", "std::vector<float>", &el_dxyError, buffersize);
   if( doWrite("el_dzError") ) tree->Branch("el_dzError", "std::vector<float>", &el_dzError, buffersize);
   
   if( doWrite("el_numberOfHits") ) tree->Branch("el_numberOfHits", "std::vector<int>", &el_numberOfHits, buffersize);
   
   if( doWrite("el_sigmaIetaIeta") ) tree->Branch("el_sigmaIetaIeta", "std::vector<float>", &el_sigmaIetaIeta, buffersize);
   if( doWrite("el_sigmaIphiIphi") ) tree->Branch("el_sigmaIphiIphi", "std::vector<float>", &el_sigmaIphiIphi, buffersize);
   if( doWrite("el_hadronicOverEm") ) tree->Branch("el_hadronicOverEm", "std::vector<float>", &el_hadronicOverEm, buffersize);
   if( doWrite("el_dr03TkSumPt") ) tree->Branch("el_dr03TkSumPt", "std::vector<float>", &el_dr03TkSumPt, buffersize);
   if( doWrite("el_dr03EcalRecHitSumEt") ) tree->Branch("el_dr03EcalRecHitSumEt", "std::vector<float>", &el_dr03EcalRecHitSumEt, buffersize);
   if( doWrite("el_dr03HcalTowerSumEt") ) tree->Branch("el_dr03HcalTowerSumEt", "std::vector<float>", &el_dr03HcalTowerSumEt, buffersize);
   if( doWrite("el_numberOfLostHits") ) tree->Branch("el_numberOfLostHits", "std::vector<int>", &el_numberOfLostHits, buffersize);

   if( doWrite("el_fbrem") ) tree->Branch("el_fbrem", "std::vector<float>", &el_fbrem, buffersize);
   if( doWrite("el_kf_normalizedChi2") ) tree->Branch("el_kf_normalizedChi2", "std::vector<float>", &el_kf_normalizedChi2, buffersize);
   if( doWrite("el_trackerLayersWithMeasurement") ) tree->Branch("el_trackerLayersWithMeasurement", "std::vector<int>", &el_trackerLayersWithMeasurement, buffersize);
   if( doWrite("el_gsf_normalizedChi2") ) tree->Branch("el_gsf_normalizedChi2", "std::vector<float>", &el_gsf_normalizedChi2, buffersize);
   if( doWrite("el_deltaEtaSuperClusterTrackAtVtx") ) tree->Branch("el_deltaEtaSuperClusterTrackAtVtx", "std::vector<float>", &el_deltaEtaSuperClusterTrackAtVtx, buffersize);
   if( doWrite("el_deltaPhiSuperClusterTrackAtVtx") ) tree->Branch("el_deltaPhiSuperClusterTrackAtVtx", "std::vector<float>", &el_deltaPhiSuperClusterTrackAtVtx, buffersize);
   if( doWrite("el_deltaEtaSeedClusterTrackAtCalo") ) tree->Branch("el_deltaEtaSeedClusterTrackAtCalo", "std::vector<float>", &el_deltaEtaSeedClusterTrackAtCalo, buffersize);
   if( doWrite("el_see") ) tree->Branch("el_see", "std::vector<float>", &el_see, buffersize);
   if( doWrite("el_spp") ) tree->Branch("el_spp", "std::vector<float>", &el_spp, buffersize);
   if( doWrite("el_superClusterEtaWidth") ) tree->Branch("el_superClusterEtaWidth", "std::vector<float>", &el_superClusterEtaWidth, buffersize);
   if( doWrite("el_superClusterPhiWidth") ) tree->Branch("el_superClusterPhiWidth", "std::vector<float>", &el_superClusterPhiWidth, buffersize);
   if( doWrite("el_full5x5_OneMinusE1x5E5x5") ) tree->Branch("el_full5x5_OneMinusE1x5E5x5", "std::vector<float>", &el_full5x5_OneMinusE1x5E5x5, buffersize);
   if( doWrite("el_OneMinusE1x5E5x5") ) tree->Branch("el_OneMinusE1x5E5x5", "std::vector<float>", &el_OneMinusE1x5E5x5, buffersize);
   if( doWrite("el_full5x5_r9") ) tree->Branch("el_full5x5_r9", "std::vector<float>", &el_full5x5_r9, buffersize);
   if( doWrite("el_r9") ) tree->Branch("el_r9", "std::vector<float>", &el_r9, buffersize);
   if( doWrite("el_eSuperClusterOverP") ) tree->Branch("el_eSuperClusterOverP", "std::vector<float>", &el_eSuperClusterOverP, buffersize);
   if( doWrite("el_IoEmIoP") ) tree->Branch("el_IoEmIoP", "std::vector<float>", &el_IoEmIoP, buffersize);
   if( doWrite("el_ooEmooP") ) tree->Branch("el_ooEmooP", "std::vector<float>", &el_ooEmooP, buffersize);
   if( doWrite("el_eleEoPout") ) tree->Branch("el_eleEoPout", "std::vector<float>", &el_eleEoPout, buffersize);
   if( doWrite("el_PreShowerOverRaw") ) tree->Branch("el_PreShowerOverRaw", "std::vector<float>", &el_PreShowerOverRaw, buffersize);
   if( doWrite("el_ecalEnergy") ) tree->Branch("el_ecalEnergy", "std::vector<float>", &el_ecalEnergy, buffersize);

   if( doWrite("el_mvaNonTrigV0") ) tree->Branch("el_mvaNonTrigV0", "std::vector<float>", &el_mvaNonTrigV0, buffersize);

   if( doWrite("el_lepMVA") ) tree->Branch("el_lepMVA", "std::vector<float>", &el_lepMVA, buffersize);

   if( doWrite("el_lepMVA_neuRelIso") ) tree->Branch("el_lepMVA_neuRelIso", "std::vector<float>", &el_lepMVA_neuRelIso, buffersize);
   if( doWrite("el_lepMVA_chRelIso") ) tree->Branch("el_lepMVA_chRelIso", "std::vector<float>", &el_lepMVA_chRelIso, buffersize);
   if( doWrite("el_lepMVA_jetDR") ) tree->Branch("el_lepMVA_jetDR", "std::vector<float>", &el_lepMVA_jetDR, buffersize);
   if( doWrite("el_lepMVA_jetPtRatio") ) tree->Branch("el_lepMVA_jetPtRatio", "std::vector<float>", &el_lepMVA_jetPtRatio, buffersize);
   if( doWrite("el_lepMVA_jetBTagCSV") ) tree->Branch("el_lepMVA_jetBTagCSV", "std::vector<float>", &el_lepMVA_jetBTagCSV, buffersize);
   if( doWrite("el_lepMVA_sip3d") ) tree->Branch("el_lepMVA_sip3d", "std::vector<float>", &el_lepMVA_sip3d, buffersize);
   if( doWrite("el_lepMVA_dxy") ) tree->Branch("el_lepMVA_dxy", "std::vector<float>", &el_lepMVA_dxy, buffersize);
   if( doWrite("el_lepMVA_dz") ) tree->Branch("el_lepMVA_dz", "std::vector<float>", &el_lepMVA_dz, buffersize);
   if( doWrite("el_lepMVA_mvaId") ) tree->Branch("el_lepMVA_mvaId", "std::vector<float>", &el_lepMVA_mvaId, buffersize);

   if( doWrite("el_hasMCMatch") ) tree->Branch("el_hasMCMatch", "std::vector<int>", &el_hasMCMatch, buffersize);
   if( doWrite("el_gen_pt") ) tree->Branch("el_gen_pt", "std::vector<float>", &el_gen_pt, buffersize);
   if( doWrite("el_gen_eta") ) tree->Branch("el_gen_eta", "std::vector<float>", &el_gen_eta, buffersize);
   if( doWrite("el_gen_phi") ) tree->Branch("el_gen_phi", "std::vector<float>", &el_gen_phi, buffersize);
   if( doWrite("el_gen_m") ) tree->Branch("el_gen_m", "std::vector<float>", &el_gen_m, buffersize);
   if( doWrite("el_gen_status") ) tree->Branch("el_gen_status", "std::vector<int>", &el_gen_status, buffersize);
   if( doWrite("el_gen_id") ) tree->Branch("el_gen_id", "std::vector<int>", &el_gen_id, buffersize);
   if( doWrite("el_gen_charge") ) tree->Branch("el_gen_charge", "std::vector<int>", &el_gen_charge, buffersize);
   if( doWrite("el_gen_dr") ) tree->Branch("el_gen_dr", "std::vector<float>", &el_gen_dr, buffersize);

   if( doWrite("el_hasMatchedConversion") ) tree->Branch("el_hasMatchedConversion", "std::vector<bool>", &el_hasMatchedConversion, buffersize);
   if( doWrite("el_expectedMissingInnerHits") ) tree->Branch("el_expectedMissingInnerHits", "std::vector<int>", &el_expectedMissingInnerHits, buffersize);

   if( doWrite("mu_n") ) tree->Branch("mu_n", &mu_n, "mu_n/I", buffersize);
   if( doWrite("mu_pt") ) tree->Branch("mu_pt", "std::vector<float>", &mu_pt, buffersize);
   if( doWrite("mu_eta") ) tree->Branch("mu_eta", "std::vector<float>", &mu_eta, buffersize);
   if( doWrite("mu_phi") ) tree->Branch("mu_phi", "std::vector<float>", &mu_phi, buffersize);
   if( doWrite("mu_m") ) tree->Branch("mu_m", "std::vector<float>", &mu_m, buffersize);
   if( doWrite("mu_E") ) tree->Branch("mu_E", "std::vector<float>", &mu_E, buffersize);
   if( doWrite("mu_id") ) tree->Branch("mu_id", "std::vector<int>", &mu_id, buffersize);
   if( doWrite("mu_charge") ) tree->Branch("mu_charge", "std::vector<int>", &mu_charge, buffersize);

   if( doWrite("mu_dB3D") ) tree->Branch("mu_dB3D", "std::vector<float>", &mu_dB3D, buffersize);
   if( doWrite("mu_edB3D") ) tree->Branch("mu_edB3D", "std::vector<float>", &mu_edB3D, buffersize);

   if( doWrite("mu_dB") ) tree->Branch("mu_dB", "std::vector<float>", &mu_dB, buffersize);
   if( doWrite("mu_edB") ) tree->Branch("mu_edB", "std::vector<float>", &mu_edB, buffersize);
   
   if( doWrite("mu_muonBest_dz") ) tree->Branch("mu_muonBest_dz", "std::vector<float>", &mu_muonBest_dz, buffersize);
   
   if( doWrite("mu_neutralHadronIso") ) tree->Branch("mu_neutralHadronIso", "std::vector<float>", &mu_neutralHadronIso, buffersize);
   if( doWrite("mu_chargedHadronIso") ) tree->Branch("mu_chargedHadronIso", "std::vector<float>", &mu_chargedHadronIso, buffersize);
   if( doWrite("mu_puChargedHadronIso") ) tree->Branch("mu_puChargedHadronIso", "std::vector<float>", &mu_puChargedHadronIso, buffersize);
   if( doWrite("mu_ecalIso") ) tree->Branch("mu_ecalIso", "std::vector<float>", &mu_ecalIso, buffersize);
   if( doWrite("mu_hcalIso") ) tree->Branch("mu_hcalIso", "std::vector<float>", &mu_hcalIso, buffersize);
   if( doWrite("mu_photonIso") ) tree->Branch("mu_photonIso", "std::vector<float>", &mu_photonIso, buffersize);
   if( doWrite("mu_trackIso") ) tree->Branch("mu_trackIso", "std::vector<float>", &mu_trackIso, buffersize);

   if( doWrite("mu_pfIso03_sumChargedHadronPt") ) tree->Branch("mu_pfIso03_sumChargedHadronPt", "std::vector<float>", &mu_pfIso03_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfIso03_sumNeutralHadronEt") ) tree->Branch("mu_pfIso03_sumNeutralHadronEt", "std::vector<float>", &mu_pfIso03_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfIso03_sumPhotonEt") ) tree->Branch("mu_pfIso03_sumPhotonEt", "std::vector<float>", &mu_pfIso03_sumPhotonEt, buffersize);
   if( doWrite("mu_pfIso03_sumPUPt") ) tree->Branch("mu_pfIso03_sumPUPt", "std::vector<float>", &mu_pfIso03_sumPUPt, buffersize);

   if( doWrite("mu_pfIso04_sumChargedHadronPt") ) tree->Branch("mu_pfIso04_sumChargedHadronPt", "std::vector<float>", &mu_pfIso04_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfIso04_sumNeutralHadronEt") ) tree->Branch("mu_pfIso04_sumNeutralHadronEt", "std::vector<float>", &mu_pfIso04_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfIso04_sumPhotonEt") ) tree->Branch("mu_pfIso04_sumPhotonEt", "std::vector<float>", &mu_pfIso04_sumPhotonEt, buffersize);
   if( doWrite("mu_pfIso04_sumPUPt") ) tree->Branch("mu_pfIso04_sumPUPt", "std::vector<float>", &mu_pfIso04_sumPUPt, buffersize);
   
   if( doWrite("mu_miniIso") ) tree->Branch("mu_miniIso", "std::vector<float>", &mu_miniIso, buffersize);

   if( doWrite("mu_isGlobalMuon") ) tree->Branch("mu_isGlobalMuon", "std::vector<int>", &mu_isGlobalMuon, buffersize);
   if( doWrite("mu_isTrackerMuon") ) tree->Branch("mu_isTrackerMuon", "std::vector<int>", &mu_isTrackerMuon, buffersize);
   if( doWrite("mu_isStandAloneMuon") ) tree->Branch("mu_isStandAloneMuon", "std::vector<int>", &mu_isStandAloneMuon, buffersize);
   if( doWrite("mu_isCaloMuon") ) tree->Branch("mu_isCaloMuon", "std::vector<int>", &mu_isCaloMuon, buffersize);
   if( doWrite("mu_isPFMuon") ) tree->Branch("mu_isPFMuon", "std::vector<int>", &mu_isPFMuon, buffersize);

   if( doWrite("mu_vx") ) tree->Branch("mu_vx", "std::vector<float>", &mu_vx, buffersize);
   if( doWrite("mu_vy") ) tree->Branch("mu_vy", "std::vector<float>", &mu_vy, buffersize);
   if( doWrite("mu_vz") ) tree->Branch("mu_vz", "std::vector<float>", &mu_vz, buffersize);
   
   if( doWrite("mu_segmentCompatibility") ) tree->Branch("mu_segmentCompatibility", "std::vector<float>", &mu_segmentCompatibility, buffersize);

   if( doWrite("mu_isTightMuon") ) tree->Branch("mu_isTightMuon", "std::vector<bool>", &mu_isTightMuon, buffersize);

   if( doWrite("mu_hasTrack") ) tree->Branch("mu_hasTrack", "std::vector<int>", &mu_hasTrack, buffersize);
   if( doWrite("mu_track_trackerLayersWithMeasurement") ) tree->Branch("mu_track_trackerLayersWithMeasurement", "std::vector<int>", &mu_track_trackerLayersWithMeasurement, buffersize);
   
   if( doWrite("mu_hasGlobalTrack") ) tree->Branch("mu_hasGlobalTrack", "std::vector<int>", &mu_hasGlobalTrack, buffersize);
   if( doWrite("mu_globalTrack_dxy") ) tree->Branch("mu_globalTrack_dxy", "std::vector<float>", &mu_globalTrack_dxy, buffersize);
   if( doWrite("mu_globalTrack_dz") ) tree->Branch("mu_globalTrack_dz", "std::vector<float>", &mu_globalTrack_dz, buffersize);
   if( doWrite("mu_globalTrack_dxyError") ) tree->Branch("mu_globalTrack_dxyError", "std::vector<float>", &mu_globalTrack_dxyError, buffersize);
   if( doWrite("mu_globalTrack_dzError") ) tree->Branch("mu_globalTrack_dzError", "std::vector<float>", &mu_globalTrack_dzError, buffersize);   
   if( doWrite("mu_globalTrack_normalizedChi2") ) tree->Branch("mu_globalTrack_normalizedChi2", "std::vector<float>", &mu_globalTrack_normalizedChi2, buffersize);
   
   if( doWrite("mu_combinedQuality_chi2LocalPosition") ) tree->Branch("mu_combinedQuality_chi2LocalPosition", "std::vector<float>", &mu_combinedQuality_chi2LocalPosition, buffersize);
   if( doWrite("mu_combinedQuality_trkKink") ) tree->Branch("mu_combinedQuality_trkKink", "std::vector<float>", &mu_combinedQuality_trkKink, buffersize);

   if( doWrite("mu_hasInnerTrack") ) tree->Branch("mu_hasInnerTrack", "std::vector<int>", &mu_hasInnerTrack, buffersize);
   if( doWrite("mu_innerTrack_dxy") ) tree->Branch("mu_innerTrack_dxy", "std::vector<float>", &mu_innerTrack_dxy, buffersize);
   if( doWrite("mu_innerTrack_dz") ) tree->Branch("mu_innerTrack_dz", "std::vector<float>", &mu_innerTrack_dz, buffersize);
   if( doWrite("mu_innerTrack_dxyError") ) tree->Branch("mu_innerTrack_dxyError", "std::vector<float>", &mu_innerTrack_dxyError, buffersize);
   if( doWrite("mu_innerTrack_dzError") ) tree->Branch("mu_innerTrack_dzError", "std::vector<float>", &mu_innerTrack_dzError, buffersize);
   if( doWrite("mu_innerTrack_normalizedChi2") ) tree->Branch("mu_innerTrack_normalizedChi2", "std::vector<float>", &mu_innerTrack_normalizedChi2, buffersize);
   
   if( doWrite("mu_innerTrack_validFraction") ) tree->Branch("mu_innerTrack_validFraction", "std::vector<float>", &mu_innerTrack_validFraction, buffersize);

   if( doWrite("mu_bestTrackType") ) tree->Branch("mu_bestTrackType", "std::vector<int>", &mu_bestTrackType, buffersize);
   if( doWrite("mu_bestTrack_dxy") ) tree->Branch("mu_bestTrack_dxy", "std::vector<float>", &mu_bestTrack_dxy, buffersize);
   if( doWrite("mu_bestTrack_dz") ) tree->Branch("mu_bestTrack_dz", "std::vector<float>", &mu_bestTrack_dz, buffersize);
   if( doWrite("mu_bestTrack_dxyError") ) tree->Branch("mu_bestTrack_dxyError", "std::vector<float>", &mu_bestTrack_dxyError, buffersize);
   if( doWrite("mu_bestTrack_dzError") ) tree->Branch("mu_bestTrack_dzError", "std::vector<float>", &mu_bestTrack_dzError, buffersize);
   if( doWrite("mu_bestTrack_normalizedChi2") ) tree->Branch("mu_bestTrack_normalizedChi2", "std::vector<float>", &mu_bestTrack_normalizedChi2, buffersize);

   if( doWrite("mu_innerTrack_pt") ) tree->Branch("mu_innerTrack_pt", "std::vector<float>", &mu_innerTrack_pt, buffersize);
   if( doWrite("mu_innerTrack_ptError") ) tree->Branch("mu_innerTrack_ptError", "std::vector<float>", &mu_innerTrack_ptError, buffersize);
   if( doWrite("mu_innerTrack_numberOfValidPixelHits") ) tree->Branch("mu_innerTrack_numberOfValidPixelHits", "std::vector<int>", &mu_innerTrack_numberOfValidPixelHits, buffersize);

   if( doWrite("mu_numberOfMatches") ) tree->Branch("mu_numberOfMatches", "std::vector<int>", &mu_numberOfMatches, buffersize);
   if( doWrite("mu_numberOfMatchedStations") ) tree->Branch("mu_numberOfMatchedStations", "std::vector<int>", &mu_numberOfMatchedStations, buffersize);
   if( doWrite("mu_numberOfValidMuonHits") ) tree->Branch("mu_numberOfValidMuonHits", "std::vector<int>", &mu_numberOfValidMuonHits, buffersize);

   if( doWrite("mu_lepMVA") ) tree->Branch("mu_lepMVA", "std::vector<float>", &mu_lepMVA, buffersize);

   if( doWrite("mu_lepMVA_neuRelIso") ) tree->Branch("mu_lepMVA_neuRelIso", "std::vector<float>", &mu_lepMVA_neuRelIso, buffersize);
   if( doWrite("mu_lepMVA_chRelIso") ) tree->Branch("mu_lepMVA_chRelIso", "std::vector<float>", &mu_lepMVA_chRelIso, buffersize);
   if( doWrite("mu_lepMVA_jetDR") ) tree->Branch("mu_lepMVA_jetDR", "std::vector<float>", &mu_lepMVA_jetDR, buffersize);
   if( doWrite("mu_lepMVA_jetPtRatio") ) tree->Branch("mu_lepMVA_jetPtRatio", "std::vector<float>", &mu_lepMVA_jetPtRatio, buffersize);
   if( doWrite("mu_lepMVA_jetBTagCSV") ) tree->Branch("mu_lepMVA_jetBTagCSV", "std::vector<float>", &mu_lepMVA_jetBTagCSV, buffersize);
   if( doWrite("mu_lepMVA_sip3d") ) tree->Branch("mu_lepMVA_sip3d", "std::vector<float>", &mu_lepMVA_sip3d, buffersize);
   if( doWrite("mu_lepMVA_dxy") ) tree->Branch("mu_lepMVA_dxy", "std::vector<float>", &mu_lepMVA_dxy, buffersize);
   if( doWrite("mu_lepMVA_dz") ) tree->Branch("mu_lepMVA_dz", "std::vector<float>", &mu_lepMVA_dz, buffersize);
   if( doWrite("mu_lepMVA_mvaId") ) tree->Branch("mu_lepMVA_mvaId", "std::vector<float>", &mu_lepMVA_mvaId, buffersize);

   if( doWrite("mu_hasMCMatch") ) tree->Branch("mu_hasMCMatch", "std::vector<int>", &mu_hasMCMatch, buffersize);
   if( doWrite("mu_gen_pt") ) tree->Branch("mu_gen_pt", "std::vector<float>", &mu_gen_pt, buffersize);
   if( doWrite("mu_gen_eta") ) tree->Branch("mu_gen_eta", "std::vector<float>", &mu_gen_eta, buffersize);
   if( doWrite("mu_gen_phi") ) tree->Branch("mu_gen_phi", "std::vector<float>", &mu_gen_phi, buffersize);
   if( doWrite("mu_gen_m") ) tree->Branch("mu_gen_m", "std::vector<float>", &mu_gen_m, buffersize);
   if( doWrite("mu_gen_status") ) tree->Branch("mu_gen_status", "std::vector<int>", &mu_gen_status, buffersize);
   if( doWrite("mu_gen_id") ) tree->Branch("mu_gen_id", "std::vector<int>", &mu_gen_id, buffersize);
   if( doWrite("mu_gen_charge") ) tree->Branch("mu_gen_charge", "std::vector<int>", &mu_gen_charge, buffersize);
   if( doWrite("mu_gen_dr") ) tree->Branch("mu_gen_dr", "std::vector<float>", &mu_gen_dr, buffersize);

   if( doWrite("tau_n") ) tree->Branch("tau_n", &tau_n, "tau_n/I", buffersize);
   if( doWrite("tau_pt") ) tree->Branch("tau_pt", "std::vector<float>", &tau_pt, buffersize);
   if( doWrite("tau_eta") ) tree->Branch("tau_eta", "std::vector<float>", &tau_eta, buffersize);
   if( doWrite("tau_phi") ) tree->Branch("tau_phi", "std::vector<float>", &tau_phi, buffersize);
   if( doWrite("tau_m") ) tree->Branch("tau_m", "std::vector<float>", &tau_m, buffersize);
   if( doWrite("tau_E") ) tree->Branch("tau_E", "std::vector<float>", &tau_E, buffersize);
   if( doWrite("tau_id") ) tree->Branch("tau_id", "std::vector<int>", &tau_id, buffersize);
   if( doWrite("tau_charge") ) tree->Branch("tau_charge", "std::vector<int>", &tau_charge, buffersize);
   
   if( doWrite("tau_hasLeadChargedHadrCand") ) tree->Branch("tau_hasLeadChargedHadrCand", "std::vector<bool>", &tau_hasLeadChargedHadrCand, buffersize);
   if( doWrite("tau_leadingTrackPt") ) tree->Branch("tau_leadingTrackPt", "std::vector<float>", &tau_leadingTrackPt, buffersize);
   if( doWrite("tau_leadingTrackCharge") ) tree->Branch("tau_leadingTrackCharge", "std::vector<int>", &tau_leadingTrackCharge, buffersize);
   
   if( doWrite("tau_decayMode") ) tree->Branch("tau_decayMode", "std::vector<int>", &tau_decayMode, buffersize);
//   if( doWrite("tau_decayModeFindingOldDMs") ) tree->Branch("tau_decayModeFindingOldDMs", "std::vector<float>", &tau_decayModeFindingOldDMs, buffersize);
   if( doWrite("tau_decayModeFindingNewDMs") ) tree->Branch("tau_decayModeFindingNewDMs", "std::vector<float>", &tau_decayModeFindingNewDMs, buffersize);
   
   if( doWrite("tau_puCorrPtSum") ) tree->Branch("tau_puCorrPtSum", "std::vector<float>", &tau_puCorrPtSum, buffersize);
   if( doWrite("tau_neutralIsoPtSum") ) tree->Branch("tau_neutralIsoPtSum", "std::vector<float>", &tau_neutralIsoPtSum, buffersize);
   if( doWrite("tau_chargedIsoPtSum") ) tree->Branch("tau_chargedIsoPtSum", "std::vector<float>", &tau_chargedIsoPtSum, buffersize);
   if( doWrite("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits") ) tree->Branch("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits", "std::vector<float>", &tau_byCombinedIsolationDeltaBetaCorrRaw3Hits, buffersize);
   
   if( doWrite("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits") ) tree->Branch("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", "std::vector<float>", &tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, buffersize);
   if( doWrite("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits") ) tree->Branch("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", "std::vector<float>", &tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, buffersize);
   if( doWrite("tau_byTightCombinedIsolationDeltaBetaCorr3Hits") ) tree->Branch("tau_byTightCombinedIsolationDeltaBetaCorr3Hits", "std::vector<float>", &tau_byTightCombinedIsolationDeltaBetaCorr3Hits, buffersize);
   
   if( doWrite("tau_againstMuonLoose3") ) tree->Branch("tau_againstMuonLoose3", "std::vector<float>", &tau_againstMuonLoose3, buffersize);
   if( doWrite("tau_againstMuonTight3") ) tree->Branch("tau_againstMuonTight3", "std::vector<float>", &tau_againstMuonTight3, buffersize);
   
   if( doWrite("tau_againstElectronVLooseMVA5") ) tree->Branch("tau_againstElectronVLooseMVA5", "std::vector<float>", &tau_againstElectronVLooseMVA5, buffersize);
   if( doWrite("tau_againstElectronLooseMVA5") ) tree->Branch("tau_againstElectronLooseMVA5", "std::vector<float>", &tau_againstElectronLooseMVA5, buffersize);
   if( doWrite("tau_againstElectronMediumMVA5") ) tree->Branch("tau_againstElectronMediumMVA5", "std::vector<float>", &tau_againstElectronMediumMVA5, buffersize);
   
   if( doWrite("tau_pfEssential_jet_pt") ) tree->Branch("tau_pfEssential_jet_pt", "std::vector<float>", &tau_pfEssential_jet_pt, buffersize);
   if( doWrite("tau_pfEssential_jet_eta") ) tree->Branch("tau_pfEssential_jet_eta", "std::vector<float>", &tau_pfEssential_jet_eta, buffersize);
   if( doWrite("tau_pfEssential_jet_phi") ) tree->Branch("tau_pfEssential_jet_phi", "std::vector<float>", &tau_pfEssential_jet_phi, buffersize);
   if( doWrite("tau_pfEssential_jet_m") ) tree->Branch("tau_pfEssential_jet_m", "std::vector<float>", &tau_pfEssential_jet_m, buffersize);

   if( doWrite("tau_pfEssential_jetCorr_pt") ) tree->Branch("tau_pfEssential_jetCorr_pt", "std::vector<float>", &tau_pfEssential_jetCorr_pt, buffersize);
   if( doWrite("tau_pfEssential_jetCorr_eta") ) tree->Branch("tau_pfEssential_jetCorr_eta", "std::vector<float>", &tau_pfEssential_jetCorr_eta, buffersize);
   if( doWrite("tau_pfEssential_jetCorr_phi") ) tree->Branch("tau_pfEssential_jetCorr_phi", "std::vector<float>", &tau_pfEssential_jetCorr_phi, buffersize);
   if( doWrite("tau_pfEssential_jetCorr_m") ) tree->Branch("tau_pfEssential_jetCorr_m", "std::vector<float>", &tau_pfEssential_jetCorr_m, buffersize);
   
   if( doWrite("tau_pfEssential_hasSV") ) tree->Branch("tau_pfEssential_hasSV", "std::vector<bool>", &tau_pfEssential_hasSV, buffersize);
   if( doWrite("tau_pfEssential_sv_x") ) tree->Branch("tau_pfEssential_sv_x", "std::vector<float>", &tau_pfEssential_sv_x, buffersize);
   if( doWrite("tau_pfEssential_sv_y") ) tree->Branch("tau_pfEssential_sv_y", "std::vector<float>", &tau_pfEssential_sv_y, buffersize);
   if( doWrite("tau_pfEssential_sv_z") ) tree->Branch("tau_pfEssential_sv_z", "std::vector<float>", &tau_pfEssential_sv_z, buffersize);
   
   if( doWrite("tau_pfEssential_flightLengthSig") ) tree->Branch("tau_pfEssential_flightLengthSig", "std::vector<float>", &tau_pfEssential_flightLengthSig, buffersize);
   if( doWrite("tau_pfEssential_dxy") ) tree->Branch("tau_pfEssential_dxy", "std::vector<float>", &tau_pfEssential_dxy, buffersize);
   if( doWrite("tau_pfEssential_dxy_error") ) tree->Branch("tau_pfEssential_dxy_error", "std::vector<float>", &tau_pfEssential_dxy_error, buffersize);
   if( doWrite("tau_pfEssential_dxy_Sig") ) tree->Branch("tau_pfEssential_dxy_Sig", "std::vector<float>", &tau_pfEssential_dxy_Sig, buffersize);
   
   if( doWrite("jet_n") ) tree->Branch("jet_n", &jet_n, "jet_n/I", buffersize);
   if( doWrite("jet_pt") ) tree->Branch("jet_pt", "std::vector<float>", &jet_pt, buffersize);
   if( doWrite("jet_eta") ) tree->Branch("jet_eta", "std::vector<float>", &jet_eta, buffersize);
   if( doWrite("jet_phi") ) tree->Branch("jet_phi", "std::vector<float>", &jet_phi, buffersize);
   if( doWrite("jet_m") ) tree->Branch("jet_m", "std::vector<float>", &jet_m, buffersize);
   if( doWrite("jet_E") ) tree->Branch("jet_E", "std::vector<float>", &jet_E, buffersize);

   if( doWrite("jet_ntrk") ) tree->Branch("jet_ntrk", "std::vector<int>", &jet_ntrk, buffersize);

   if( doWrite("jet_JBP") ) tree->Branch("jet_JBP", "std::vector<float>", &jet_JBP, buffersize);
   if( doWrite("jet_JP") ) tree->Branch("jet_JP", "std::vector<float>", &jet_JP, buffersize);
   if( doWrite("jet_TCHP") ) tree->Branch("jet_TCHP", "std::vector<float>", &jet_TCHP, buffersize);
   if( doWrite("jet_TCHE") ) tree->Branch("jet_TCHE", "std::vector<float>", &jet_TCHE, buffersize);
   if( doWrite("jet_SSVHE") ) tree->Branch("jet_SSVHE", "std::vector<float>", &jet_SSVHE, buffersize);
   if( doWrite("jet_SSVHP") ) tree->Branch("jet_SSVHP", "std::vector<float>", &jet_SSVHP, buffersize);
   if( doWrite("jet_CMVA") ) tree->Branch("jet_CMVA", "std::vector<float>", &jet_CMVA, buffersize);
   if( doWrite("jet_CSV") ) tree->Branch("jet_CSV", "std::vector<float>", &jet_CSV, buffersize);
   if( doWrite("jet_CSVv2") ) tree->Branch("jet_CSVv2", "std::vector<float>", &jet_CSVv2, buffersize);
   if( doWrite("jet_partonFlavour") ) tree->Branch("jet_partonFlavour", "std::vector<int>", &jet_partonFlavour, buffersize);
   if( doWrite("jet_hadronFlavour") ) tree->Branch("jet_hadronFlavour", "std::vector<int>", &jet_hadronFlavour, buffersize);

   if( doWrite("jet_neutralHadronEnergy") ) tree->Branch("jet_neutralHadronEnergy", "std::vector<float>", &jet_neutralHadronEnergy, buffersize);
   if( doWrite("jet_neutralEmEnergy") ) tree->Branch("jet_neutralEmEnergy", "std::vector<float>", &jet_neutralEmEnergy, buffersize);
   if( doWrite("jet_chargedHadronEnergy") ) tree->Branch("jet_chargedHadronEnergy", "std::vector<float>", &jet_chargedHadronEnergy, buffersize);
   if( doWrite("jet_chargedEmEnergy") ) tree->Branch("jet_chargedEmEnergy", "std::vector<float>", &jet_chargedEmEnergy, buffersize);
   if( doWrite("jet_electronEnergy") ) tree->Branch("jet_electronEnergy", "std::vector<float>", &jet_electronEnergy, buffersize);
   if( doWrite("jet_muonEnergy") ) tree->Branch("jet_muonEnergy", "std::vector<float>", &jet_muonEnergy, buffersize);
   if( doWrite("jet_photonEnergy") ) tree->Branch("jet_photonEnergy", "std::vector<float>", &jet_photonEnergy, buffersize);

   if( doWrite("jet_chargedMultiplicity") ) tree->Branch("jet_chargedMultiplicity", "std::vector<int>", &jet_chargedMultiplicity, buffersize);
   if( doWrite("jet_neutralMultiplicity") ) tree->Branch("jet_neutralMultiplicity", "std::vector<int>", &jet_neutralMultiplicity, buffersize);
   if( doWrite("jet_chargedHadronMultiplicity") ) tree->Branch("jet_chargedHadronMultiplicity", "std::vector<int>", &jet_chargedHadronMultiplicity, buffersize);
   
   if( doWrite("jet_jecFactorUncorrected") ) tree->Branch("jet_jecFactorUncorrected", "std::vector<float>", &jet_jecFactorUncorrected, buffersize);
   if( doWrite("jet_jecFactorL1FastJet") ) tree->Branch("jet_jecFactorL1FastJet", "std::vector<float>", &jet_jecFactorL1FastJet, buffersize);
   if( doWrite("jet_jecFactorL2Relative") ) tree->Branch("jet_jecFactorL2Relative", "std::vector<float>", &jet_jecFactorL2Relative, buffersize);
   if( doWrite("jet_jecFactorL3Absolute") ) tree->Branch("jet_jecFactorL3Absolute", "std::vector<float>", &jet_jecFactorL3Absolute, buffersize);
   
   if( doWrite("jet_Unc") ) tree->Branch("jet_Unc", "std::vector<float>", &jet_Unc, buffersize);
   
   if( doWrite("jet_pileupJetId") ) tree->Branch("jet_pileupJetId", "std::vector<float>", &jet_pileupJetId, buffersize);

   if( doWrite("jet_hasGenJet") ) tree->Branch("jet_hasGenJet", "std::vector<bool>", &jet_hasGenJet, buffersize);   
   if( doWrite("jet_genJet_pt") ) tree->Branch("jet_genJet_pt", "std::vector<float>", &jet_genJet_pt, buffersize);
   if( doWrite("jet_genJet_eta") ) tree->Branch("jet_genJet_eta", "std::vector<float>", &jet_genJet_eta, buffersize);
   if( doWrite("jet_genJet_phi") ) tree->Branch("jet_genJet_phi", "std::vector<float>", &jet_genJet_phi, buffersize);
   if( doWrite("jet_genJet_m") ) tree->Branch("jet_genJet_m", "std::vector<float>", &jet_genJet_m, buffersize);
   if( doWrite("jet_genJet_E") ) tree->Branch("jet_genJet_E", "std::vector<float>", &jet_genJet_E, buffersize);
   if( doWrite("jet_genJet_status") ) tree->Branch("jet_genJet_status", "std::vector<int>", &jet_genJet_status, buffersize);
   if( doWrite("jet_genJet_id") ) tree->Branch("jet_genJet_id", "std::vector<int>", &jet_genJet_id, buffersize);

   if( doWrite("jet_hasGenParton") ) tree->Branch("jet_hasGenParton", "std::vector<bool>", &jet_hasGenParton, buffersize);   
   if( doWrite("jet_genParton_pt") ) tree->Branch("jet_genParton_pt", "std::vector<float>", &jet_genParton_pt, buffersize);
   if( doWrite("jet_genParton_eta") ) tree->Branch("jet_genParton_eta", "std::vector<float>", &jet_genParton_eta, buffersize);
   if( doWrite("jet_genParton_phi") ) tree->Branch("jet_genParton_phi", "std::vector<float>", &jet_genParton_phi, buffersize);
   if( doWrite("jet_genParton_m") ) tree->Branch("jet_genParton_m", "std::vector<float>", &jet_genParton_m, buffersize);
   if( doWrite("jet_genParton_E") ) tree->Branch("jet_genParton_E", "std::vector<float>", &jet_genParton_E, buffersize);
   if( doWrite("jet_genParton_status") ) tree->Branch("jet_genParton_status", "std::vector<int>", &jet_genParton_status, buffersize);
   if( doWrite("jet_genParton_id") ) tree->Branch("jet_genParton_id", "std::vector<int>", &jet_genParton_id, buffersize);
   
   if( doWrite("jetPuppi_n") ) tree->Branch("jetPuppi_n", &jetPuppi_n, "jetPuppi_n/I", buffersize);
   if( doWrite("jetPuppi_pt") ) tree->Branch("jetPuppi_pt", "std::vector<float>", &jetPuppi_pt, buffersize);
   if( doWrite("jetPuppi_eta") ) tree->Branch("jetPuppi_eta", "std::vector<float>", &jetPuppi_eta, buffersize);
   if( doWrite("jetPuppi_phi") ) tree->Branch("jetPuppi_phi", "std::vector<float>", &jetPuppi_phi, buffersize);
   if( doWrite("jetPuppi_m") ) tree->Branch("jetPuppi_m", "std::vector<float>", &jetPuppi_m, buffersize);
   if( doWrite("jetPuppi_E") ) tree->Branch("jetPuppi_E", "std::vector<float>", &jetPuppi_E, buffersize);

   if( doWrite("jetPuppi_ntrk") ) tree->Branch("jetPuppi_ntrk", "std::vector<int>", &jetPuppi_ntrk, buffersize);

   if( doWrite("jetPuppi_JBP") ) tree->Branch("jetPuppi_JBP", "std::vector<float>", &jetPuppi_JBP, buffersize);
   if( doWrite("jetPuppi_JP") ) tree->Branch("jetPuppi_JP", "std::vector<float>", &jetPuppi_JP, buffersize);
   if( doWrite("jetPuppi_TCHP") ) tree->Branch("jetPuppi_TCHP", "std::vector<float>", &jetPuppi_TCHP, buffersize);
   if( doWrite("jetPuppi_TCHE") ) tree->Branch("jetPuppi_TCHE", "std::vector<float>", &jetPuppi_TCHE, buffersize);
   if( doWrite("jetPuppi_SSVHE") ) tree->Branch("jetPuppi_SSVHE", "std::vector<float>", &jetPuppi_SSVHE, buffersize);
   if( doWrite("jetPuppi_SSVHP") ) tree->Branch("jetPuppi_SSVHP", "std::vector<float>", &jetPuppi_SSVHP, buffersize);
   if( doWrite("jetPuppi_CMVA") ) tree->Branch("jetPuppi_CMVA", "std::vector<float>", &jetPuppi_CMVA, buffersize);
   if( doWrite("jetPuppi_CSV") ) tree->Branch("jetPuppi_CSV", "std::vector<float>", &jetPuppi_CSV, buffersize);
   if( doWrite("jetPuppi_CSVv2") ) tree->Branch("jetPuppi_CSVv2", "std::vector<float>", &jetPuppi_CSVv2, buffersize);
   if( doWrite("jetPuppi_partonFlavour") ) tree->Branch("jetPuppi_partonFlavour", "std::vector<int>", &jetPuppi_partonFlavour, buffersize);
   if( doWrite("jetPuppi_hadronFlavour") ) tree->Branch("jetPuppi_hadronFlavour", "std::vector<int>", &jetPuppi_hadronFlavour, buffersize);

   if( doWrite("jetPuppi_neutralHadronEnergy") ) tree->Branch("jetPuppi_neutralHadronEnergy", "std::vector<float>", &jetPuppi_neutralHadronEnergy, buffersize);
   if( doWrite("jetPuppi_neutralEmEnergy") ) tree->Branch("jetPuppi_neutralEmEnergy", "std::vector<float>", &jetPuppi_neutralEmEnergy, buffersize);
   if( doWrite("jetPuppi_chargedHadronEnergy") ) tree->Branch("jetPuppi_chargedHadronEnergy", "std::vector<float>", &jetPuppi_chargedHadronEnergy, buffersize);
   if( doWrite("jetPuppi_chargedEmEnergy") ) tree->Branch("jetPuppi_chargedEmEnergy", "std::vector<float>", &jetPuppi_chargedEmEnergy, buffersize);
   if( doWrite("jetPuppi_electronEnergy") ) tree->Branch("jetPuppi_electronEnergy", "std::vector<float>", &jetPuppi_electronEnergy, buffersize);
   if( doWrite("jetPuppi_muonEnergy") ) tree->Branch("jetPuppi_muonEnergy", "std::vector<float>", &jetPuppi_muonEnergy, buffersize);
   if( doWrite("jetPuppi_photonEnergy") ) tree->Branch("jetPuppi_photonEnergy", "std::vector<float>", &jetPuppi_photonEnergy, buffersize);

   if( doWrite("jetPuppi_chargedMultiplicity") ) tree->Branch("jetPuppi_chargedMultiplicity", "std::vector<int>", &jetPuppi_chargedMultiplicity, buffersize);
   if( doWrite("jetPuppi_neutralMultiplicity") ) tree->Branch("jetPuppi_neutralMultiplicity", "std::vector<int>", &jetPuppi_neutralMultiplicity, buffersize);
   if( doWrite("jetPuppi_chargedHadronMultiplicity") ) tree->Branch("jetPuppi_chargedHadronMultiplicity", "std::vector<int>", &jetPuppi_chargedHadronMultiplicity, buffersize);
   
   if( doWrite("jetPuppi_jecFactorUncorrected") ) tree->Branch("jetPuppi_jecFactorUncorrected", "std::vector<float>", &jetPuppi_jecFactorUncorrected, buffersize);
   if( doWrite("jetPuppi_jecFactorL1FastJet") ) tree->Branch("jetPuppi_jecFactorL1FastJet", "std::vector<float>", &jetPuppi_jecFactorL1FastJet, buffersize);
   if( doWrite("jetPuppi_jecFactorL2Relative") ) tree->Branch("jetPuppi_jecFactorL2Relative", "std::vector<float>", &jetPuppi_jecFactorL2Relative, buffersize);
   if( doWrite("jetPuppi_jecFactorL3Absolute") ) tree->Branch("jetPuppi_jecFactorL3Absolute", "std::vector<float>", &jetPuppi_jecFactorL3Absolute, buffersize);
   
   if( doWrite("jetPuppi_pileupJetId") ) tree->Branch("jetPuppi_pileupJetId", "std::vector<float>", &jetPuppi_pileupJetId, buffersize);

   if( doWrite("jetPuppi_hasGenJet") ) tree->Branch("jetPuppi_hasGenJet", "std::vector<bool>", &jetPuppi_hasGenJet, buffersize);   
   if( doWrite("jetPuppi_genJet_pt") ) tree->Branch("jetPuppi_genJet_pt", "std::vector<float>", &jetPuppi_genJet_pt, buffersize);
   if( doWrite("jetPuppi_genJet_eta") ) tree->Branch("jetPuppi_genJet_eta", "std::vector<float>", &jetPuppi_genJet_eta, buffersize);
   if( doWrite("jetPuppi_genJet_phi") ) tree->Branch("jetPuppi_genJet_phi", "std::vector<float>", &jetPuppi_genJet_phi, buffersize);
   if( doWrite("jetPuppi_genJet_m") ) tree->Branch("jetPuppi_genJet_m", "std::vector<float>", &jetPuppi_genJet_m, buffersize);
   if( doWrite("jetPuppi_genJet_E") ) tree->Branch("jetPuppi_genJet_E", "std::vector<float>", &jetPuppi_genJet_E, buffersize);
   if( doWrite("jetPuppi_genJet_status") ) tree->Branch("jetPuppi_genJet_status", "std::vector<int>", &jetPuppi_genJet_status, buffersize);
   if( doWrite("jetPuppi_genJet_id") ) tree->Branch("jetPuppi_genJet_id", "std::vector<int>", &jetPuppi_genJet_id, buffersize);       

   if( doWrite("jetPuppi_hasGenParton") ) tree->Branch("jetPuppi_hasGenJet", "std::vector<bool>", &jetPuppi_hasGenParton, buffersize);   
   if( doWrite("jetPuppi_genParton_pt") ) tree->Branch("jetPuppi_genParton_pt", "std::vector<float>", &jetPuppi_genParton_pt, buffersize);
   if( doWrite("jetPuppi_genParton_eta") ) tree->Branch("jetPuppi_genParton_eta", "std::vector<float>", &jetPuppi_genParton_eta, buffersize);
   if( doWrite("jetPuppi_genParton_phi") ) tree->Branch("jetPuppi_genParton_phi", "std::vector<float>", &jetPuppi_genParton_phi, buffersize);
   if( doWrite("jetPuppi_genParton_m") ) tree->Branch("jetPuppi_genParton_m", "std::vector<float>", &jetPuppi_genParton_m, buffersize);
   if( doWrite("jetPuppi_genParton_E") ) tree->Branch("jetPuppi_genParton_E", "std::vector<float>", &jetPuppi_genParton_E, buffersize);
   if( doWrite("jetPuppi_genParton_status") ) tree->Branch("jetPuppi_genParton_status", "std::vector<int>", &jetPuppi_genParton_status, buffersize);
   if( doWrite("jetPuppi_genParton_id") ) tree->Branch("jetPuppi_genParton_id", "std::vector<int>", &jetPuppi_genParton_id, buffersize);       
   
   if( doWrite("genJet_n") ) tree->Branch("genJet_n", &genJet_n, "genJet_n/I", buffersize);
   if( doWrite("genJet_pt") ) tree->Branch("genJet_pt", "std::vector<float>", &genJet_pt, buffersize);
   if( doWrite("genJet_eta") ) tree->Branch("genJet_eta", "std::vector<float>", &genJet_eta, buffersize);
   if( doWrite("genJet_phi") ) tree->Branch("genJet_phi", "std::vector<float>", &genJet_phi, buffersize);
   if( doWrite("genJet_m") ) tree->Branch("genJet_m", "std::vector<float>", &genJet_m, buffersize);
   if( doWrite("genJet_E") ) tree->Branch("genJet_E", "std::vector<float>", &genJet_E, buffersize);
   if( doWrite("genJet_emEnergy") ) tree->Branch("genJet_emEnergy", "std::vector<float>", &genJet_emEnergy, buffersize);
   if( doWrite("genJet_hadEnergy") ) tree->Branch("genJet_hadEnergy", "std::vector<float>", &genJet_hadEnergy, buffersize);
   if( doWrite("genJet_invisibleEnergy") ) tree->Branch("genJet_invisibleEnergy", "std::vector<float>", &genJet_invisibleEnergy, buffersize);
   if( doWrite("genJet_auxiliaryEnergy") ) tree->Branch("genJet_auxiliaryEnergy", "std::vector<float>", &genJet_auxiliaryEnergy, buffersize);
   if( doWrite("genJet_flavour") ) tree->Branch("genJet_flavour", "std::vector<int>", &genJet_flavour, buffersize);
   
   if( doWrite("mc_truth_tth") )
     {
	tree->Branch("mc_truth_tth_channel", &mc_truth_tth_channel, "mc_truth_tth_channel/I", buffersize);

	if( doWrite("mc_truth_p4") )
	  {	
	     tree->Branch("mc_truth_h0_p4", "TLorentzVector", &mc_truth_h0_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0W1_p4", "TLorentzVector", &mc_truth_h0W1_p4, buffersize);
	     tree->Branch("mc_truth_h0W2_p4", "TLorentzVector", &mc_truth_h0W2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wl1_p4", "TLorentzVector", &mc_truth_h0Wl1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnu1_p4", "TLorentzVector", &mc_truth_h0Wnu1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtau1_p4", "TLorentzVector", &mc_truth_h0Wtau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnutau1_p4", "TLorentzVector", &mc_truth_h0Wnutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaul1_p4", "TLorentzVector", &mc_truth_h0Wtaul1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunu1_p4", "TLorentzVector", &mc_truth_h0Wtaunu1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunutau1_p4", "TLorentzVector", &mc_truth_h0Wtaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wl2_p4", "TLorentzVector", &mc_truth_h0Wl2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnu2_p4", "TLorentzVector", &mc_truth_h0Wnu2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtau2_p4", "TLorentzVector", &mc_truth_h0Wtau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnutau2_p4", "TLorentzVector", &mc_truth_h0Wnutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaul2_p4", "TLorentzVector", &mc_truth_h0Wtaul2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunu2_p4", "TLorentzVector", &mc_truth_h0Wtaunu2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunutau2_p4", "TLorentzVector", &mc_truth_h0Wtaunutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq11_p4", "TLorentzVector", &mc_truth_h0Wq11_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq21_p4", "TLorentzVector", &mc_truth_h0Wq21_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq12_p4", "TLorentzVector", &mc_truth_h0Wq12_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq22_p4", "TLorentzVector", &mc_truth_h0Wq22_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0Z1_p4", "TLorentzVector", &mc_truth_h0Z1_p4, buffersize);
	     tree->Branch("mc_truth_h0Z2_p4", "TLorentzVector", &mc_truth_h0Z2_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl11_p4", "TLorentzVector", &mc_truth_h0Zl11_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl21_p4", "TLorentzVector", &mc_truth_h0Zl21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl12_p4", "TLorentzVector", &mc_truth_h0Zl12_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl22_p4", "TLorentzVector", &mc_truth_h0Zl22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau11_p4", "TLorentzVector", &mc_truth_h0Ztau11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau21_p4", "TLorentzVector", &mc_truth_h0Ztau21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul11_p4", "TLorentzVector", &mc_truth_h0Ztaul11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul21_p4", "TLorentzVector", &mc_truth_h0Ztaul21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu11_p4", "TLorentzVector", &mc_truth_h0Ztaunu11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu21_p4", "TLorentzVector", &mc_truth_h0Ztaunu21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau11_p4", "TLorentzVector", &mc_truth_h0Ztaunutau11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau21_p4", "TLorentzVector", &mc_truth_h0Ztaunutau21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq11_p4", "TLorentzVector", &mc_truth_h0Zq11_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq21_p4", "TLorentzVector", &mc_truth_h0Zq21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq12_p4", "TLorentzVector", &mc_truth_h0Zq12_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq22_p4", "TLorentzVector", &mc_truth_h0Zq22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau12_p4", "TLorentzVector", &mc_truth_h0Ztau12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau22_p4", "TLorentzVector", &mc_truth_h0Ztau22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul12_p4", "TLorentzVector", &mc_truth_h0Ztaul12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul22_p4", "TLorentzVector", &mc_truth_h0Ztaul22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu12_p4", "TLorentzVector", &mc_truth_h0Ztaunu12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu22_p4", "TLorentzVector", &mc_truth_h0Ztaunu22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau12_p4", "TLorentzVector", &mc_truth_h0Ztaunutau12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau22_p4", "TLorentzVector", &mc_truth_h0Ztaunutau22_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq12_p4", "TLorentzVector", &mc_truth_h0Zq12_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq22_p4", "TLorentzVector", &mc_truth_h0Zq22_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu11_p4", "TLorentzVector", &mc_truth_h0Znu11_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu21_p4", "TLorentzVector", &mc_truth_h0Znu21_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu12_p4", "TLorentzVector", &mc_truth_h0Znu12_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu22_p4", "TLorentzVector", &mc_truth_h0Znu22_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0tau1_p4", "TLorentzVector", &mc_truth_h0tau1_p4, buffersize);
	     tree->Branch("mc_truth_h0tau2_p4", "TLorentzVector", &mc_truth_h0tau2_p4, buffersize);
	     tree->Branch("mc_truth_h0taul1_p4", "TLorentzVector", &mc_truth_h0taul1_p4, buffersize);
	     tree->Branch("mc_truth_h0taunutau1_p4", "TLorentzVector", &mc_truth_h0taunutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0taunu1_p4", "TLorentzVector", &mc_truth_h0taunu1_p4, buffersize);
	     tree->Branch("mc_truth_h0taul2_p4", "TLorentzVector", &mc_truth_h0taul2_p4, buffersize);
	     tree->Branch("mc_truth_h0taunutau2_p4", "TLorentzVector", &mc_truth_h0taunutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0taunu2_p4", "TLorentzVector", &mc_truth_h0taunu2_p4, buffersize);
	     
	     tree->Branch("mc_truth_t1_p4", "TLorentzVector", &mc_truth_t1_p4, buffersize);
	     tree->Branch("mc_truth_t2_p4", "TLorentzVector", &mc_truth_t2_p4, buffersize);
	     tree->Branch("mc_truth_tb1_p4", "TLorentzVector", &mc_truth_tb1_p4, buffersize);
	     tree->Branch("mc_truth_tb2_p4", "TLorentzVector", &mc_truth_tb2_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW1_p4", "TLorentzVector", &mc_truth_tW1_p4, buffersize);
	     tree->Branch("mc_truth_tWnu1_p4", "TLorentzVector", &mc_truth_tWnu1_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau1_p4", "TLorentzVector", &mc_truth_tWnutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWl1_p4", "TLorentzVector", &mc_truth_tWl1_p4, buffersize);
	     tree->Branch("mc_truth_tWtau1_p4", "TLorentzVector", &mc_truth_tWtau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu1_p4", "TLorentzVector", &mc_truth_tWtaunu1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau1_p4", "TLorentzVector", &mc_truth_tWtaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul1_p4", "TLorentzVector", &mc_truth_tWtaul1_p4, buffersize);
	     tree->Branch("mc_truth_tWq11_p4", "TLorentzVector", &mc_truth_tWq11_p4, buffersize);
	     tree->Branch("mc_truth_tWq21_p4", "TLorentzVector", &mc_truth_tWq21_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW2_p4", "TLorentzVector", &mc_truth_tW2_p4, buffersize);
	     tree->Branch("mc_truth_tWnu2_p4", "TLorentzVector", &mc_truth_tWnu2_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau2_p4", "TLorentzVector", &mc_truth_tWnutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWl2_p4", "TLorentzVector", &mc_truth_tWl2_p4, buffersize);
	     tree->Branch("mc_truth_tWtau2_p4", "TLorentzVector", &mc_truth_tWtau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu2_p4", "TLorentzVector", &mc_truth_tWtaunu2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau2_p4", "TLorentzVector", &mc_truth_tWtaunutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul2_p4", "TLorentzVector", &mc_truth_tWtaul2_p4, buffersize);
	     tree->Branch("mc_truth_tWq12_p4", "TLorentzVector", &mc_truth_tWq12_p4, buffersize);
	     tree->Branch("mc_truth_tWq22_p4", "TLorentzVector", &mc_truth_tWq22_p4, buffersize);
	     
	     tree->Branch("mc_truth_j1_p4", "TLorentzVector", &mc_truth_j1_p4, buffersize);
	     tree->Branch("mc_truth_j2_p4", "TLorentzVector", &mc_truth_j2_p4, buffersize);
	     tree->Branch("mc_truth_j3_p4", "TLorentzVector", &mc_truth_j3_p4, buffersize);
	  }

	tree->Branch("mc_truth_h0_pt", &mc_truth_h0_pt, "mc_truth_h0_pt/F", buffersize);

	tree->Branch("mc_truth_h0W1_pt", &mc_truth_h0W1_pt, "mc_truth_h0W1_pt/F", buffersize);
	tree->Branch("mc_truth_h0W2_pt", &mc_truth_h0W2_pt, "mc_truth_h0W2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_pt", &mc_truth_h0Wl1_pt, "mc_truth_h0Wl1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_pt", &mc_truth_h0Wnu1_pt, "mc_truth_h0Wnu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_pt", &mc_truth_h0Wtau1_pt, "mc_truth_h0Wtau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_pt", &mc_truth_h0Wnutau1_pt, "mc_truth_h0Wnutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_pt", &mc_truth_h0Wtaul1_pt, "mc_truth_h0Wtaul1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_pt", &mc_truth_h0Wtaunu1_pt, "mc_truth_h0Wtaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_pt", &mc_truth_h0Wtaunutau1_pt, "mc_truth_h0Wtaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_pt", &mc_truth_h0Wl2_pt, "mc_truth_h0Wl2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_pt", &mc_truth_h0Wnu2_pt, "mc_truth_h0Wnu2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_pt", &mc_truth_h0Wtau2_pt, "mc_truth_h0Wtau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_pt", &mc_truth_h0Wnutau2_pt, "mc_truth_h0Wnutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_pt", &mc_truth_h0Wtaul2_pt, "mc_truth_h0Wtaul2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_pt", &mc_truth_h0Wtaunu2_pt, "mc_truth_h0Wtaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_pt", &mc_truth_h0Wtaunutau2_pt, "mc_truth_h0Wtaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_pt", &mc_truth_h0Wq11_pt, "mc_truth_h0Wq11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_pt", &mc_truth_h0Wq21_pt, "mc_truth_h0Wq21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_pt", &mc_truth_h0Wq12_pt, "mc_truth_h0Wq12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_pt", &mc_truth_h0Wq22_pt, "mc_truth_h0Wq22_pt/F", buffersize);

	tree->Branch("mc_truth_h0Z1_pt", &mc_truth_h0Z1_pt, "mc_truth_h0Z1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Z2_pt", &mc_truth_h0Z2_pt, "mc_truth_h0Z2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_pt", &mc_truth_h0Zl11_pt, "mc_truth_h0Zl11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_pt", &mc_truth_h0Zl21_pt, "mc_truth_h0Zl21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_pt", &mc_truth_h0Ztau11_pt, "mc_truth_h0Ztau11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_pt", &mc_truth_h0Ztau21_pt, "mc_truth_h0Ztau21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_pt", &mc_truth_h0Ztaul11_pt, "mc_truth_h0Ztaul11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_pt", &mc_truth_h0Ztaul21_pt, "mc_truth_h0Ztaul21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_pt", &mc_truth_h0Ztaunu11_pt, "mc_truth_h0Ztaunu11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_pt", &mc_truth_h0Ztaunu21_pt, "mc_truth_h0Ztaunu21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_pt", &mc_truth_h0Ztaunutau11_pt, "mc_truth_h0Ztaunutau11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_pt", &mc_truth_h0Ztaunutau21_pt, "mc_truth_h0Ztaunutau21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_pt", &mc_truth_h0Zq11_pt, "mc_truth_h0Zq11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_pt", &mc_truth_h0Zq21_pt, "mc_truth_h0Zq21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_pt", &mc_truth_h0Zl12_pt, "mc_truth_h0Zl12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_pt", &mc_truth_h0Zl22_pt, "mc_truth_h0Zl22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_pt", &mc_truth_h0Ztau12_pt, "mc_truth_h0Ztau12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_pt", &mc_truth_h0Ztau22_pt, "mc_truth_h0Ztau22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_pt", &mc_truth_h0Ztaul12_pt, "mc_truth_h0Ztaul12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_pt", &mc_truth_h0Ztaul22_pt, "mc_truth_h0Ztaul22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_pt", &mc_truth_h0Ztaunu12_pt, "mc_truth_h0Ztaunu12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_pt", &mc_truth_h0Ztaunu22_pt, "mc_truth_h0Ztaunu22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_pt", &mc_truth_h0Ztaunutau12_pt, "mc_truth_h0Ztaunutau12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_pt", &mc_truth_h0Ztaunutau22_pt, "mc_truth_h0Ztaunutau22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_pt", &mc_truth_h0Zq12_pt, "mc_truth_h0Zq12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_pt", &mc_truth_h0Zq22_pt, "mc_truth_h0Zq22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_pt", &mc_truth_h0Znu11_pt, "mc_truth_h0Znu11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_pt", &mc_truth_h0Znu21_pt, "mc_truth_h0Znu21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_pt", &mc_truth_h0Znu12_pt, "mc_truth_h0Znu12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_pt", &mc_truth_h0Znu22_pt, "mc_truth_h0Znu22_pt/F", buffersize);

	tree->Branch("mc_truth_h0tau1_pt", &mc_truth_h0tau1_pt, "mc_truth_h0tau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0tau2_pt", &mc_truth_h0tau2_pt, "mc_truth_h0tau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taul1_pt", &mc_truth_h0taul1_pt, "mc_truth_h0taul1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_pt", &mc_truth_h0taunutau1_pt, "mc_truth_h0taunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_pt", &mc_truth_h0taunu1_pt, "mc_truth_h0taunu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taul2_pt", &mc_truth_h0taul2_pt, "mc_truth_h0taul2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_pt", &mc_truth_h0taunutau2_pt, "mc_truth_h0taunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_pt", &mc_truth_h0taunu2_pt, "mc_truth_h0taunu2_pt/F", buffersize);

	tree->Branch("mc_truth_t1_pt", &mc_truth_t1_pt, "mc_truth_t1_pt/F", buffersize);
	tree->Branch("mc_truth_t2_pt", &mc_truth_t2_pt, "mc_truth_t2_pt/F", buffersize);
	tree->Branch("mc_truth_tb1_pt", &mc_truth_tb1_pt, "mc_truth_tb1_pt/F", buffersize);
	tree->Branch("mc_truth_tb2_pt", &mc_truth_tb2_pt, "mc_truth_tb2_pt/F", buffersize);

	tree->Branch("mc_truth_tW1_pt", &mc_truth_tW1_pt, "mc_truth_tW1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu1_pt", &mc_truth_tWnu1_pt, "mc_truth_tWnu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_pt", &mc_truth_tWnutau1_pt, "mc_truth_tWnutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWl1_pt", &mc_truth_tWl1_pt, "mc_truth_tWl1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau1_pt", &mc_truth_tWtau1_pt, "mc_truth_tWtau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_pt", &mc_truth_tWtaunu1_pt, "mc_truth_tWtaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_pt", &mc_truth_tWtaunutau1_pt, "mc_truth_tWtaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_pt", &mc_truth_tWtaul1_pt, "mc_truth_tWtaul1_pt/F", buffersize);
	tree->Branch("mc_truth_tWq11_pt", &mc_truth_tWq11_pt, "mc_truth_tWq11_pt/F", buffersize);
	tree->Branch("mc_truth_tWq21_pt", &mc_truth_tWq21_pt, "mc_truth_tWq21_pt/F", buffersize);

	tree->Branch("mc_truth_tW2_pt", &mc_truth_tW2_pt, "mc_truth_tW2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu2_pt", &mc_truth_tWnu2_pt, "mc_truth_tWnu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_pt", &mc_truth_tWnutau2_pt, "mc_truth_tWnutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWl2_pt", &mc_truth_tWl2_pt, "mc_truth_tWl2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau2_pt", &mc_truth_tWtau2_pt, "mc_truth_tWtau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_pt", &mc_truth_tWtaunu2_pt, "mc_truth_tWtaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_pt", &mc_truth_tWtaunutau2_pt, "mc_truth_tWtaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_pt", &mc_truth_tWtaul2_pt, "mc_truth_tWtaul2_pt/F", buffersize);
	tree->Branch("mc_truth_tWq12_pt", &mc_truth_tWq12_pt, "mc_truth_tWq12_pt/F", buffersize);
	tree->Branch("mc_truth_tWq22_pt", &mc_truth_tWq22_pt, "mc_truth_tWq22_pt/F", buffersize);

	tree->Branch("mc_truth_j1_pt", &mc_truth_j1_pt, "mc_truth_j1_pt/F", buffersize);
	tree->Branch("mc_truth_j2_pt", &mc_truth_j2_pt, "mc_truth_j2_pt/F", buffersize);
	tree->Branch("mc_truth_j3_pt", &mc_truth_j3_pt, "mc_truth_j3_pt/F", buffersize);
		
	tree->Branch("mc_truth_h0_eta", &mc_truth_h0_eta, "mc_truth_h0_eta/F", buffersize);

	tree->Branch("mc_truth_h0W1_eta", &mc_truth_h0W1_eta, "mc_truth_h0W1_eta/F", buffersize);
	tree->Branch("mc_truth_h0W2_eta", &mc_truth_h0W2_eta, "mc_truth_h0W2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_eta", &mc_truth_h0Wl1_eta, "mc_truth_h0Wl1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_eta", &mc_truth_h0Wnu1_eta, "mc_truth_h0Wnu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_eta", &mc_truth_h0Wtau1_eta, "mc_truth_h0Wtau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_eta", &mc_truth_h0Wnutau1_eta, "mc_truth_h0Wnutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_eta", &mc_truth_h0Wtaul1_eta, "mc_truth_h0Wtaul1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_eta", &mc_truth_h0Wtaunu1_eta, "mc_truth_h0Wtaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_eta", &mc_truth_h0Wtaunutau1_eta, "mc_truth_h0Wtaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_eta", &mc_truth_h0Wl2_eta, "mc_truth_h0Wl2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_eta", &mc_truth_h0Wnu2_eta, "mc_truth_h0Wnu2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_eta", &mc_truth_h0Wtau2_eta, "mc_truth_h0Wtau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_eta", &mc_truth_h0Wnutau2_eta, "mc_truth_h0Wnutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_eta", &mc_truth_h0Wtaul2_eta, "mc_truth_h0Wtaul2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_eta", &mc_truth_h0Wtaunu2_eta, "mc_truth_h0Wtaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_eta", &mc_truth_h0Wtaunutau2_eta, "mc_truth_h0Wtaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_eta", &mc_truth_h0Wq11_eta, "mc_truth_h0Wq11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_eta", &mc_truth_h0Wq21_eta, "mc_truth_h0Wq21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_eta", &mc_truth_h0Wq12_eta, "mc_truth_h0Wq12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_eta", &mc_truth_h0Wq22_eta, "mc_truth_h0Wq22_eta/F", buffersize);

	tree->Branch("mc_truth_h0Z1_eta", &mc_truth_h0Z1_eta, "mc_truth_h0Z1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Z2_eta", &mc_truth_h0Z2_eta, "mc_truth_h0Z2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_eta", &mc_truth_h0Zl11_eta, "mc_truth_h0Zl11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_eta", &mc_truth_h0Zl21_eta, "mc_truth_h0Zl21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_eta", &mc_truth_h0Ztau11_eta, "mc_truth_h0Ztau11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_eta", &mc_truth_h0Ztau21_eta, "mc_truth_h0Ztau21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_eta", &mc_truth_h0Ztaul11_eta, "mc_truth_h0Ztaul11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_eta", &mc_truth_h0Ztaul21_eta, "mc_truth_h0Ztaul21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_eta", &mc_truth_h0Ztaunu11_eta, "mc_truth_h0Ztaunu11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_eta", &mc_truth_h0Ztaunu21_eta, "mc_truth_h0Ztaunu21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_eta", &mc_truth_h0Ztaunutau11_eta, "mc_truth_h0Ztaunutau11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_eta", &mc_truth_h0Ztaunutau21_eta, "mc_truth_h0Ztaunutau21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_eta", &mc_truth_h0Zq11_eta, "mc_truth_h0Zq11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_eta", &mc_truth_h0Zq21_eta, "mc_truth_h0Zq21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_eta", &mc_truth_h0Zl12_eta, "mc_truth_h0Zl12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_eta", &mc_truth_h0Zl22_eta, "mc_truth_h0Zl22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_eta", &mc_truth_h0Ztau12_eta, "mc_truth_h0Ztau12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_eta", &mc_truth_h0Ztau22_eta, "mc_truth_h0Ztau22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_eta", &mc_truth_h0Ztaul12_eta, "mc_truth_h0Ztaul12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_eta", &mc_truth_h0Ztaul22_eta, "mc_truth_h0Ztaul22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_eta", &mc_truth_h0Ztaunu12_eta, "mc_truth_h0Ztaunu12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_eta", &mc_truth_h0Ztaunu22_eta, "mc_truth_h0Ztaunu22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_eta", &mc_truth_h0Ztaunutau12_eta, "mc_truth_h0Ztaunutau12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_eta", &mc_truth_h0Ztaunutau22_eta, "mc_truth_h0Ztaunutau22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_eta", &mc_truth_h0Zq12_eta, "mc_truth_h0Zq12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_eta", &mc_truth_h0Zq22_eta, "mc_truth_h0Zq22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_eta", &mc_truth_h0Znu11_eta, "mc_truth_h0Znu11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_eta", &mc_truth_h0Znu21_eta, "mc_truth_h0Znu21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_eta", &mc_truth_h0Znu12_eta, "mc_truth_h0Znu12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_eta", &mc_truth_h0Znu22_eta, "mc_truth_h0Znu22_eta/F", buffersize);

	tree->Branch("mc_truth_h0tau1_eta", &mc_truth_h0tau1_eta, "mc_truth_h0tau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0tau2_eta", &mc_truth_h0tau2_eta, "mc_truth_h0tau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taul1_eta", &mc_truth_h0taul1_eta, "mc_truth_h0taul1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_eta", &mc_truth_h0taunutau1_eta, "mc_truth_h0taunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_eta", &mc_truth_h0taunu1_eta, "mc_truth_h0taunu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taul2_eta", &mc_truth_h0taul2_eta, "mc_truth_h0taul2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_eta", &mc_truth_h0taunutau2_eta, "mc_truth_h0taunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_eta", &mc_truth_h0taunu2_eta, "mc_truth_h0taunu2_eta/F", buffersize);

	tree->Branch("mc_truth_t1_eta", &mc_truth_t1_eta, "mc_truth_t1_eta/F", buffersize);
	tree->Branch("mc_truth_t2_eta", &mc_truth_t2_eta, "mc_truth_t2_eta/F", buffersize);
	tree->Branch("mc_truth_tb1_eta", &mc_truth_tb1_eta, "mc_truth_tb1_eta/F", buffersize);
	tree->Branch("mc_truth_tb2_eta", &mc_truth_tb2_eta, "mc_truth_tb2_eta/F", buffersize);

	tree->Branch("mc_truth_tW1_eta", &mc_truth_tW1_eta, "mc_truth_tW1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu1_eta", &mc_truth_tWnu1_eta, "mc_truth_tWnu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_eta", &mc_truth_tWnutau1_eta, "mc_truth_tWnutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWl1_eta", &mc_truth_tWl1_eta, "mc_truth_tWl1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau1_eta", &mc_truth_tWtau1_eta, "mc_truth_tWtau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_eta", &mc_truth_tWtaunu1_eta, "mc_truth_tWtaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_eta", &mc_truth_tWtaunutau1_eta, "mc_truth_tWtaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_eta", &mc_truth_tWtaul1_eta, "mc_truth_tWtaul1_eta/F", buffersize);
	tree->Branch("mc_truth_tWq11_eta", &mc_truth_tWq11_eta, "mc_truth_tWq11_eta/F", buffersize);
	tree->Branch("mc_truth_tWq21_eta", &mc_truth_tWq21_eta, "mc_truth_tWq21_eta/F", buffersize);

	tree->Branch("mc_truth_tW2_eta", &mc_truth_tW2_eta, "mc_truth_tW2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu2_eta", &mc_truth_tWnu2_eta, "mc_truth_tWnu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_eta", &mc_truth_tWnutau2_eta, "mc_truth_tWnutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWl2_eta", &mc_truth_tWl2_eta, "mc_truth_tWl2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau2_eta", &mc_truth_tWtau2_eta, "mc_truth_tWtau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_eta", &mc_truth_tWtaunu2_eta, "mc_truth_tWtaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_eta", &mc_truth_tWtaunutau2_eta, "mc_truth_tWtaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_eta", &mc_truth_tWtaul2_eta, "mc_truth_tWtaul2_eta/F", buffersize);
	tree->Branch("mc_truth_tWq12_eta", &mc_truth_tWq12_eta, "mc_truth_tWq12_eta/F", buffersize);
	tree->Branch("mc_truth_tWq22_eta", &mc_truth_tWq22_eta, "mc_truth_tWq22_eta/F", buffersize);

	tree->Branch("mc_truth_j1_eta", &mc_truth_j1_eta, "mc_truth_j1_eta/F", buffersize);
	tree->Branch("mc_truth_j2_eta", &mc_truth_j2_eta, "mc_truth_j2_eta/F", buffersize);
	tree->Branch("mc_truth_j3_eta", &mc_truth_j3_eta, "mc_truth_j3_eta/F", buffersize);
		
	tree->Branch("mc_truth_h0_phi", &mc_truth_h0_phi, "mc_truth_h0_phi/F", buffersize);

	tree->Branch("mc_truth_h0W1_phi", &mc_truth_h0W1_phi, "mc_truth_h0W1_phi/F", buffersize);
	tree->Branch("mc_truth_h0W2_phi", &mc_truth_h0W2_phi, "mc_truth_h0W2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_phi", &mc_truth_h0Wl1_phi, "mc_truth_h0Wl1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_phi", &mc_truth_h0Wnu1_phi, "mc_truth_h0Wnu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_phi", &mc_truth_h0Wtau1_phi, "mc_truth_h0Wtau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_phi", &mc_truth_h0Wnutau1_phi, "mc_truth_h0Wnutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_phi", &mc_truth_h0Wtaul1_phi, "mc_truth_h0Wtaul1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_phi", &mc_truth_h0Wtaunu1_phi, "mc_truth_h0Wtaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_phi", &mc_truth_h0Wtaunutau1_phi, "mc_truth_h0Wtaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_phi", &mc_truth_h0Wl2_phi, "mc_truth_h0Wl2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_phi", &mc_truth_h0Wnu2_phi, "mc_truth_h0Wnu2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_phi", &mc_truth_h0Wtau2_phi, "mc_truth_h0Wtau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_phi", &mc_truth_h0Wnutau2_phi, "mc_truth_h0Wnutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_phi", &mc_truth_h0Wtaul2_phi, "mc_truth_h0Wtaul2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_phi", &mc_truth_h0Wtaunu2_phi, "mc_truth_h0Wtaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_phi", &mc_truth_h0Wtaunutau2_phi, "mc_truth_h0Wtaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_phi", &mc_truth_h0Wq11_phi, "mc_truth_h0Wq11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_phi", &mc_truth_h0Wq21_phi, "mc_truth_h0Wq21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_phi", &mc_truth_h0Wq12_phi, "mc_truth_h0Wq12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_phi", &mc_truth_h0Wq22_phi, "mc_truth_h0Wq22_phi/F", buffersize);

	tree->Branch("mc_truth_h0Z1_phi", &mc_truth_h0Z1_phi, "mc_truth_h0Z1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Z2_phi", &mc_truth_h0Z2_phi, "mc_truth_h0Z2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_phi", &mc_truth_h0Zl11_phi, "mc_truth_h0Zl11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_phi", &mc_truth_h0Zl21_phi, "mc_truth_h0Zl21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_phi", &mc_truth_h0Ztau11_phi, "mc_truth_h0Ztau11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_phi", &mc_truth_h0Ztau21_phi, "mc_truth_h0Ztau21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_phi", &mc_truth_h0Ztaul11_phi, "mc_truth_h0Ztaul11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_phi", &mc_truth_h0Ztaul21_phi, "mc_truth_h0Ztaul21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_phi", &mc_truth_h0Ztaunu11_phi, "mc_truth_h0Ztaunu11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_phi", &mc_truth_h0Ztaunu21_phi, "mc_truth_h0Ztaunu21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_phi", &mc_truth_h0Ztaunutau11_phi, "mc_truth_h0Ztaunutau11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_phi", &mc_truth_h0Ztaunutau21_phi, "mc_truth_h0Ztaunutau21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_phi", &mc_truth_h0Zq11_phi, "mc_truth_h0Zq11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_phi", &mc_truth_h0Zq21_phi, "mc_truth_h0Zq21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_phi", &mc_truth_h0Zl12_phi, "mc_truth_h0Zl12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_phi", &mc_truth_h0Zl22_phi, "mc_truth_h0Zl22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_phi", &mc_truth_h0Ztau12_phi, "mc_truth_h0Ztau12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_phi", &mc_truth_h0Ztau22_phi, "mc_truth_h0Ztau22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_phi", &mc_truth_h0Ztaul12_phi, "mc_truth_h0Ztaul12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_phi", &mc_truth_h0Ztaul22_phi, "mc_truth_h0Ztaul22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_phi", &mc_truth_h0Ztaunu12_phi, "mc_truth_h0Ztaunu12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_phi", &mc_truth_h0Ztaunu22_phi, "mc_truth_h0Ztaunu22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_phi", &mc_truth_h0Ztaunutau12_phi, "mc_truth_h0Ztaunutau12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_phi", &mc_truth_h0Ztaunutau22_phi, "mc_truth_h0Ztaunutau22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_phi", &mc_truth_h0Zq12_phi, "mc_truth_h0Zq12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_phi", &mc_truth_h0Zq22_phi, "mc_truth_h0Zq22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_phi", &mc_truth_h0Znu11_phi, "mc_truth_h0Znu11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_phi", &mc_truth_h0Znu21_phi, "mc_truth_h0Znu21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_phi", &mc_truth_h0Znu12_phi, "mc_truth_h0Znu12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_phi", &mc_truth_h0Znu22_phi, "mc_truth_h0Znu22_phi/F", buffersize);

	tree->Branch("mc_truth_h0tau1_phi", &mc_truth_h0tau1_phi, "mc_truth_h0tau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0tau2_phi", &mc_truth_h0tau2_phi, "mc_truth_h0tau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taul1_phi", &mc_truth_h0taul1_phi, "mc_truth_h0taul1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_phi", &mc_truth_h0taunutau1_phi, "mc_truth_h0taunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_phi", &mc_truth_h0taunu1_phi, "mc_truth_h0taunu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taul2_phi", &mc_truth_h0taul2_phi, "mc_truth_h0taul2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_phi", &mc_truth_h0taunutau2_phi, "mc_truth_h0taunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_phi", &mc_truth_h0taunu2_phi, "mc_truth_h0taunu2_phi/F", buffersize);

	tree->Branch("mc_truth_t1_phi", &mc_truth_t1_phi, "mc_truth_t1_phi/F", buffersize);
	tree->Branch("mc_truth_t2_phi", &mc_truth_t2_phi, "mc_truth_t2_phi/F", buffersize);
	tree->Branch("mc_truth_tb1_phi", &mc_truth_tb1_phi, "mc_truth_tb1_phi/F", buffersize);
	tree->Branch("mc_truth_tb2_phi", &mc_truth_tb2_phi, "mc_truth_tb2_phi/F", buffersize);

	tree->Branch("mc_truth_tW1_phi", &mc_truth_tW1_phi, "mc_truth_tW1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu1_phi", &mc_truth_tWnu1_phi, "mc_truth_tWnu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_phi", &mc_truth_tWnutau1_phi, "mc_truth_tWnutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWl1_phi", &mc_truth_tWl1_phi, "mc_truth_tWl1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau1_phi", &mc_truth_tWtau1_phi, "mc_truth_tWtau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_phi", &mc_truth_tWtaunu1_phi, "mc_truth_tWtaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_phi", &mc_truth_tWtaunutau1_phi, "mc_truth_tWtaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_phi", &mc_truth_tWtaul1_phi, "mc_truth_tWtaul1_phi/F", buffersize);
	tree->Branch("mc_truth_tWq11_phi", &mc_truth_tWq11_phi, "mc_truth_tWq11_phi/F", buffersize);
	tree->Branch("mc_truth_tWq21_phi", &mc_truth_tWq21_phi, "mc_truth_tWq21_phi/F", buffersize);

	tree->Branch("mc_truth_tW2_phi", &mc_truth_tW2_phi, "mc_truth_tW2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu2_phi", &mc_truth_tWnu2_phi, "mc_truth_tWnu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_phi", &mc_truth_tWnutau2_phi, "mc_truth_tWnutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWl2_phi", &mc_truth_tWl2_phi, "mc_truth_tWl2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau2_phi", &mc_truth_tWtau2_phi, "mc_truth_tWtau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_phi", &mc_truth_tWtaunu2_phi, "mc_truth_tWtaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_phi", &mc_truth_tWtaunutau2_phi, "mc_truth_tWtaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_phi", &mc_truth_tWtaul2_phi, "mc_truth_tWtaul2_phi/F", buffersize);
	tree->Branch("mc_truth_tWq12_phi", &mc_truth_tWq12_phi, "mc_truth_tWq12_phi/F", buffersize);
	tree->Branch("mc_truth_tWq22_phi", &mc_truth_tWq22_phi, "mc_truth_tWq22_phi/F", buffersize);

	tree->Branch("mc_truth_j1_phi", &mc_truth_j1_phi, "mc_truth_j1_phi/F", buffersize);
	tree->Branch("mc_truth_j2_phi", &mc_truth_j2_phi, "mc_truth_j2_phi/F", buffersize);
	tree->Branch("mc_truth_j3_phi", &mc_truth_j3_phi, "mc_truth_j3_phi/F", buffersize);
	
	tree->Branch("mc_truth_h0_E", &mc_truth_h0_E, "mc_truth_h0_E/F", buffersize);

	tree->Branch("mc_truth_h0W1_E", &mc_truth_h0W1_E, "mc_truth_h0W1_E/F", buffersize);
	tree->Branch("mc_truth_h0W2_E", &mc_truth_h0W2_E, "mc_truth_h0W2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_E", &mc_truth_h0Wl1_E, "mc_truth_h0Wl1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_E", &mc_truth_h0Wnu1_E, "mc_truth_h0Wnu1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_E", &mc_truth_h0Wtau1_E, "mc_truth_h0Wtau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_E", &mc_truth_h0Wnutau1_E, "mc_truth_h0Wnutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_E", &mc_truth_h0Wtaul1_E, "mc_truth_h0Wtaul1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_E", &mc_truth_h0Wtaunu1_E, "mc_truth_h0Wtaunu1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_E", &mc_truth_h0Wtaunutau1_E, "mc_truth_h0Wtaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_E", &mc_truth_h0Wl2_E, "mc_truth_h0Wl2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_E", &mc_truth_h0Wnu2_E, "mc_truth_h0Wnu2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_E", &mc_truth_h0Wtau2_E, "mc_truth_h0Wtau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_E", &mc_truth_h0Wnutau2_E, "mc_truth_h0Wnutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_E", &mc_truth_h0Wtaul2_E, "mc_truth_h0Wtaul2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_E", &mc_truth_h0Wtaunu2_E, "mc_truth_h0Wtaunu2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_E", &mc_truth_h0Wtaunutau2_E, "mc_truth_h0Wtaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_E", &mc_truth_h0Wq11_E, "mc_truth_h0Wq11_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_E", &mc_truth_h0Wq21_E, "mc_truth_h0Wq21_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_E", &mc_truth_h0Wq12_E, "mc_truth_h0Wq12_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_E", &mc_truth_h0Wq22_E, "mc_truth_h0Wq22_E/F", buffersize);

	tree->Branch("mc_truth_h0Z1_E", &mc_truth_h0Z1_E, "mc_truth_h0Z1_E/F", buffersize);
	tree->Branch("mc_truth_h0Z2_E", &mc_truth_h0Z2_E, "mc_truth_h0Z2_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_E", &mc_truth_h0Zl11_E, "mc_truth_h0Zl11_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_E", &mc_truth_h0Zl21_E, "mc_truth_h0Zl21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_E", &mc_truth_h0Ztau11_E, "mc_truth_h0Ztau11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_E", &mc_truth_h0Ztau21_E, "mc_truth_h0Ztau21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_E", &mc_truth_h0Ztaul11_E, "mc_truth_h0Ztaul11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_E", &mc_truth_h0Ztaul21_E, "mc_truth_h0Ztaul21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_E", &mc_truth_h0Ztaunu11_E, "mc_truth_h0Ztaunu11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_E", &mc_truth_h0Ztaunu21_E, "mc_truth_h0Ztaunu21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_E", &mc_truth_h0Ztaunutau11_E, "mc_truth_h0Ztaunutau11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_E", &mc_truth_h0Ztaunutau21_E, "mc_truth_h0Ztaunutau21_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_E", &mc_truth_h0Zq11_E, "mc_truth_h0Zq11_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_E", &mc_truth_h0Zq21_E, "mc_truth_h0Zq21_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_E", &mc_truth_h0Zl12_E, "mc_truth_h0Zl12_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_E", &mc_truth_h0Zl22_E, "mc_truth_h0Zl22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_E", &mc_truth_h0Ztau12_E, "mc_truth_h0Ztau12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_E", &mc_truth_h0Ztau22_E, "mc_truth_h0Ztau22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_E", &mc_truth_h0Ztaul12_E, "mc_truth_h0Ztaul12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_E", &mc_truth_h0Ztaul22_E, "mc_truth_h0Ztaul22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_E", &mc_truth_h0Ztaunu12_E, "mc_truth_h0Ztaunu12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_E", &mc_truth_h0Ztaunu22_E, "mc_truth_h0Ztaunu22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_E", &mc_truth_h0Ztaunutau12_E, "mc_truth_h0Ztaunutau12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_E", &mc_truth_h0Ztaunutau22_E, "mc_truth_h0Ztaunutau22_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_E", &mc_truth_h0Zq12_E, "mc_truth_h0Zq12_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_E", &mc_truth_h0Zq22_E, "mc_truth_h0Zq22_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_E", &mc_truth_h0Znu11_E, "mc_truth_h0Znu11_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_E", &mc_truth_h0Znu21_E, "mc_truth_h0Znu21_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_E", &mc_truth_h0Znu12_E, "mc_truth_h0Znu12_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_E", &mc_truth_h0Znu22_E, "mc_truth_h0Znu22_E/F", buffersize);

	tree->Branch("mc_truth_h0tau1_E", &mc_truth_h0tau1_E, "mc_truth_h0tau1_E/F", buffersize);
	tree->Branch("mc_truth_h0tau2_E", &mc_truth_h0tau2_E, "mc_truth_h0tau2_E/F", buffersize);
	tree->Branch("mc_truth_h0taul1_E", &mc_truth_h0taul1_E, "mc_truth_h0taul1_E/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_E", &mc_truth_h0taunutau1_E, "mc_truth_h0taunutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_E", &mc_truth_h0taunu1_E, "mc_truth_h0taunu1_E/F", buffersize);
	tree->Branch("mc_truth_h0taul2_E", &mc_truth_h0taul2_E, "mc_truth_h0taul2_E/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_E", &mc_truth_h0taunutau2_E, "mc_truth_h0taunutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_E", &mc_truth_h0taunu2_E, "mc_truth_h0taunu2_E/F", buffersize);

	tree->Branch("mc_truth_t1_E", &mc_truth_t1_E, "mc_truth_t1_E/F", buffersize);
	tree->Branch("mc_truth_t2_E", &mc_truth_t2_E, "mc_truth_t2_E/F", buffersize);
	tree->Branch("mc_truth_tb1_E", &mc_truth_tb1_E, "mc_truth_tb1_E/F", buffersize);
	tree->Branch("mc_truth_tb2_E", &mc_truth_tb2_E, "mc_truth_tb2_E/F", buffersize);

	tree->Branch("mc_truth_tW1_E", &mc_truth_tW1_E, "mc_truth_tW1_E/F", buffersize);
	tree->Branch("mc_truth_tWnu1_E", &mc_truth_tWnu1_E, "mc_truth_tWnu1_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_E", &mc_truth_tWnutau1_E, "mc_truth_tWnutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWl1_E", &mc_truth_tWl1_E, "mc_truth_tWl1_E/F", buffersize);
	tree->Branch("mc_truth_tWtau1_E", &mc_truth_tWtau1_E, "mc_truth_tWtau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_E", &mc_truth_tWtaunu1_E, "mc_truth_tWtaunu1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_E", &mc_truth_tWtaunutau1_E, "mc_truth_tWtaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_E", &mc_truth_tWtaul1_E, "mc_truth_tWtaul1_E/F", buffersize);
	tree->Branch("mc_truth_tWq11_E", &mc_truth_tWq11_E, "mc_truth_tWq11_E/F", buffersize);
	tree->Branch("mc_truth_tWq21_E", &mc_truth_tWq21_E, "mc_truth_tWq21_E/F", buffersize);

	tree->Branch("mc_truth_tW2_E", &mc_truth_tW2_E, "mc_truth_tW2_E/F", buffersize);
	tree->Branch("mc_truth_tWnu2_E", &mc_truth_tWnu2_E, "mc_truth_tWnu2_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_E", &mc_truth_tWnutau2_E, "mc_truth_tWnutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWl2_E", &mc_truth_tWl2_E, "mc_truth_tWl2_E/F", buffersize);
	tree->Branch("mc_truth_tWtau2_E", &mc_truth_tWtau2_E, "mc_truth_tWtau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_E", &mc_truth_tWtaunu2_E, "mc_truth_tWtaunu2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_E", &mc_truth_tWtaunutau2_E, "mc_truth_tWtaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_E", &mc_truth_tWtaul2_E, "mc_truth_tWtaul2_E/F", buffersize);
	tree->Branch("mc_truth_tWq12_E", &mc_truth_tWq12_E, "mc_truth_tWq12_E/F", buffersize);
	tree->Branch("mc_truth_tWq22_E", &mc_truth_tWq22_E, "mc_truth_tWq22_E/F", buffersize);

	tree->Branch("mc_truth_j1_E", &mc_truth_j1_E, "mc_truth_j1_E/F", buffersize);
	tree->Branch("mc_truth_j2_E", &mc_truth_j2_E, "mc_truth_j2_E/F", buffersize);
	tree->Branch("mc_truth_j3_E", &mc_truth_j3_E, "mc_truth_j3_E/F", buffersize);
		
	tree->Branch("mc_truth_h0_id", &mc_truth_h0_id, "mc_truth_h0_id/I", buffersize);

	tree->Branch("mc_truth_h0W1_id", &mc_truth_h0W1_id, "mc_truth_h0W1_id/I", buffersize);
	tree->Branch("mc_truth_h0W2_id", &mc_truth_h0W2_id, "mc_truth_h0W2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wl1_id", &mc_truth_h0Wl1_id, "mc_truth_h0Wl1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnu1_id", &mc_truth_h0Wnu1_id, "mc_truth_h0Wnu1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtau1_id", &mc_truth_h0Wtau1_id, "mc_truth_h0Wtau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_id", &mc_truth_h0Wnutau1_id, "mc_truth_h0Wnutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_id", &mc_truth_h0Wtaul1_id, "mc_truth_h0Wtaul1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_id", &mc_truth_h0Wtaunu1_id, "mc_truth_h0Wtaunu1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_id", &mc_truth_h0Wtaunutau1_id, "mc_truth_h0Wtaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wl2_id", &mc_truth_h0Wl2_id, "mc_truth_h0Wl2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnu2_id", &mc_truth_h0Wnu2_id, "mc_truth_h0Wnu2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtau2_id", &mc_truth_h0Wtau2_id, "mc_truth_h0Wtau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_id", &mc_truth_h0Wnutau2_id, "mc_truth_h0Wnutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_id", &mc_truth_h0Wtaul2_id, "mc_truth_h0Wtaul2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_id", &mc_truth_h0Wtaunu2_id, "mc_truth_h0Wtaunu2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_id", &mc_truth_h0Wtaunutau2_id, "mc_truth_h0Wtaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_id", &mc_truth_h0Wq11_id, "mc_truth_h0Wq11_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_id", &mc_truth_h0Wq21_id, "mc_truth_h0Wq21_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_id", &mc_truth_h0Wq12_id, "mc_truth_h0Wq12_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_id", &mc_truth_h0Wq22_id, "mc_truth_h0Wq22_id/I", buffersize);

	tree->Branch("mc_truth_h0Z1_id", &mc_truth_h0Z1_id, "mc_truth_h0Z1_id/I", buffersize);
	tree->Branch("mc_truth_h0Z2_id", &mc_truth_h0Z2_id, "mc_truth_h0Z2_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl11_id", &mc_truth_h0Zl11_id, "mc_truth_h0Zl11_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl21_id", &mc_truth_h0Zl21_id, "mc_truth_h0Zl21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau11_id", &mc_truth_h0Ztau11_id, "mc_truth_h0Ztau11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau21_id", &mc_truth_h0Ztau21_id, "mc_truth_h0Ztau21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_id", &mc_truth_h0Ztaul11_id, "mc_truth_h0Ztaul11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_id", &mc_truth_h0Ztaul21_id, "mc_truth_h0Ztaul21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_id", &mc_truth_h0Ztaunu11_id, "mc_truth_h0Ztaunu11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_id", &mc_truth_h0Ztaunu21_id, "mc_truth_h0Ztaunu21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_id", &mc_truth_h0Ztaunutau11_id, "mc_truth_h0Ztaunutau11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_id", &mc_truth_h0Ztaunutau21_id, "mc_truth_h0Ztaunutau21_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_id", &mc_truth_h0Zq11_id, "mc_truth_h0Zq11_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_id", &mc_truth_h0Zq21_id, "mc_truth_h0Zq21_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl12_id", &mc_truth_h0Zl12_id, "mc_truth_h0Zl12_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl22_id", &mc_truth_h0Zl22_id, "mc_truth_h0Zl22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau12_id", &mc_truth_h0Ztau12_id, "mc_truth_h0Ztau12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau22_id", &mc_truth_h0Ztau22_id, "mc_truth_h0Ztau22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_id", &mc_truth_h0Ztaul12_id, "mc_truth_h0Ztaul12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_id", &mc_truth_h0Ztaul22_id, "mc_truth_h0Ztaul22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_id", &mc_truth_h0Ztaunu12_id, "mc_truth_h0Ztaunu12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_id", &mc_truth_h0Ztaunu22_id, "mc_truth_h0Ztaunu22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_id", &mc_truth_h0Ztaunutau12_id, "mc_truth_h0Ztaunutau12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_id", &mc_truth_h0Ztaunutau22_id, "mc_truth_h0Ztaunutau22_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_id", &mc_truth_h0Zq12_id, "mc_truth_h0Zq12_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_id", &mc_truth_h0Zq22_id, "mc_truth_h0Zq22_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu11_id", &mc_truth_h0Znu11_id, "mc_truth_h0Znu11_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu21_id", &mc_truth_h0Znu21_id, "mc_truth_h0Znu21_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu12_id", &mc_truth_h0Znu12_id, "mc_truth_h0Znu12_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu22_id", &mc_truth_h0Znu22_id, "mc_truth_h0Znu22_id/I", buffersize);

	tree->Branch("mc_truth_h0tau1_id", &mc_truth_h0tau1_id, "mc_truth_h0tau1_id/I", buffersize);
	tree->Branch("mc_truth_h0tau2_id", &mc_truth_h0tau2_id, "mc_truth_h0tau2_id/I", buffersize);
	tree->Branch("mc_truth_h0taul1_id", &mc_truth_h0taul1_id, "mc_truth_h0taul1_id/I", buffersize);
	tree->Branch("mc_truth_h0taunutau1_id", &mc_truth_h0taunutau1_id, "mc_truth_h0taunutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0taunu1_id", &mc_truth_h0taunu1_id, "mc_truth_h0taunu1_id/I", buffersize);
	tree->Branch("mc_truth_h0taul2_id", &mc_truth_h0taul2_id, "mc_truth_h0taul2_id/I", buffersize);
	tree->Branch("mc_truth_h0taunutau2_id", &mc_truth_h0taunutau2_id, "mc_truth_h0taunutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0taunu2_id", &mc_truth_h0taunu2_id, "mc_truth_h0taunu2_id/I", buffersize);

	tree->Branch("mc_truth_t1_id", &mc_truth_t1_id, "mc_truth_t1_id/I", buffersize);
	tree->Branch("mc_truth_t2_id", &mc_truth_t2_id, "mc_truth_t2_id/I", buffersize);
	tree->Branch("mc_truth_tb1_id", &mc_truth_tb1_id, "mc_truth_tb1_id/I", buffersize);
	tree->Branch("mc_truth_tb2_id", &mc_truth_tb2_id, "mc_truth_tb2_id/I", buffersize);

	tree->Branch("mc_truth_tW1_id", &mc_truth_tW1_id, "mc_truth_tW1_id/I", buffersize);
	tree->Branch("mc_truth_tWnu1_id", &mc_truth_tWnu1_id, "mc_truth_tWnu1_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_id", &mc_truth_tWnutau1_id, "mc_truth_tWnutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWl1_id", &mc_truth_tWl1_id, "mc_truth_tWl1_id/I", buffersize);
	tree->Branch("mc_truth_tWtau1_id", &mc_truth_tWtau1_id, "mc_truth_tWtau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_id", &mc_truth_tWtaunu1_id, "mc_truth_tWtaunu1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_id", &mc_truth_tWtaunutau1_id, "mc_truth_tWtaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_id", &mc_truth_tWtaul1_id, "mc_truth_tWtaul1_id/I", buffersize);
	tree->Branch("mc_truth_tWq11_id", &mc_truth_tWq11_id, "mc_truth_tWq11_id/I", buffersize);
	tree->Branch("mc_truth_tWq21_id", &mc_truth_tWq21_id, "mc_truth_tWq21_id/I", buffersize);

	tree->Branch("mc_truth_tW2_id", &mc_truth_tW2_id, "mc_truth_tW2_id/I", buffersize);
	tree->Branch("mc_truth_tWnu2_id", &mc_truth_tWnu2_id, "mc_truth_tWnu2_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_id", &mc_truth_tWnutau2_id, "mc_truth_tWnutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWl2_id", &mc_truth_tWl2_id, "mc_truth_tWl2_id/I", buffersize);
	tree->Branch("mc_truth_tWtau2_id", &mc_truth_tWtau2_id, "mc_truth_tWtau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_id", &mc_truth_tWtaunu2_id, "mc_truth_tWtaunu2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_id", &mc_truth_tWtaunutau2_id, "mc_truth_tWtaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_id", &mc_truth_tWtaul2_id, "mc_truth_tWtaul2_id/I", buffersize);
	tree->Branch("mc_truth_tWq12_id", &mc_truth_tWq12_id, "mc_truth_tWq12_id/I", buffersize);
	tree->Branch("mc_truth_tWq22_id", &mc_truth_tWq22_id, "mc_truth_tWq22_id/I", buffersize);

	tree->Branch("mc_truth_j1_id", &mc_truth_j1_id, "mc_truth_j1_id/I", buffersize);
	tree->Branch("mc_truth_j2_id", &mc_truth_j2_id, "mc_truth_j2_id/I", buffersize);
	tree->Branch("mc_truth_j3_id", &mc_truth_j3_id, "mc_truth_j3_id/I", buffersize);

	tree->Branch("mc_truth_h0_status", &mc_truth_h0_status, "mc_truth_h0_status/I", buffersize);

	tree->Branch("mc_truth_h0W1_status", &mc_truth_h0W1_status, "mc_truth_h0W1_status/I", buffersize);
	tree->Branch("mc_truth_h0W2_status", &mc_truth_h0W2_status, "mc_truth_h0W2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wl1_status", &mc_truth_h0Wl1_status, "mc_truth_h0Wl1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnu1_status", &mc_truth_h0Wnu1_status, "mc_truth_h0Wnu1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtau1_status", &mc_truth_h0Wtau1_status, "mc_truth_h0Wtau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_status", &mc_truth_h0Wnutau1_status, "mc_truth_h0Wnutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_status", &mc_truth_h0Wtaul1_status, "mc_truth_h0Wtaul1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_status", &mc_truth_h0Wtaunu1_status, "mc_truth_h0Wtaunu1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_status", &mc_truth_h0Wtaunutau1_status, "mc_truth_h0Wtaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wl2_status", &mc_truth_h0Wl2_status, "mc_truth_h0Wl2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnu2_status", &mc_truth_h0Wnu2_status, "mc_truth_h0Wnu2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtau2_status", &mc_truth_h0Wtau2_status, "mc_truth_h0Wtau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_status", &mc_truth_h0Wnutau2_status, "mc_truth_h0Wnutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_status", &mc_truth_h0Wtaul2_status, "mc_truth_h0Wtaul2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_status", &mc_truth_h0Wtaunu2_status, "mc_truth_h0Wtaunu2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_status", &mc_truth_h0Wtaunutau2_status, "mc_truth_h0Wtaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_status", &mc_truth_h0Wq11_status, "mc_truth_h0Wq11_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_status", &mc_truth_h0Wq21_status, "mc_truth_h0Wq21_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_status", &mc_truth_h0Wq12_status, "mc_truth_h0Wq12_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_status", &mc_truth_h0Wq22_status, "mc_truth_h0Wq22_status/I", buffersize);

	tree->Branch("mc_truth_h0Z1_status", &mc_truth_h0Z1_status, "mc_truth_h0Z1_status/I", buffersize);
	tree->Branch("mc_truth_h0Z2_status", &mc_truth_h0Z2_status, "mc_truth_h0Z2_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl11_status", &mc_truth_h0Zl11_status, "mc_truth_h0Zl11_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl21_status", &mc_truth_h0Zl21_status, "mc_truth_h0Zl21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau11_status", &mc_truth_h0Ztau11_status, "mc_truth_h0Ztau11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau21_status", &mc_truth_h0Ztau21_status, "mc_truth_h0Ztau21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_status", &mc_truth_h0Ztaul11_status, "mc_truth_h0Ztaul11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_status", &mc_truth_h0Ztaul21_status, "mc_truth_h0Ztaul21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_status", &mc_truth_h0Ztaunu11_status, "mc_truth_h0Ztaunu11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_status", &mc_truth_h0Ztaunu21_status, "mc_truth_h0Ztaunu21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_status", &mc_truth_h0Ztaunutau11_status, "mc_truth_h0Ztaunutau11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_status", &mc_truth_h0Ztaunutau21_status, "mc_truth_h0Ztaunutau21_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_status", &mc_truth_h0Zq11_status, "mc_truth_h0Zq11_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_status", &mc_truth_h0Zq21_status, "mc_truth_h0Zq21_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl12_status", &mc_truth_h0Zl12_status, "mc_truth_h0Zl12_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl22_status", &mc_truth_h0Zl22_status, "mc_truth_h0Zl22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau12_status", &mc_truth_h0Ztau12_status, "mc_truth_h0Ztau12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau22_status", &mc_truth_h0Ztau22_status, "mc_truth_h0Ztau22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_status", &mc_truth_h0Ztaul12_status, "mc_truth_h0Ztaul12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_status", &mc_truth_h0Ztaul22_status, "mc_truth_h0Ztaul22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_status", &mc_truth_h0Ztaunu12_status, "mc_truth_h0Ztaunu12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_status", &mc_truth_h0Ztaunu22_status, "mc_truth_h0Ztaunu22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_status", &mc_truth_h0Ztaunutau12_status, "mc_truth_h0Ztaunutau12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_status", &mc_truth_h0Ztaunutau22_status, "mc_truth_h0Ztaunutau22_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_status", &mc_truth_h0Zq12_status, "mc_truth_h0Zq12_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_status", &mc_truth_h0Zq22_status, "mc_truth_h0Zq22_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu11_status", &mc_truth_h0Znu11_status, "mc_truth_h0Znu11_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu21_status", &mc_truth_h0Znu21_status, "mc_truth_h0Znu21_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu12_status", &mc_truth_h0Znu12_status, "mc_truth_h0Znu12_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu22_status", &mc_truth_h0Znu22_status, "mc_truth_h0Znu22_status/I", buffersize);

	tree->Branch("mc_truth_h0tau1_status", &mc_truth_h0tau1_status, "mc_truth_h0tau1_status/I", buffersize);
	tree->Branch("mc_truth_h0tau2_status", &mc_truth_h0tau2_status, "mc_truth_h0tau2_status/I", buffersize);
	tree->Branch("mc_truth_h0taul1_status", &mc_truth_h0taul1_status, "mc_truth_h0taul1_status/I", buffersize);
	tree->Branch("mc_truth_h0taunutau1_status", &mc_truth_h0taunutau1_status, "mc_truth_h0taunutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0taunu1_status", &mc_truth_h0taunu1_status, "mc_truth_h0taunu1_status/I", buffersize);
	tree->Branch("mc_truth_h0taul2_status", &mc_truth_h0taul2_status, "mc_truth_h0taul2_status/I", buffersize);
	tree->Branch("mc_truth_h0taunutau2_status", &mc_truth_h0taunutau2_status, "mc_truth_h0taunutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0taunu2_status", &mc_truth_h0taunu2_status, "mc_truth_h0taunu2_status/I", buffersize);

	tree->Branch("mc_truth_t1_status", &mc_truth_t1_status, "mc_truth_t1_status/I", buffersize);
	tree->Branch("mc_truth_t2_status", &mc_truth_t2_status, "mc_truth_t2_status/I", buffersize);
	tree->Branch("mc_truth_tb1_status", &mc_truth_tb1_status, "mc_truth_tb1_status/I", buffersize);
	tree->Branch("mc_truth_tb2_status", &mc_truth_tb2_status, "mc_truth_tb2_status/I", buffersize);

	tree->Branch("mc_truth_tW1_status", &mc_truth_tW1_status, "mc_truth_tW1_status/I", buffersize);
	tree->Branch("mc_truth_tWnu1_status", &mc_truth_tWnu1_status, "mc_truth_tWnu1_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_status", &mc_truth_tWnutau1_status, "mc_truth_tWnutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWl1_status", &mc_truth_tWl1_status, "mc_truth_tWl1_status/I", buffersize);
	tree->Branch("mc_truth_tWtau1_status", &mc_truth_tWtau1_status, "mc_truth_tWtau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_status", &mc_truth_tWtaunu1_status, "mc_truth_tWtaunu1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_status", &mc_truth_tWtaunutau1_status, "mc_truth_tWtaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_status", &mc_truth_tWtaul1_status, "mc_truth_tWtaul1_status/I", buffersize);
	tree->Branch("mc_truth_tWq11_status", &mc_truth_tWq11_status, "mc_truth_tWq11_status/I", buffersize);
	tree->Branch("mc_truth_tWq21_status", &mc_truth_tWq21_status, "mc_truth_tWq21_status/I", buffersize);

	tree->Branch("mc_truth_tW2_status", &mc_truth_tW2_status, "mc_truth_tW2_status/I", buffersize);
	tree->Branch("mc_truth_tWnu2_status", &mc_truth_tWnu2_status, "mc_truth_tWnu2_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_status", &mc_truth_tWnutau2_status, "mc_truth_tWnutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWl2_status", &mc_truth_tWl2_status, "mc_truth_tWl2_status/I", buffersize);
	tree->Branch("mc_truth_tWtau2_status", &mc_truth_tWtau2_status, "mc_truth_tWtau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_status", &mc_truth_tWtaunu2_status, "mc_truth_tWtaunu2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_status", &mc_truth_tWtaunutau2_status, "mc_truth_tWtaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_status", &mc_truth_tWtaul2_status, "mc_truth_tWtaul2_status/I", buffersize);
	tree->Branch("mc_truth_tWq12_status", &mc_truth_tWq12_status, "mc_truth_tWq12_status/I", buffersize);
	tree->Branch("mc_truth_tWq22_status", &mc_truth_tWq22_status, "mc_truth_tWq22_status/I", buffersize);

	tree->Branch("mc_truth_j1_status", &mc_truth_j1_status, "mc_truth_j1_status/I", buffersize);
	tree->Branch("mc_truth_j2_status", &mc_truth_j2_status, "mc_truth_j2_status/I", buffersize);
	tree->Branch("mc_truth_j3_status", &mc_truth_j3_status, "mc_truth_j3_status/I", buffersize);
     }

   if( doWrite("mc_truth_ttz") )
     {
	if( doWrite("mc_truth_p4") )
	  {
	     tree->Branch("mc_truth_Z_p4", "TLorentzVector", &mc_truth_Z_p4, buffersize);
	     tree->Branch("mc_truth_Zl1_p4", "TLorentzVector", &mc_truth_Zl1_p4, buffersize);
	     tree->Branch("mc_truth_Zl2_p4", "TLorentzVector", &mc_truth_Zl2_p4, buffersize);
	     tree->Branch("mc_truth_Ztau1_p4", "TLorentzVector", &mc_truth_Ztau1_p4, buffersize);
	     tree->Branch("mc_truth_Ztau2_p4", "TLorentzVector", &mc_truth_Ztau2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaul1_p4", "TLorentzVector", &mc_truth_Ztaul1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaul2_p4", "TLorentzVector", &mc_truth_Ztaul2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunu1_p4", "TLorentzVector", &mc_truth_Ztaunu1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunu2_p4", "TLorentzVector", &mc_truth_Ztaunu2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunutau1_p4", "TLorentzVector", &mc_truth_Ztaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunutau2_p4", "TLorentzVector", &mc_truth_Ztaunutau2_p4, buffersize);	
	     tree->Branch("mc_truth_Zq1_p4", "TLorentzVector", &mc_truth_Zq1_p4, buffersize);
	     tree->Branch("mc_truth_Zq2_p4", "TLorentzVector", &mc_truth_Zq2_p4, buffersize);
	     tree->Branch("mc_truth_Znu1_p4", "TLorentzVector", &mc_truth_Znu1_p4, buffersize);
	     tree->Branch("mc_truth_Znu2_p4", "TLorentzVector", &mc_truth_Znu2_p4, buffersize);

	     tree->Branch("mc_truth_gammal1_p4", "TLorentzVector", &mc_truth_gammal1_p4, buffersize);
	     tree->Branch("mc_truth_gammal2_p4", "TLorentzVector", &mc_truth_gammal2_p4, buffersize);
	     tree->Branch("mc_truth_gammatau1_p4", "TLorentzVector", &mc_truth_gammatau1_p4, buffersize);
	     tree->Branch("mc_truth_gammatau2_p4", "TLorentzVector", &mc_truth_gammatau2_p4, buffersize);
	     tree->Branch("mc_truth_gammataul1_p4", "TLorentzVector", &mc_truth_gammataul1_p4, buffersize);
	     tree->Branch("mc_truth_gammataul2_p4", "TLorentzVector", &mc_truth_gammataul2_p4, buffersize);
	     tree->Branch("mc_truth_gammataunu1_p4", "TLorentzVector", &mc_truth_gammataunu1_p4, buffersize);
	     tree->Branch("mc_truth_gammataunu2_p4", "TLorentzVector", &mc_truth_gammataunu2_p4, buffersize);
	     tree->Branch("mc_truth_gammataunutau1_p4", "TLorentzVector", &mc_truth_gammataunutau1_p4, buffersize);
	     tree->Branch("mc_truth_gammataunutau2_p4", "TLorentzVector", &mc_truth_gammataunutau2_p4, buffersize);	
	     
	     tree->Branch("mc_truth_t1_p4", "TLorentzVector", &mc_truth_t1_p4, buffersize);
	     tree->Branch("mc_truth_t2_p4", "TLorentzVector", &mc_truth_t2_p4, buffersize);
	     tree->Branch("mc_truth_tb1_p4", "TLorentzVector", &mc_truth_tb1_p4, buffersize);
	     tree->Branch("mc_truth_tb2_p4", "TLorentzVector", &mc_truth_tb2_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW1_p4", "TLorentzVector", &mc_truth_tW1_p4, buffersize);
	     tree->Branch("mc_truth_tWnu1_p4", "TLorentzVector", &mc_truth_tWnu1_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau1_p4", "TLorentzVector", &mc_truth_tWnutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWl1_p4", "TLorentzVector", &mc_truth_tWl1_p4, buffersize);
	     tree->Branch("mc_truth_tWtau1_p4", "TLorentzVector", &mc_truth_tWtau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu1_p4", "TLorentzVector", &mc_truth_tWtaunu1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau1_p4", "TLorentzVector", &mc_truth_tWtaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul1_p4", "TLorentzVector", &mc_truth_tWtaul1_p4, buffersize);
	     tree->Branch("mc_truth_tWq11_p4", "TLorentzVector", &mc_truth_tWq11_p4, buffersize);
	     tree->Branch("mc_truth_tWq21_p4", "TLorentzVector", &mc_truth_tWq21_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW2_p4", "TLorentzVector", &mc_truth_tW2_p4, buffersize);
	     tree->Branch("mc_truth_tWnu2_p4", "TLorentzVector", &mc_truth_tWnu2_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau2_p4", "TLorentzVector", &mc_truth_tWnutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWl2_p4", "TLorentzVector", &mc_truth_tWl2_p4, buffersize);
	     tree->Branch("mc_truth_tWtau2_p4", "TLorentzVector", &mc_truth_tWtau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu2_p4", "TLorentzVector", &mc_truth_tWtaunu2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau2_p4", "TLorentzVector", &mc_truth_tWtaunutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul2_p4", "TLorentzVector", &mc_truth_tWtaul2_p4, buffersize);
	     tree->Branch("mc_truth_tWq12_p4", "TLorentzVector", &mc_truth_tWq12_p4, buffersize);
	     tree->Branch("mc_truth_tWq22_p4", "TLorentzVector", &mc_truth_tWq22_p4, buffersize);
	     
	     tree->Branch("mc_truth_j1_p4", "TLorentzVector", &mc_truth_j1_p4, buffersize);
	     tree->Branch("mc_truth_j2_p4", "TLorentzVector", &mc_truth_j2_p4, buffersize);
	     tree->Branch("mc_truth_j3_p4", "TLorentzVector", &mc_truth_j3_p4, buffersize);
	  }

	tree->Branch("mc_truth_Z_pt", &mc_truth_Z_pt, "mc_truth_Z_pt/F", buffersize);
	tree->Branch("mc_truth_Zl1_pt", &mc_truth_Zl1_pt, "mc_truth_Zl1_pt/F", buffersize);
	tree->Branch("mc_truth_Zl2_pt", &mc_truth_Zl2_pt, "mc_truth_Zl2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztau1_pt", &mc_truth_Ztau1_pt, "mc_truth_Ztau1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztau2_pt", &mc_truth_Ztau2_pt, "mc_truth_Ztau2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_pt", &mc_truth_Ztaul1_pt, "mc_truth_Ztaul1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_pt", &mc_truth_Ztaul2_pt, "mc_truth_Ztaul2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_pt", &mc_truth_Ztaunu1_pt, "mc_truth_Ztaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_pt", &mc_truth_Ztaunu2_pt, "mc_truth_Ztaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_pt", &mc_truth_Ztaunutau1_pt, "mc_truth_Ztaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_pt", &mc_truth_Ztaunutau2_pt, "mc_truth_Ztaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_Zq1_pt", &mc_truth_Zq1_pt, "mc_truth_Zq1_pt/F", buffersize);
	tree->Branch("mc_truth_Zq2_pt", &mc_truth_Zq2_pt, "mc_truth_Zq2_pt/F", buffersize);
	tree->Branch("mc_truth_Znu1_pt", &mc_truth_Znu1_pt, "mc_truth_Znu1_pt/F", buffersize);
	tree->Branch("mc_truth_Znu2_pt", &mc_truth_Znu2_pt, "mc_truth_Znu2_pt/F", buffersize);

	tree->Branch("mc_truth_gammal1_pt", &mc_truth_gammal1_pt, "mc_truth_gammal1_pt/F", buffersize);
	tree->Branch("mc_truth_gammal2_pt", &mc_truth_gammal2_pt, "mc_truth_gammal2_pt/F", buffersize);
	tree->Branch("mc_truth_gammatau1_pt", &mc_truth_gammatau1_pt, "mc_truth_gammatau1_pt/F", buffersize);
	tree->Branch("mc_truth_gammatau2_pt", &mc_truth_gammatau2_pt, "mc_truth_gammatau2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataul1_pt", &mc_truth_gammataul1_pt, "mc_truth_gammataul1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataul2_pt", &mc_truth_gammataul2_pt, "mc_truth_gammataul2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_pt", &mc_truth_gammataunu1_pt, "mc_truth_gammataunu1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_pt", &mc_truth_gammataunu2_pt, "mc_truth_gammataunu2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_pt", &mc_truth_gammataunutau1_pt, "mc_truth_gammataunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_pt", &mc_truth_gammataunutau2_pt, "mc_truth_gammataunutau2_pt/F", buffersize);
	
	tree->Branch("mc_truth_t1_pt", &mc_truth_t1_pt, "mc_truth_t1_pt/F", buffersize);
	tree->Branch("mc_truth_t2_pt", &mc_truth_t2_pt, "mc_truth_t2_pt/F", buffersize);
	tree->Branch("mc_truth_tb1_pt", &mc_truth_tb1_pt, "mc_truth_tb1_pt/F", buffersize);
	tree->Branch("mc_truth_tb2_pt", &mc_truth_tb2_pt, "mc_truth_tb2_pt/F", buffersize);

	tree->Branch("mc_truth_tW1_pt", &mc_truth_tW1_pt, "mc_truth_tW1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu1_pt", &mc_truth_tWnu1_pt, "mc_truth_tWnu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_pt", &mc_truth_tWnutau1_pt, "mc_truth_tWnutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWl1_pt", &mc_truth_tWl1_pt, "mc_truth_tWl1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau1_pt", &mc_truth_tWtau1_pt, "mc_truth_tWtau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_pt", &mc_truth_tWtaunu1_pt, "mc_truth_tWtaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_pt", &mc_truth_tWtaunutau1_pt, "mc_truth_tWtaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_pt", &mc_truth_tWtaul1_pt, "mc_truth_tWtaul1_pt/F", buffersize);
	tree->Branch("mc_truth_tWq11_pt", &mc_truth_tWq11_pt, "mc_truth_tWq11_pt/F", buffersize);
	tree->Branch("mc_truth_tWq21_pt", &mc_truth_tWq21_pt, "mc_truth_tWq21_pt/F", buffersize);

	tree->Branch("mc_truth_tW2_pt", &mc_truth_tW2_pt, "mc_truth_tW2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu2_pt", &mc_truth_tWnu2_pt, "mc_truth_tWnu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_pt", &mc_truth_tWnutau2_pt, "mc_truth_tWnutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWl2_pt", &mc_truth_tWl2_pt, "mc_truth_tWl2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau2_pt", &mc_truth_tWtau2_pt, "mc_truth_tWtau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_pt", &mc_truth_tWtaunu2_pt, "mc_truth_tWtaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_pt", &mc_truth_tWtaunutau2_pt, "mc_truth_tWtaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_pt", &mc_truth_tWtaul2_pt, "mc_truth_tWtaul2_pt/F", buffersize);
	tree->Branch("mc_truth_tWq12_pt", &mc_truth_tWq12_pt, "mc_truth_tWq12_pt/F", buffersize);
	tree->Branch("mc_truth_tWq22_pt", &mc_truth_tWq22_pt, "mc_truth_tWq22_pt/F", buffersize);

	tree->Branch("mc_truth_j1_pt", &mc_truth_j1_pt, "mc_truth_j1_pt/F", buffersize);
	tree->Branch("mc_truth_j2_pt", &mc_truth_j2_pt, "mc_truth_j2_pt/F", buffersize);
	tree->Branch("mc_truth_j3_pt", &mc_truth_j3_pt, "mc_truth_j3_pt/F", buffersize);
		
	tree->Branch("mc_truth_Z_eta", &mc_truth_Z_eta, "mc_truth_Z_eta/F", buffersize);
	tree->Branch("mc_truth_Zl1_eta", &mc_truth_Zl1_eta, "mc_truth_Zl1_eta/F", buffersize);
	tree->Branch("mc_truth_Zl2_eta", &mc_truth_Zl2_eta, "mc_truth_Zl2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztau1_eta", &mc_truth_Ztau1_eta, "mc_truth_Ztau1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztau2_eta", &mc_truth_Ztau2_eta, "mc_truth_Ztau2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_eta", &mc_truth_Ztaul1_eta, "mc_truth_Ztaul1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_eta", &mc_truth_Ztaul2_eta, "mc_truth_Ztaul2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_eta", &mc_truth_Ztaunu1_eta, "mc_truth_Ztaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_eta", &mc_truth_Ztaunu2_eta, "mc_truth_Ztaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_eta", &mc_truth_Ztaunutau1_eta, "mc_truth_Ztaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_eta", &mc_truth_Ztaunutau2_eta, "mc_truth_Ztaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_Zq1_eta", &mc_truth_Zq1_eta, "mc_truth_Zq1_eta/F", buffersize);
	tree->Branch("mc_truth_Zq2_eta", &mc_truth_Zq2_eta, "mc_truth_Zq2_eta/F", buffersize);
	tree->Branch("mc_truth_Znu1_eta", &mc_truth_Znu1_eta, "mc_truth_Znu1_eta/F", buffersize);
	tree->Branch("mc_truth_Znu2_eta", &mc_truth_Znu2_eta, "mc_truth_Znu2_eta/F", buffersize);

	tree->Branch("mc_truth_gammal1_eta", &mc_truth_gammal1_eta, "mc_truth_gammal1_eta/F", buffersize);
	tree->Branch("mc_truth_gammal2_eta", &mc_truth_gammal2_eta, "mc_truth_gammal2_eta/F", buffersize);
	tree->Branch("mc_truth_gammatau1_eta", &mc_truth_gammatau1_eta, "mc_truth_gammatau1_eta/F", buffersize);
	tree->Branch("mc_truth_gammatau2_eta", &mc_truth_gammatau2_eta, "mc_truth_gammatau2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataul1_eta", &mc_truth_gammataul1_eta, "mc_truth_gammataul1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataul2_eta", &mc_truth_gammataul2_eta, "mc_truth_gammataul2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_eta", &mc_truth_gammataunu1_eta, "mc_truth_gammataunu1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_eta", &mc_truth_gammataunu2_eta, "mc_truth_gammataunu2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_eta", &mc_truth_gammataunutau1_eta, "mc_truth_gammataunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_eta", &mc_truth_gammataunutau2_eta, "mc_truth_gammataunutau2_eta/F", buffersize);
	
	tree->Branch("mc_truth_t1_eta", &mc_truth_t1_eta, "mc_truth_t1_eta/F", buffersize);
	tree->Branch("mc_truth_t2_eta", &mc_truth_t2_eta, "mc_truth_t2_eta/F", buffersize);
	tree->Branch("mc_truth_tb1_eta", &mc_truth_tb1_eta, "mc_truth_tb1_eta/F", buffersize);
	tree->Branch("mc_truth_tb2_eta", &mc_truth_tb2_eta, "mc_truth_tb2_eta/F", buffersize);

	tree->Branch("mc_truth_tW1_eta", &mc_truth_tW1_eta, "mc_truth_tW1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu1_eta", &mc_truth_tWnu1_eta, "mc_truth_tWnu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_eta", &mc_truth_tWnutau1_eta, "mc_truth_tWnutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWl1_eta", &mc_truth_tWl1_eta, "mc_truth_tWl1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau1_eta", &mc_truth_tWtau1_eta, "mc_truth_tWtau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_eta", &mc_truth_tWtaunu1_eta, "mc_truth_tWtaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_eta", &mc_truth_tWtaunutau1_eta, "mc_truth_tWtaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_eta", &mc_truth_tWtaul1_eta, "mc_truth_tWtaul1_eta/F", buffersize);
	tree->Branch("mc_truth_tWq11_eta", &mc_truth_tWq11_eta, "mc_truth_tWq11_eta/F", buffersize);
	tree->Branch("mc_truth_tWq21_eta", &mc_truth_tWq21_eta, "mc_truth_tWq21_eta/F", buffersize);

	tree->Branch("mc_truth_tW2_eta", &mc_truth_tW2_eta, "mc_truth_tW2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu2_eta", &mc_truth_tWnu2_eta, "mc_truth_tWnu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_eta", &mc_truth_tWnutau2_eta, "mc_truth_tWnutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWl2_eta", &mc_truth_tWl2_eta, "mc_truth_tWl2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau2_eta", &mc_truth_tWtau2_eta, "mc_truth_tWtau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_eta", &mc_truth_tWtaunu2_eta, "mc_truth_tWtaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_eta", &mc_truth_tWtaunutau2_eta, "mc_truth_tWtaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_eta", &mc_truth_tWtaul2_eta, "mc_truth_tWtaul2_eta/F", buffersize);
	tree->Branch("mc_truth_tWq12_eta", &mc_truth_tWq12_eta, "mc_truth_tWq12_eta/F", buffersize);
	tree->Branch("mc_truth_tWq22_eta", &mc_truth_tWq22_eta, "mc_truth_tWq22_eta/F", buffersize);

	tree->Branch("mc_truth_j1_eta", &mc_truth_j1_eta, "mc_truth_j1_eta/F", buffersize);
	tree->Branch("mc_truth_j2_eta", &mc_truth_j2_eta, "mc_truth_j2_eta/F", buffersize);
	tree->Branch("mc_truth_j3_eta", &mc_truth_j3_eta, "mc_truth_j3_eta/F", buffersize);
		
	tree->Branch("mc_truth_Z_phi", &mc_truth_Z_phi, "mc_truth_Z_phi/F", buffersize);
	tree->Branch("mc_truth_Zl1_phi", &mc_truth_Zl1_phi, "mc_truth_Zl1_phi/F", buffersize);
	tree->Branch("mc_truth_Zl2_phi", &mc_truth_Zl2_phi, "mc_truth_Zl2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztau1_phi", &mc_truth_Ztau1_phi, "mc_truth_Ztau1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztau2_phi", &mc_truth_Ztau2_phi, "mc_truth_Ztau2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_phi", &mc_truth_Ztaul1_phi, "mc_truth_Ztaul1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_phi", &mc_truth_Ztaul2_phi, "mc_truth_Ztaul2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_phi", &mc_truth_Ztaunu1_phi, "mc_truth_Ztaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_phi", &mc_truth_Ztaunu2_phi, "mc_truth_Ztaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_phi", &mc_truth_Ztaunutau1_phi, "mc_truth_Ztaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_phi", &mc_truth_Ztaunutau2_phi, "mc_truth_Ztaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_Zq1_phi", &mc_truth_Zq1_phi, "mc_truth_Zq1_phi/F", buffersize);
	tree->Branch("mc_truth_Zq2_phi", &mc_truth_Zq2_phi, "mc_truth_Zq2_phi/F", buffersize);
	tree->Branch("mc_truth_Znu1_phi", &mc_truth_Znu1_phi, "mc_truth_Znu1_phi/F", buffersize);
	tree->Branch("mc_truth_Znu2_phi", &mc_truth_Znu2_phi, "mc_truth_Znu2_phi/F", buffersize);

	tree->Branch("mc_truth_gammal1_phi", &mc_truth_gammal1_phi, "mc_truth_gammal1_phi/F", buffersize);
	tree->Branch("mc_truth_gammal2_phi", &mc_truth_gammal2_phi, "mc_truth_gammal2_phi/F", buffersize);
	tree->Branch("mc_truth_gammatau1_phi", &mc_truth_gammatau1_phi, "mc_truth_gammatau1_phi/F", buffersize);
	tree->Branch("mc_truth_gammatau2_phi", &mc_truth_gammatau2_phi, "mc_truth_gammatau2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataul1_phi", &mc_truth_gammataul1_phi, "mc_truth_gammataul1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataul2_phi", &mc_truth_gammataul2_phi, "mc_truth_gammataul2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_phi", &mc_truth_gammataunu1_phi, "mc_truth_gammataunu1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_phi", &mc_truth_gammataunu2_phi, "mc_truth_gammataunu2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_phi", &mc_truth_gammataunutau1_phi, "mc_truth_gammataunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_phi", &mc_truth_gammataunutau2_phi, "mc_truth_gammataunutau2_phi/F", buffersize);
	
	tree->Branch("mc_truth_t1_phi", &mc_truth_t1_phi, "mc_truth_t1_phi/F", buffersize);
	tree->Branch("mc_truth_t2_phi", &mc_truth_t2_phi, "mc_truth_t2_phi/F", buffersize);
	tree->Branch("mc_truth_tb1_phi", &mc_truth_tb1_phi, "mc_truth_tb1_phi/F", buffersize);
	tree->Branch("mc_truth_tb2_phi", &mc_truth_tb2_phi, "mc_truth_tb2_phi/F", buffersize);

	tree->Branch("mc_truth_tW1_phi", &mc_truth_tW1_phi, "mc_truth_tW1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu1_phi", &mc_truth_tWnu1_phi, "mc_truth_tWnu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_phi", &mc_truth_tWnutau1_phi, "mc_truth_tWnutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWl1_phi", &mc_truth_tWl1_phi, "mc_truth_tWl1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau1_phi", &mc_truth_tWtau1_phi, "mc_truth_tWtau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_phi", &mc_truth_tWtaunu1_phi, "mc_truth_tWtaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_phi", &mc_truth_tWtaunutau1_phi, "mc_truth_tWtaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_phi", &mc_truth_tWtaul1_phi, "mc_truth_tWtaul1_phi/F", buffersize);
	tree->Branch("mc_truth_tWq11_phi", &mc_truth_tWq11_phi, "mc_truth_tWq11_phi/F", buffersize);
	tree->Branch("mc_truth_tWq21_phi", &mc_truth_tWq21_phi, "mc_truth_tWq21_phi/F", buffersize);

	tree->Branch("mc_truth_tW2_phi", &mc_truth_tW2_phi, "mc_truth_tW2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu2_phi", &mc_truth_tWnu2_phi, "mc_truth_tWnu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_phi", &mc_truth_tWnutau2_phi, "mc_truth_tWnutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWl2_phi", &mc_truth_tWl2_phi, "mc_truth_tWl2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau2_phi", &mc_truth_tWtau2_phi, "mc_truth_tWtau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_phi", &mc_truth_tWtaunu2_phi, "mc_truth_tWtaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_phi", &mc_truth_tWtaunutau2_phi, "mc_truth_tWtaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_phi", &mc_truth_tWtaul2_phi, "mc_truth_tWtaul2_phi/F", buffersize);
	tree->Branch("mc_truth_tWq12_phi", &mc_truth_tWq12_phi, "mc_truth_tWq12_phi/F", buffersize);
	tree->Branch("mc_truth_tWq22_phi", &mc_truth_tWq22_phi, "mc_truth_tWq22_phi/F", buffersize);

	tree->Branch("mc_truth_j1_phi", &mc_truth_j1_phi, "mc_truth_j1_phi/F", buffersize);
	tree->Branch("mc_truth_j2_phi", &mc_truth_j2_phi, "mc_truth_j2_phi/F", buffersize);
	tree->Branch("mc_truth_j3_phi", &mc_truth_j3_phi, "mc_truth_j3_phi/F", buffersize);
	
	tree->Branch("mc_truth_Z_E", &mc_truth_Z_E, "mc_truth_Z_E/F", buffersize);
	tree->Branch("mc_truth_Zl1_E", &mc_truth_Zl1_E, "mc_truth_Zl1_E/F", buffersize);
	tree->Branch("mc_truth_Zl2_E", &mc_truth_Zl2_E, "mc_truth_Zl2_E/F", buffersize);
	tree->Branch("mc_truth_Ztau1_E", &mc_truth_Ztau1_E, "mc_truth_Ztau1_E/F", buffersize);
	tree->Branch("mc_truth_Ztau2_E", &mc_truth_Ztau2_E, "mc_truth_Ztau2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_E", &mc_truth_Ztaul1_E, "mc_truth_Ztaul1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_E", &mc_truth_Ztaul2_E, "mc_truth_Ztaul2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_E", &mc_truth_Ztaunu1_E, "mc_truth_Ztaunu1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_E", &mc_truth_Ztaunu2_E, "mc_truth_Ztaunu2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_E", &mc_truth_Ztaunutau1_E, "mc_truth_Ztaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_E", &mc_truth_Ztaunutau2_E, "mc_truth_Ztaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_Zq1_E", &mc_truth_Zq1_E, "mc_truth_Zq1_E/F", buffersize);
	tree->Branch("mc_truth_Zq2_E", &mc_truth_Zq2_E, "mc_truth_Zq2_E/F", buffersize);
	tree->Branch("mc_truth_Znu1_E", &mc_truth_Znu1_E, "mc_truth_Znu1_E/F", buffersize);
	tree->Branch("mc_truth_Znu2_E", &mc_truth_Znu2_E, "mc_truth_Znu2_E/F", buffersize);

	tree->Branch("mc_truth_gammal1_E", &mc_truth_gammal1_E, "mc_truth_gammal1_E/F", buffersize);
	tree->Branch("mc_truth_gammal2_E", &mc_truth_gammal2_E, "mc_truth_gammal2_E/F", buffersize);
	tree->Branch("mc_truth_gammatau1_E", &mc_truth_gammatau1_E, "mc_truth_gammatau1_E/F", buffersize);
	tree->Branch("mc_truth_gammatau2_E", &mc_truth_gammatau2_E, "mc_truth_gammatau2_E/F", buffersize);
	tree->Branch("mc_truth_gammataul1_E", &mc_truth_gammataul1_E, "mc_truth_gammataul1_E/F", buffersize);
	tree->Branch("mc_truth_gammataul2_E", &mc_truth_gammataul2_E, "mc_truth_gammataul2_E/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_E", &mc_truth_gammataunu1_E, "mc_truth_gammataunu1_E/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_E", &mc_truth_gammataunu2_E, "mc_truth_gammataunu2_E/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_E", &mc_truth_gammataunutau1_E, "mc_truth_gammataunutau1_E/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_E", &mc_truth_gammataunutau2_E, "mc_truth_gammataunutau2_E/F", buffersize);
	
	tree->Branch("mc_truth_t1_E", &mc_truth_t1_E, "mc_truth_t1_E/F", buffersize);
	tree->Branch("mc_truth_t2_E", &mc_truth_t2_E, "mc_truth_t2_E/F", buffersize);
	tree->Branch("mc_truth_tb1_E", &mc_truth_tb1_E, "mc_truth_tb1_E/F", buffersize);
	tree->Branch("mc_truth_tb2_E", &mc_truth_tb2_E, "mc_truth_tb2_E/F", buffersize);

	tree->Branch("mc_truth_tW1_E", &mc_truth_tW1_E, "mc_truth_tW1_E/F", buffersize);
	tree->Branch("mc_truth_tWnu1_E", &mc_truth_tWnu1_E, "mc_truth_tWnu1_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_E", &mc_truth_tWnutau1_E, "mc_truth_tWnutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWl1_E", &mc_truth_tWl1_E, "mc_truth_tWl1_E/F", buffersize);
	tree->Branch("mc_truth_tWtau1_E", &mc_truth_tWtau1_E, "mc_truth_tWtau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_E", &mc_truth_tWtaunu1_E, "mc_truth_tWtaunu1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_E", &mc_truth_tWtaunutau1_E, "mc_truth_tWtaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_E", &mc_truth_tWtaul1_E, "mc_truth_tWtaul1_E/F", buffersize);
	tree->Branch("mc_truth_tWq11_E", &mc_truth_tWq11_E, "mc_truth_tWq11_E/F", buffersize);
	tree->Branch("mc_truth_tWq21_E", &mc_truth_tWq21_E, "mc_truth_tWq21_E/F", buffersize);

	tree->Branch("mc_truth_tW2_E", &mc_truth_tW2_E, "mc_truth_tW2_E/F", buffersize);
	tree->Branch("mc_truth_tWnu2_E", &mc_truth_tWnu2_E, "mc_truth_tWnu2_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_E", &mc_truth_tWnutau2_E, "mc_truth_tWnutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWl2_E", &mc_truth_tWl2_E, "mc_truth_tWl2_E/F", buffersize);
	tree->Branch("mc_truth_tWtau2_E", &mc_truth_tWtau2_E, "mc_truth_tWtau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_E", &mc_truth_tWtaunu2_E, "mc_truth_tWtaunu2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_E", &mc_truth_tWtaunutau2_E, "mc_truth_tWtaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_E", &mc_truth_tWtaul2_E, "mc_truth_tWtaul2_E/F", buffersize);
	tree->Branch("mc_truth_tWq12_E", &mc_truth_tWq12_E, "mc_truth_tWq12_E/F", buffersize);
	tree->Branch("mc_truth_tWq22_E", &mc_truth_tWq22_E, "mc_truth_tWq22_E/F", buffersize);

	tree->Branch("mc_truth_j1_E", &mc_truth_j1_E, "mc_truth_j1_E/F", buffersize);
	tree->Branch("mc_truth_j2_E", &mc_truth_j2_E, "mc_truth_j2_E/F", buffersize);
	tree->Branch("mc_truth_j3_E", &mc_truth_j3_E, "mc_truth_j3_E/F", buffersize);
		
	tree->Branch("mc_truth_Z_id", &mc_truth_Z_id, "mc_truth_Z_id/I", buffersize);
	tree->Branch("mc_truth_Zl1_id", &mc_truth_Zl1_id, "mc_truth_Zl1_id/I", buffersize);
	tree->Branch("mc_truth_Zl2_id", &mc_truth_Zl2_id, "mc_truth_Zl2_id/I", buffersize);
	tree->Branch("mc_truth_Ztau1_id", &mc_truth_Ztau1_id, "mc_truth_Ztau1_id/I", buffersize);
	tree->Branch("mc_truth_Ztau2_id", &mc_truth_Ztau2_id, "mc_truth_Ztau2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaul1_id", &mc_truth_Ztaul1_id, "mc_truth_Ztaul1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaul2_id", &mc_truth_Ztaul2_id, "mc_truth_Ztaul2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunu1_id", &mc_truth_Ztaunu1_id, "mc_truth_Ztaunu1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunu2_id", &mc_truth_Ztaunu2_id, "mc_truth_Ztaunu2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_id", &mc_truth_Ztaunutau1_id, "mc_truth_Ztaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_id", &mc_truth_Ztaunutau2_id, "mc_truth_Ztaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_Zq1_id", &mc_truth_Zq1_id, "mc_truth_Zq1_id/I", buffersize);
	tree->Branch("mc_truth_Zq2_id", &mc_truth_Zq2_id, "mc_truth_Zq2_id/I", buffersize);
	tree->Branch("mc_truth_Znu1_id", &mc_truth_Znu1_id, "mc_truth_Znu1_id/I", buffersize);
	tree->Branch("mc_truth_Znu2_id", &mc_truth_Znu2_id, "mc_truth_Znu2_id/I", buffersize);

	tree->Branch("mc_truth_gammal1_id", &mc_truth_gammal1_id, "mc_truth_gammal1_id/I", buffersize);
	tree->Branch("mc_truth_gammal2_id", &mc_truth_gammal2_id, "mc_truth_gammal2_id/I", buffersize);
	tree->Branch("mc_truth_gammatau1_id", &mc_truth_gammatau1_id, "mc_truth_gammatau1_id/I", buffersize);
	tree->Branch("mc_truth_gammatau2_id", &mc_truth_gammatau2_id, "mc_truth_gammatau2_id/I", buffersize);
	tree->Branch("mc_truth_gammataul1_id", &mc_truth_gammataul1_id, "mc_truth_gammataul1_id/I", buffersize);
	tree->Branch("mc_truth_gammataul2_id", &mc_truth_gammataul2_id, "mc_truth_gammataul2_id/I", buffersize);
	tree->Branch("mc_truth_gammataunu1_id", &mc_truth_gammataunu1_id, "mc_truth_gammataunu1_id/I", buffersize);
	tree->Branch("mc_truth_gammataunu2_id", &mc_truth_gammataunu2_id, "mc_truth_gammataunu2_id/I", buffersize);
	tree->Branch("mc_truth_gammataunutau1_id", &mc_truth_gammataunutau1_id, "mc_truth_gammataunutau1_id/I", buffersize);
	tree->Branch("mc_truth_gammataunutau2_id", &mc_truth_gammataunutau2_id, "mc_truth_gammataunutau2_id/I", buffersize);
	
	tree->Branch("mc_truth_t1_id", &mc_truth_t1_id, "mc_truth_t1_id/I", buffersize);
	tree->Branch("mc_truth_t2_id", &mc_truth_t2_id, "mc_truth_t2_id/I", buffersize);
	tree->Branch("mc_truth_tb1_id", &mc_truth_tb1_id, "mc_truth_tb1_id/I", buffersize);
	tree->Branch("mc_truth_tb2_id", &mc_truth_tb2_id, "mc_truth_tb2_id/I", buffersize);

	tree->Branch("mc_truth_tW1_id", &mc_truth_tW1_id, "mc_truth_tW1_id/I", buffersize);
	tree->Branch("mc_truth_tWnu1_id", &mc_truth_tWnu1_id, "mc_truth_tWnu1_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_id", &mc_truth_tWnutau1_id, "mc_truth_tWnutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWl1_id", &mc_truth_tWl1_id, "mc_truth_tWl1_id/I", buffersize);
	tree->Branch("mc_truth_tWtau1_id", &mc_truth_tWtau1_id, "mc_truth_tWtau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_id", &mc_truth_tWtaunu1_id, "mc_truth_tWtaunu1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_id", &mc_truth_tWtaunutau1_id, "mc_truth_tWtaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_id", &mc_truth_tWtaul1_id, "mc_truth_tWtaul1_id/I", buffersize);
	tree->Branch("mc_truth_tWq11_id", &mc_truth_tWq11_id, "mc_truth_tWq11_id/I", buffersize);
	tree->Branch("mc_truth_tWq21_id", &mc_truth_tWq21_id, "mc_truth_tWq21_id/I", buffersize);

	tree->Branch("mc_truth_tW2_id", &mc_truth_tW2_id, "mc_truth_tW2_id/I", buffersize);
	tree->Branch("mc_truth_tWnu2_id", &mc_truth_tWnu2_id, "mc_truth_tWnu2_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_id", &mc_truth_tWnutau2_id, "mc_truth_tWnutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWl2_id", &mc_truth_tWl2_id, "mc_truth_tWl2_id/I", buffersize);
	tree->Branch("mc_truth_tWtau2_id", &mc_truth_tWtau2_id, "mc_truth_tWtau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_id", &mc_truth_tWtaunu2_id, "mc_truth_tWtaunu2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_id", &mc_truth_tWtaunutau2_id, "mc_truth_tWtaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_id", &mc_truth_tWtaul2_id, "mc_truth_tWtaul2_id/I", buffersize);
	tree->Branch("mc_truth_tWq12_id", &mc_truth_tWq12_id, "mc_truth_tWq12_id/I", buffersize);
	tree->Branch("mc_truth_tWq22_id", &mc_truth_tWq22_id, "mc_truth_tWq22_id/I", buffersize);

	tree->Branch("mc_truth_j1_id", &mc_truth_j1_id, "mc_truth_j1_id/I", buffersize);
	tree->Branch("mc_truth_j2_id", &mc_truth_j2_id, "mc_truth_j2_id/I", buffersize);
	tree->Branch("mc_truth_j3_id", &mc_truth_j3_id, "mc_truth_j3_id/I", buffersize);

	tree->Branch("mc_truth_Z_status", &mc_truth_Z_status, "mc_truth_Z_status/I", buffersize);
	tree->Branch("mc_truth_Zl1_status", &mc_truth_Zl1_status, "mc_truth_Zl1_status/I", buffersize);
	tree->Branch("mc_truth_Zl2_status", &mc_truth_Zl2_status, "mc_truth_Zl2_status/I", buffersize);
	tree->Branch("mc_truth_Ztau1_status", &mc_truth_Ztau1_status, "mc_truth_Ztau1_status/I", buffersize);
	tree->Branch("mc_truth_Ztau2_status", &mc_truth_Ztau2_status, "mc_truth_Ztau2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaul1_status", &mc_truth_Ztaul1_status, "mc_truth_Ztaul1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaul2_status", &mc_truth_Ztaul2_status, "mc_truth_Ztaul2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunu1_status", &mc_truth_Ztaunu1_status, "mc_truth_Ztaunu1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunu2_status", &mc_truth_Ztaunu2_status, "mc_truth_Ztaunu2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_status", &mc_truth_Ztaunutau1_status, "mc_truth_Ztaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_status", &mc_truth_Ztaunutau2_status, "mc_truth_Ztaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_Zq1_status", &mc_truth_Zq1_status, "mc_truth_Zq1_status/I", buffersize);
	tree->Branch("mc_truth_Zq2_status", &mc_truth_Zq2_status, "mc_truth_Zq2_status/I", buffersize);
	tree->Branch("mc_truth_Znu1_status", &mc_truth_Znu1_status, "mc_truth_Znu1_status/I", buffersize);
	tree->Branch("mc_truth_Znu2_status", &mc_truth_Znu2_status, "mc_truth_Znu2_status/I", buffersize);

	tree->Branch("mc_truth_gammal1_status", &mc_truth_gammal1_status, "mc_truth_gammal1_status/I", buffersize);
	tree->Branch("mc_truth_gammal2_status", &mc_truth_gammal2_status, "mc_truth_gammal2_status/I", buffersize);
	tree->Branch("mc_truth_gammatau1_status", &mc_truth_gammatau1_status, "mc_truth_gammatau1_status/I", buffersize);
	tree->Branch("mc_truth_gammatau2_status", &mc_truth_gammatau2_status, "mc_truth_gammatau2_status/I", buffersize);
	tree->Branch("mc_truth_gammataul1_status", &mc_truth_gammataul1_status, "mc_truth_gammataul1_status/I", buffersize);
	tree->Branch("mc_truth_gammataul2_status", &mc_truth_gammataul2_status, "mc_truth_gammataul2_status/I", buffersize);
	tree->Branch("mc_truth_gammataunu1_status", &mc_truth_gammataunu1_status, "mc_truth_gammataunu1_status/I", buffersize);
	tree->Branch("mc_truth_gammataunu2_status", &mc_truth_gammataunu2_status, "mc_truth_gammataunu2_status/I", buffersize);
	tree->Branch("mc_truth_gammataunutau1_status", &mc_truth_gammataunutau1_status, "mc_truth_gammataunutau1_status/I", buffersize);
	tree->Branch("mc_truth_gammataunutau2_status", &mc_truth_gammataunutau2_status, "mc_truth_gammataunutau2_status/I", buffersize);
	
	tree->Branch("mc_truth_t1_status", &mc_truth_t1_status, "mc_truth_t1_status/I", buffersize);
	tree->Branch("mc_truth_t2_status", &mc_truth_t2_status, "mc_truth_t2_status/I", buffersize);
	tree->Branch("mc_truth_tb1_status", &mc_truth_tb1_status, "mc_truth_tb1_status/I", buffersize);
	tree->Branch("mc_truth_tb2_status", &mc_truth_tb2_status, "mc_truth_tb2_status/I", buffersize);

	tree->Branch("mc_truth_tW1_status", &mc_truth_tW1_status, "mc_truth_tW1_status/I", buffersize);
	tree->Branch("mc_truth_tWnu1_status", &mc_truth_tWnu1_status, "mc_truth_tWnu1_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_status", &mc_truth_tWnutau1_status, "mc_truth_tWnutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWl1_status", &mc_truth_tWl1_status, "mc_truth_tWl1_status/I", buffersize);
	tree->Branch("mc_truth_tWtau1_status", &mc_truth_tWtau1_status, "mc_truth_tWtau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_status", &mc_truth_tWtaunu1_status, "mc_truth_tWtaunu1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_status", &mc_truth_tWtaunutau1_status, "mc_truth_tWtaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_status", &mc_truth_tWtaul1_status, "mc_truth_tWtaul1_status/I", buffersize);
	tree->Branch("mc_truth_tWq11_status", &mc_truth_tWq11_status, "mc_truth_tWq11_status/I", buffersize);
	tree->Branch("mc_truth_tWq21_status", &mc_truth_tWq21_status, "mc_truth_tWq21_status/I", buffersize);

	tree->Branch("mc_truth_tW2_status", &mc_truth_tW2_status, "mc_truth_tW2_status/I", buffersize);
	tree->Branch("mc_truth_tWnu2_status", &mc_truth_tWnu2_status, "mc_truth_tWnu2_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_status", &mc_truth_tWnutau2_status, "mc_truth_tWnutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWl2_status", &mc_truth_tWl2_status, "mc_truth_tWl2_status/I", buffersize);
	tree->Branch("mc_truth_tWtau2_status", &mc_truth_tWtau2_status, "mc_truth_tWtau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_status", &mc_truth_tWtaunu2_status, "mc_truth_tWtaunu2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_status", &mc_truth_tWtaunutau2_status, "mc_truth_tWtaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_status", &mc_truth_tWtaul2_status, "mc_truth_tWtaul2_status/I", buffersize);
	tree->Branch("mc_truth_tWq12_status", &mc_truth_tWq12_status, "mc_truth_tWq12_status/I", buffersize);
	tree->Branch("mc_truth_tWq22_status", &mc_truth_tWq22_status, "mc_truth_tWq22_status/I", buffersize);

	tree->Branch("mc_truth_j1_status", &mc_truth_j1_status, "mc_truth_j1_status/I", buffersize);
	tree->Branch("mc_truth_j2_status", &mc_truth_j2_status, "mc_truth_j2_status/I", buffersize);
	tree->Branch("mc_truth_j3_status", &mc_truth_j3_status, "mc_truth_j3_status/I", buffersize);
     }

   if( doWrite("mc_truth_ttw") )
     {
	if( doWrite("mc_truth_p4") )
	  {
	     tree->Branch("mc_truth_W_p4", "TLorentzVector", &mc_truth_W_p4, buffersize);
	     tree->Branch("mc_truth_Wnu_p4", "TLorentzVector", &mc_truth_Wnu_p4, buffersize);
	     tree->Branch("mc_truth_Wnutau_p4", "TLorentzVector", &mc_truth_Wnutau_p4, buffersize);
	     tree->Branch("mc_truth_Wl_p4", "TLorentzVector", &mc_truth_Wl_p4, buffersize);
	     tree->Branch("mc_truth_Wtau_p4", "TLorentzVector", &mc_truth_Wtau_p4, buffersize);
	     tree->Branch("mc_truth_Wtaunu_p4", "TLorentzVector", &mc_truth_Wtaunu_p4, buffersize);
	     tree->Branch("mc_truth_Wtaunutau_p4", "TLorentzVector", &mc_truth_Wtaunutau_p4, buffersize);
	     tree->Branch("mc_truth_Wtaul_p4", "TLorentzVector", &mc_truth_Wtaul_p4, buffersize);
	     tree->Branch("mc_truth_Wq1_p4", "TLorentzVector", &mc_truth_Wq1_p4, buffersize);
	     tree->Branch("mc_truth_Wq2_p4", "TLorentzVector", &mc_truth_Wq2_p4, buffersize);
	     
	     tree->Branch("mc_truth_gammal1_p4", "TLorentzVector", &mc_truth_gammal1_p4, buffersize);
	     tree->Branch("mc_truth_gammal2_p4", "TLorentzVector", &mc_truth_gammal2_p4, buffersize);
	     tree->Branch("mc_truth_gammatau1_p4", "TLorentzVector", &mc_truth_gammatau1_p4, buffersize);
	     tree->Branch("mc_truth_gammatau2_p4", "TLorentzVector", &mc_truth_gammatau2_p4, buffersize);
	     tree->Branch("mc_truth_gammataul1_p4", "TLorentzVector", &mc_truth_gammataul1_p4, buffersize);
	     tree->Branch("mc_truth_gammataul2_p4", "TLorentzVector", &mc_truth_gammataul2_p4, buffersize);
	     tree->Branch("mc_truth_gammataunu1_p4", "TLorentzVector", &mc_truth_gammataunu1_p4, buffersize);
	     tree->Branch("mc_truth_gammataunu2_p4", "TLorentzVector", &mc_truth_gammataunu2_p4, buffersize);
	     tree->Branch("mc_truth_gammataunutau1_p4", "TLorentzVector", &mc_truth_gammataunutau1_p4, buffersize);
	     tree->Branch("mc_truth_gammataunutau2_p4", "TLorentzVector", &mc_truth_gammataunutau2_p4, buffersize);	
	     
	     tree->Branch("mc_truth_t1_p4", "TLorentzVector", &mc_truth_t1_p4, buffersize);
	     tree->Branch("mc_truth_t2_p4", "TLorentzVector", &mc_truth_t2_p4, buffersize);
	     tree->Branch("mc_truth_tb1_p4", "TLorentzVector", &mc_truth_tb1_p4, buffersize);
	     tree->Branch("mc_truth_tb2_p4", "TLorentzVector", &mc_truth_tb2_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW1_p4", "TLorentzVector", &mc_truth_tW1_p4, buffersize);
	     tree->Branch("mc_truth_tWnu1_p4", "TLorentzVector", &mc_truth_tWnu1_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau1_p4", "TLorentzVector", &mc_truth_tWnutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWl1_p4", "TLorentzVector", &mc_truth_tWl1_p4, buffersize);
	     tree->Branch("mc_truth_tWtau1_p4", "TLorentzVector", &mc_truth_tWtau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu1_p4", "TLorentzVector", &mc_truth_tWtaunu1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau1_p4", "TLorentzVector", &mc_truth_tWtaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul1_p4", "TLorentzVector", &mc_truth_tWtaul1_p4, buffersize);
	     tree->Branch("mc_truth_tWq11_p4", "TLorentzVector", &mc_truth_tWq11_p4, buffersize);
	     tree->Branch("mc_truth_tWq21_p4", "TLorentzVector", &mc_truth_tWq21_p4, buffersize);

	     tree->Branch("mc_truth_tW2_p4", "TLorentzVector", &mc_truth_tW2_p4, buffersize);
	     tree->Branch("mc_truth_tWnu2_p4", "TLorentzVector", &mc_truth_tWnu2_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau2_p4", "TLorentzVector", &mc_truth_tWnutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWl2_p4", "TLorentzVector", &mc_truth_tWl2_p4, buffersize);
	     tree->Branch("mc_truth_tWtau2_p4", "TLorentzVector", &mc_truth_tWtau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu2_p4", "TLorentzVector", &mc_truth_tWtaunu2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau2_p4", "TLorentzVector", &mc_truth_tWtaunutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul2_p4", "TLorentzVector", &mc_truth_tWtaul2_p4, buffersize);
	     tree->Branch("mc_truth_tWq12_p4", "TLorentzVector", &mc_truth_tWq12_p4, buffersize);
	     tree->Branch("mc_truth_tWq22_p4", "TLorentzVector", &mc_truth_tWq22_p4, buffersize);
	     
	     tree->Branch("mc_truth_j1_p4", "TLorentzVector", &mc_truth_j1_p4, buffersize);
	     tree->Branch("mc_truth_j2_p4", "TLorentzVector", &mc_truth_j2_p4, buffersize);
	     tree->Branch("mc_truth_j3_p4", "TLorentzVector", &mc_truth_j3_p4, buffersize);
	  }	

	tree->Branch("mc_truth_W_pt", &mc_truth_W_pt, "mc_truth_W_pt/F", buffersize);
	tree->Branch("mc_truth_Wnu_pt", &mc_truth_Wnu_pt, "mc_truth_Wnu_pt/F", buffersize);
	tree->Branch("mc_truth_Wnutau_pt", &mc_truth_Wnutau_pt, "mc_truth_Wnutau_pt/F", buffersize);
	tree->Branch("mc_truth_Wl_pt", &mc_truth_Wl_pt, "mc_truth_Wl_pt/F", buffersize);
	tree->Branch("mc_truth_Wtau_pt", &mc_truth_Wtau_pt, "mc_truth_Wtau_pt/F", buffersize);
	tree->Branch("mc_truth_Wtaunu_pt", &mc_truth_Wtaunu_pt, "mc_truth_Wtaunu_pt/F", buffersize);
	tree->Branch("mc_truth_Wtaunutau_pt", &mc_truth_Wtaunutau_pt, "mc_truth_Wtaunutau_pt/F", buffersize);
	tree->Branch("mc_truth_Wtaul_pt", &mc_truth_Wtaul_pt, "mc_truth_Wtaul_pt/F", buffersize);
	tree->Branch("mc_truth_Wq1_pt", &mc_truth_Wq1_pt, "mc_truth_Wq1_pt/F", buffersize);
	tree->Branch("mc_truth_Wq2_pt", &mc_truth_Wq2_pt, "mc_truth_Wq2_pt/F", buffersize);

	tree->Branch("mc_truth_gammal1_pt", &mc_truth_gammal1_pt, "mc_truth_gammal1_pt/F", buffersize);
	tree->Branch("mc_truth_gammal2_pt", &mc_truth_gammal2_pt, "mc_truth_gammal2_pt/F", buffersize);
	tree->Branch("mc_truth_gammatau1_pt", &mc_truth_gammatau1_pt, "mc_truth_gammatau1_pt/F", buffersize);
	tree->Branch("mc_truth_gammatau2_pt", &mc_truth_gammatau2_pt, "mc_truth_gammatau2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataul1_pt", &mc_truth_gammataul1_pt, "mc_truth_gammataul1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataul2_pt", &mc_truth_gammataul2_pt, "mc_truth_gammataul2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_pt", &mc_truth_gammataunu1_pt, "mc_truth_gammataunu1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_pt", &mc_truth_gammataunu2_pt, "mc_truth_gammataunu2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_pt", &mc_truth_gammataunutau1_pt, "mc_truth_gammataunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_pt", &mc_truth_gammataunutau2_pt, "mc_truth_gammataunutau2_pt/F", buffersize);
	
	tree->Branch("mc_truth_t1_pt", &mc_truth_t1_pt, "mc_truth_t1_pt/F", buffersize);
	tree->Branch("mc_truth_t2_pt", &mc_truth_t2_pt, "mc_truth_t2_pt/F", buffersize);
	tree->Branch("mc_truth_tb1_pt", &mc_truth_tb1_pt, "mc_truth_tb1_pt/F", buffersize);
	tree->Branch("mc_truth_tb2_pt", &mc_truth_tb2_pt, "mc_truth_tb2_pt/F", buffersize);

	tree->Branch("mc_truth_tW1_pt", &mc_truth_tW1_pt, "mc_truth_tW1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu1_pt", &mc_truth_tWnu1_pt, "mc_truth_tWnu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_pt", &mc_truth_tWnutau1_pt, "mc_truth_tWnutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWl1_pt", &mc_truth_tWl1_pt, "mc_truth_tWl1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau1_pt", &mc_truth_tWtau1_pt, "mc_truth_tWtau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_pt", &mc_truth_tWtaunu1_pt, "mc_truth_tWtaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_pt", &mc_truth_tWtaunutau1_pt, "mc_truth_tWtaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_pt", &mc_truth_tWtaul1_pt, "mc_truth_tWtaul1_pt/F", buffersize);
	tree->Branch("mc_truth_tWq11_pt", &mc_truth_tWq11_pt, "mc_truth_tWq11_pt/F", buffersize);
	tree->Branch("mc_truth_tWq21_pt", &mc_truth_tWq21_pt, "mc_truth_tWq21_pt/F", buffersize);

	tree->Branch("mc_truth_tW2_pt", &mc_truth_tW2_pt, "mc_truth_tW2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu2_pt", &mc_truth_tWnu2_pt, "mc_truth_tWnu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_pt", &mc_truth_tWnutau2_pt, "mc_truth_tWnutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWl2_pt", &mc_truth_tWl2_pt, "mc_truth_tWl2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau2_pt", &mc_truth_tWtau2_pt, "mc_truth_tWtau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_pt", &mc_truth_tWtaunu2_pt, "mc_truth_tWtaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_pt", &mc_truth_tWtaunutau2_pt, "mc_truth_tWtaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_pt", &mc_truth_tWtaul2_pt, "mc_truth_tWtaul2_pt/F", buffersize);
	tree->Branch("mc_truth_tWq12_pt", &mc_truth_tWq12_pt, "mc_truth_tWq12_pt/F", buffersize);
	tree->Branch("mc_truth_tWq22_pt", &mc_truth_tWq22_pt, "mc_truth_tWq22_pt/F", buffersize);

	tree->Branch("mc_truth_j1_pt", &mc_truth_j1_pt, "mc_truth_j1_pt/F", buffersize);
	tree->Branch("mc_truth_j2_pt", &mc_truth_j2_pt, "mc_truth_j2_pt/F", buffersize);
	tree->Branch("mc_truth_j3_pt", &mc_truth_j3_pt, "mc_truth_j3_pt/F", buffersize);
			
	tree->Branch("mc_truth_W_eta", &mc_truth_W_eta, "mc_truth_W_eta/F", buffersize);
	tree->Branch("mc_truth_Wnu_eta", &mc_truth_Wnu_eta, "mc_truth_Wnu_eta/F", buffersize);
	tree->Branch("mc_truth_Wnutau_eta", &mc_truth_Wnutau_eta, "mc_truth_Wnutau_eta/F", buffersize);
	tree->Branch("mc_truth_Wl_eta", &mc_truth_Wl_eta, "mc_truth_Wl_eta/F", buffersize);
	tree->Branch("mc_truth_Wtau_eta", &mc_truth_Wtau_eta, "mc_truth_Wtau_eta/F", buffersize);
	tree->Branch("mc_truth_Wtaunu_eta", &mc_truth_Wtaunu_eta, "mc_truth_Wtaunu_eta/F", buffersize);
	tree->Branch("mc_truth_Wtaunutau_eta", &mc_truth_Wtaunutau_eta, "mc_truth_Wtaunutau_eta/F", buffersize);
	tree->Branch("mc_truth_Wtaul_eta", &mc_truth_Wtaul_eta, "mc_truth_Wtaul_eta/F", buffersize);
	tree->Branch("mc_truth_Wq1_eta", &mc_truth_Wq1_eta, "mc_truth_Wq1_eta/F", buffersize);
	tree->Branch("mc_truth_Wq2_eta", &mc_truth_Wq2_eta, "mc_truth_Wq2_eta/F", buffersize);

	tree->Branch("mc_truth_gammal1_eta", &mc_truth_gammal1_eta, "mc_truth_gammal1_eta/F", buffersize);
	tree->Branch("mc_truth_gammal2_eta", &mc_truth_gammal2_eta, "mc_truth_gammal2_eta/F", buffersize);
	tree->Branch("mc_truth_gammatau1_eta", &mc_truth_gammatau1_eta, "mc_truth_gammatau1_eta/F", buffersize);
	tree->Branch("mc_truth_gammatau2_eta", &mc_truth_gammatau2_eta, "mc_truth_gammatau2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataul1_eta", &mc_truth_gammataul1_eta, "mc_truth_gammataul1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataul2_eta", &mc_truth_gammataul2_eta, "mc_truth_gammataul2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_eta", &mc_truth_gammataunu1_eta, "mc_truth_gammataunu1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_eta", &mc_truth_gammataunu2_eta, "mc_truth_gammataunu2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_eta", &mc_truth_gammataunutau1_eta, "mc_truth_gammataunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_eta", &mc_truth_gammataunutau2_eta, "mc_truth_gammataunutau2_eta/F", buffersize);
	
	tree->Branch("mc_truth_t1_eta", &mc_truth_t1_eta, "mc_truth_t1_eta/F", buffersize);
	tree->Branch("mc_truth_t2_eta", &mc_truth_t2_eta, "mc_truth_t2_eta/F", buffersize);
	tree->Branch("mc_truth_tb1_eta", &mc_truth_tb1_eta, "mc_truth_tb1_eta/F", buffersize);
	tree->Branch("mc_truth_tb2_eta", &mc_truth_tb2_eta, "mc_truth_tb2_eta/F", buffersize);

	tree->Branch("mc_truth_tW1_eta", &mc_truth_tW1_eta, "mc_truth_tW1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu1_eta", &mc_truth_tWnu1_eta, "mc_truth_tWnu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_eta", &mc_truth_tWnutau1_eta, "mc_truth_tWnutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWl1_eta", &mc_truth_tWl1_eta, "mc_truth_tWl1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau1_eta", &mc_truth_tWtau1_eta, "mc_truth_tWtau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_eta", &mc_truth_tWtaunu1_eta, "mc_truth_tWtaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_eta", &mc_truth_tWtaunutau1_eta, "mc_truth_tWtaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_eta", &mc_truth_tWtaul1_eta, "mc_truth_tWtaul1_eta/F", buffersize);
	tree->Branch("mc_truth_tWq11_eta", &mc_truth_tWq11_eta, "mc_truth_tWq11_eta/F", buffersize);
	tree->Branch("mc_truth_tWq21_eta", &mc_truth_tWq21_eta, "mc_truth_tWq21_eta/F", buffersize);

	tree->Branch("mc_truth_tW2_eta", &mc_truth_tW2_eta, "mc_truth_tW2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu2_eta", &mc_truth_tWnu2_eta, "mc_truth_tWnu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_eta", &mc_truth_tWnutau2_eta, "mc_truth_tWnutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWl2_eta", &mc_truth_tWl2_eta, "mc_truth_tWl2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau2_eta", &mc_truth_tWtau2_eta, "mc_truth_tWtau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_eta", &mc_truth_tWtaunu2_eta, "mc_truth_tWtaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_eta", &mc_truth_tWtaunutau2_eta, "mc_truth_tWtaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_eta", &mc_truth_tWtaul2_eta, "mc_truth_tWtaul2_eta/F", buffersize);
	tree->Branch("mc_truth_tWq12_eta", &mc_truth_tWq12_eta, "mc_truth_tWq12_eta/F", buffersize);
	tree->Branch("mc_truth_tWq22_eta", &mc_truth_tWq22_eta, "mc_truth_tWq22_eta/F", buffersize);

	tree->Branch("mc_truth_j1_eta", &mc_truth_j1_eta, "mc_truth_j1_eta/F", buffersize);
	tree->Branch("mc_truth_j2_eta", &mc_truth_j2_eta, "mc_truth_j2_eta/F", buffersize);
	tree->Branch("mc_truth_j3_eta", &mc_truth_j3_eta, "mc_truth_j3_eta/F", buffersize);
		
	tree->Branch("mc_truth_W_phi", &mc_truth_W_phi, "mc_truth_W_phi/F", buffersize);
	tree->Branch("mc_truth_Wnu_phi", &mc_truth_Wnu_phi, "mc_truth_Wnu_phi/F", buffersize);
	tree->Branch("mc_truth_Wnutau_phi", &mc_truth_Wnutau_phi, "mc_truth_Wnutau_phi/F", buffersize);
	tree->Branch("mc_truth_Wl_phi", &mc_truth_Wl_phi, "mc_truth_Wl_phi/F", buffersize);
	tree->Branch("mc_truth_Wtau_phi", &mc_truth_Wtau_phi, "mc_truth_Wtau_phi/F", buffersize);
	tree->Branch("mc_truth_Wtaunu_phi", &mc_truth_Wtaunu_phi, "mc_truth_Wtaunu_phi/F", buffersize);
	tree->Branch("mc_truth_Wtaunutau_phi", &mc_truth_Wtaunutau_phi, "mc_truth_Wtaunutau_phi/F", buffersize);
	tree->Branch("mc_truth_Wtaul_phi", &mc_truth_Wtaul_phi, "mc_truth_Wtaul_phi/F", buffersize);
	tree->Branch("mc_truth_Wq1_phi", &mc_truth_Wq1_phi, "mc_truth_Wq1_phi/F", buffersize);
	tree->Branch("mc_truth_Wq2_phi", &mc_truth_Wq2_phi, "mc_truth_Wq2_phi/F", buffersize);

	tree->Branch("mc_truth_gammal1_phi", &mc_truth_gammal1_phi, "mc_truth_gammal1_phi/F", buffersize);
	tree->Branch("mc_truth_gammal2_phi", &mc_truth_gammal2_phi, "mc_truth_gammal2_phi/F", buffersize);
	tree->Branch("mc_truth_gammatau1_phi", &mc_truth_gammatau1_phi, "mc_truth_gammatau1_phi/F", buffersize);
	tree->Branch("mc_truth_gammatau2_phi", &mc_truth_gammatau2_phi, "mc_truth_gammatau2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataul1_phi", &mc_truth_gammataul1_phi, "mc_truth_gammataul1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataul2_phi", &mc_truth_gammataul2_phi, "mc_truth_gammataul2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_phi", &mc_truth_gammataunu1_phi, "mc_truth_gammataunu1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_phi", &mc_truth_gammataunu2_phi, "mc_truth_gammataunu2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_phi", &mc_truth_gammataunutau1_phi, "mc_truth_gammataunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_phi", &mc_truth_gammataunutau2_phi, "mc_truth_gammataunutau2_phi/F", buffersize);
	
	tree->Branch("mc_truth_t1_phi", &mc_truth_t1_phi, "mc_truth_t1_phi/F", buffersize);
	tree->Branch("mc_truth_t2_phi", &mc_truth_t2_phi, "mc_truth_t2_phi/F", buffersize);
	tree->Branch("mc_truth_tb1_phi", &mc_truth_tb1_phi, "mc_truth_tb1_phi/F", buffersize);
	tree->Branch("mc_truth_tb2_phi", &mc_truth_tb2_phi, "mc_truth_tb2_phi/F", buffersize);

	tree->Branch("mc_truth_tW1_phi", &mc_truth_tW1_phi, "mc_truth_tW1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu1_phi", &mc_truth_tWnu1_phi, "mc_truth_tWnu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_phi", &mc_truth_tWnutau1_phi, "mc_truth_tWnutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWl1_phi", &mc_truth_tWl1_phi, "mc_truth_tWl1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau1_phi", &mc_truth_tWtau1_phi, "mc_truth_tWtau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_phi", &mc_truth_tWtaunu1_phi, "mc_truth_tWtaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_phi", &mc_truth_tWtaunutau1_phi, "mc_truth_tWtaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_phi", &mc_truth_tWtaul1_phi, "mc_truth_tWtaul1_phi/F", buffersize);
	tree->Branch("mc_truth_tWq11_phi", &mc_truth_tWq11_phi, "mc_truth_tWq11_phi/F", buffersize);
	tree->Branch("mc_truth_tWq21_phi", &mc_truth_tWq21_phi, "mc_truth_tWq21_phi/F", buffersize);

	tree->Branch("mc_truth_tW2_phi", &mc_truth_tW2_phi, "mc_truth_tW2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu2_phi", &mc_truth_tWnu2_phi, "mc_truth_tWnu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_phi", &mc_truth_tWnutau2_phi, "mc_truth_tWnutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWl2_phi", &mc_truth_tWl2_phi, "mc_truth_tWl2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau2_phi", &mc_truth_tWtau2_phi, "mc_truth_tWtau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_phi", &mc_truth_tWtaunu2_phi, "mc_truth_tWtaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_phi", &mc_truth_tWtaunutau2_phi, "mc_truth_tWtaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_phi", &mc_truth_tWtaul2_phi, "mc_truth_tWtaul2_phi/F", buffersize);
	tree->Branch("mc_truth_tWq12_phi", &mc_truth_tWq12_phi, "mc_truth_tWq12_phi/F", buffersize);
	tree->Branch("mc_truth_tWq22_phi", &mc_truth_tWq22_phi, "mc_truth_tWq22_phi/F", buffersize);

	tree->Branch("mc_truth_j1_phi", &mc_truth_j1_phi, "mc_truth_j1_phi/F", buffersize);
	tree->Branch("mc_truth_j2_phi", &mc_truth_j2_phi, "mc_truth_j2_phi/F", buffersize);
	tree->Branch("mc_truth_j3_phi", &mc_truth_j3_phi, "mc_truth_j3_phi/F", buffersize);

	tree->Branch("mc_truth_W_E", &mc_truth_W_E, "mc_truth_W_E/F", buffersize);
	tree->Branch("mc_truth_Wnu_E", &mc_truth_Wnu_E, "mc_truth_Wnu_E/F", buffersize);
	tree->Branch("mc_truth_Wnutau_E", &mc_truth_Wnutau_E, "mc_truth_Wnutau_E/F", buffersize);
	tree->Branch("mc_truth_Wl_E", &mc_truth_Wl_E, "mc_truth_Wl_E/F", buffersize);
	tree->Branch("mc_truth_Wtau_E", &mc_truth_Wtau_E, "mc_truth_Wtau_E/F", buffersize);
	tree->Branch("mc_truth_Wtaunu_E", &mc_truth_Wtaunu_E, "mc_truth_Wtaunu_E/F", buffersize);
	tree->Branch("mc_truth_Wtaunutau_E", &mc_truth_Wtaunutau_E, "mc_truth_Wtaunutau_E/F", buffersize);
	tree->Branch("mc_truth_Wtaul_E", &mc_truth_Wtaul_E, "mc_truth_Wtaul_E/F", buffersize);
	tree->Branch("mc_truth_Wq1_E", &mc_truth_Wq1_E, "mc_truth_Wq1_E/F", buffersize);
	tree->Branch("mc_truth_Wq2_E", &mc_truth_Wq2_E, "mc_truth_Wq2_E/F", buffersize);

	tree->Branch("mc_truth_gammal1_E", &mc_truth_gammal1_E, "mc_truth_gammal1_E/F", buffersize);
	tree->Branch("mc_truth_gammal2_E", &mc_truth_gammal2_E, "mc_truth_gammal2_E/F", buffersize);
	tree->Branch("mc_truth_gammatau1_E", &mc_truth_gammatau1_E, "mc_truth_gammatau1_E/F", buffersize);
	tree->Branch("mc_truth_gammatau2_E", &mc_truth_gammatau2_E, "mc_truth_gammatau2_E/F", buffersize);
	tree->Branch("mc_truth_gammataul1_E", &mc_truth_gammataul1_E, "mc_truth_gammataul1_E/F", buffersize);
	tree->Branch("mc_truth_gammataul2_E", &mc_truth_gammataul2_E, "mc_truth_gammataul2_E/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_E", &mc_truth_gammataunu1_E, "mc_truth_gammataunu1_E/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_E", &mc_truth_gammataunu2_E, "mc_truth_gammataunu2_E/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_E", &mc_truth_gammataunutau1_E, "mc_truth_gammataunutau1_E/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_E", &mc_truth_gammataunutau2_E, "mc_truth_gammataunutau2_E/F", buffersize);
	
	tree->Branch("mc_truth_t1_E", &mc_truth_t1_E, "mc_truth_t1_E/F", buffersize);
	tree->Branch("mc_truth_t2_E", &mc_truth_t2_E, "mc_truth_t2_E/F", buffersize);
	tree->Branch("mc_truth_tb1_E", &mc_truth_tb1_E, "mc_truth_tb1_E/F", buffersize);
	tree->Branch("mc_truth_tb2_E", &mc_truth_tb2_E, "mc_truth_tb2_E/F", buffersize);

	tree->Branch("mc_truth_tW1_E", &mc_truth_tW1_E, "mc_truth_tW1_E/F", buffersize);
	tree->Branch("mc_truth_tWnu1_E", &mc_truth_tWnu1_E, "mc_truth_tWnu1_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_E", &mc_truth_tWnutau1_E, "mc_truth_tWnutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWl1_E", &mc_truth_tWl1_E, "mc_truth_tWl1_E/F", buffersize);
	tree->Branch("mc_truth_tWtau1_E", &mc_truth_tWtau1_E, "mc_truth_tWtau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_E", &mc_truth_tWtaunu1_E, "mc_truth_tWtaunu1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_E", &mc_truth_tWtaunutau1_E, "mc_truth_tWtaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_E", &mc_truth_tWtaul1_E, "mc_truth_tWtaul1_E/F", buffersize);
	tree->Branch("mc_truth_tWq11_E", &mc_truth_tWq11_E, "mc_truth_tWq11_E/F", buffersize);
	tree->Branch("mc_truth_tWq21_E", &mc_truth_tWq21_E, "mc_truth_tWq21_E/F", buffersize);

	tree->Branch("mc_truth_tW2_E", &mc_truth_tW2_E, "mc_truth_tW2_E/F", buffersize);
	tree->Branch("mc_truth_tWnu2_E", &mc_truth_tWnu2_E, "mc_truth_tWnu2_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_E", &mc_truth_tWnutau2_E, "mc_truth_tWnutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWl2_E", &mc_truth_tWl2_E, "mc_truth_tWl2_E/F", buffersize);
	tree->Branch("mc_truth_tWtau2_E", &mc_truth_tWtau2_E, "mc_truth_tWtau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_E", &mc_truth_tWtaunu2_E, "mc_truth_tWtaunu2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_E", &mc_truth_tWtaunutau2_E, "mc_truth_tWtaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_E", &mc_truth_tWtaul2_E, "mc_truth_tWtaul2_E/F", buffersize);
	tree->Branch("mc_truth_tWq12_E", &mc_truth_tWq12_E, "mc_truth_tWq12_E/F", buffersize);
	tree->Branch("mc_truth_tWq22_E", &mc_truth_tWq22_E, "mc_truth_tWq22_E/F", buffersize);

	tree->Branch("mc_truth_j1_E", &mc_truth_j1_E, "mc_truth_j1_E/F", buffersize);
	tree->Branch("mc_truth_j2_E", &mc_truth_j2_E, "mc_truth_j2_E/F", buffersize);
	tree->Branch("mc_truth_j3_E", &mc_truth_j3_E, "mc_truth_j3_E/F", buffersize);
			
	tree->Branch("mc_truth_W_id", &mc_truth_W_id, "mc_truth_W_id/I", buffersize);
	tree->Branch("mc_truth_Wnu_id", &mc_truth_Wnu_id, "mc_truth_Wnu_id/I", buffersize);
	tree->Branch("mc_truth_Wnutau_id", &mc_truth_Wnutau_id, "mc_truth_Wnutau_id/I", buffersize);
	tree->Branch("mc_truth_Wl_id", &mc_truth_Wl_id, "mc_truth_Wl_id/I", buffersize);
	tree->Branch("mc_truth_Wtau_id", &mc_truth_Wtau_id, "mc_truth_Wtau_id/I", buffersize);
	tree->Branch("mc_truth_Wtaunu_id", &mc_truth_Wtaunu_id, "mc_truth_Wtaunu_id/I", buffersize);
	tree->Branch("mc_truth_Wtaunutau_id", &mc_truth_Wtaunutau_id, "mc_truth_Wtaunutau_id/I", buffersize);
	tree->Branch("mc_truth_Wtaul_id", &mc_truth_Wtaul_id, "mc_truth_Wtaul_id/I", buffersize);
	tree->Branch("mc_truth_Wq1_id", &mc_truth_Wq1_id, "mc_truth_Wq1_id/I", buffersize);
	tree->Branch("mc_truth_Wq2_id", &mc_truth_Wq2_id, "mc_truth_Wq2_id/I", buffersize);

	tree->Branch("mc_truth_gammal1_id", &mc_truth_gammal1_id, "mc_truth_gammal1_id/I", buffersize);
	tree->Branch("mc_truth_gammal2_id", &mc_truth_gammal2_id, "mc_truth_gammal2_id/I", buffersize);
	tree->Branch("mc_truth_gammatau1_id", &mc_truth_gammatau1_id, "mc_truth_gammatau1_id/I", buffersize);
	tree->Branch("mc_truth_gammatau2_id", &mc_truth_gammatau2_id, "mc_truth_gammatau2_id/I", buffersize);
	tree->Branch("mc_truth_gammataul1_id", &mc_truth_gammataul1_id, "mc_truth_gammataul1_id/I", buffersize);
	tree->Branch("mc_truth_gammataul2_id", &mc_truth_gammataul2_id, "mc_truth_gammataul2_id/I", buffersize);
	tree->Branch("mc_truth_gammataunu1_id", &mc_truth_gammataunu1_id, "mc_truth_gammataunu1_id/I", buffersize);
	tree->Branch("mc_truth_gammataunu2_id", &mc_truth_gammataunu2_id, "mc_truth_gammataunu2_id/I", buffersize);
	tree->Branch("mc_truth_gammataunutau1_id", &mc_truth_gammataunutau1_id, "mc_truth_gammataunutau1_id/I", buffersize);
	tree->Branch("mc_truth_gammataunutau2_id", &mc_truth_gammataunutau2_id, "mc_truth_gammataunutau2_id/I", buffersize);
	
	tree->Branch("mc_truth_t1_id", &mc_truth_t1_id, "mc_truth_t1_id/I", buffersize);
	tree->Branch("mc_truth_t2_id", &mc_truth_t2_id, "mc_truth_t2_id/I", buffersize);
	tree->Branch("mc_truth_tb1_id", &mc_truth_tb1_id, "mc_truth_tb1_id/I", buffersize);
	tree->Branch("mc_truth_tb2_id", &mc_truth_tb2_id, "mc_truth_tb2_id/I", buffersize);

	tree->Branch("mc_truth_tW1_id", &mc_truth_tW1_id, "mc_truth_tW1_id/I", buffersize);
	tree->Branch("mc_truth_tWnu1_id", &mc_truth_tWnu1_id, "mc_truth_tWnu1_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_id", &mc_truth_tWnutau1_id, "mc_truth_tWnutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWl1_id", &mc_truth_tWl1_id, "mc_truth_tWl1_id/I", buffersize);
	tree->Branch("mc_truth_tWtau1_id", &mc_truth_tWtau1_id, "mc_truth_tWtau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_id", &mc_truth_tWtaunu1_id, "mc_truth_tWtaunu1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_id", &mc_truth_tWtaunutau1_id, "mc_truth_tWtaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_id", &mc_truth_tWtaul1_id, "mc_truth_tWtaul1_id/I", buffersize);
	tree->Branch("mc_truth_tWq11_id", &mc_truth_tWq11_id, "mc_truth_tWq11_id/I", buffersize);
	tree->Branch("mc_truth_tWq21_id", &mc_truth_tWq21_id, "mc_truth_tWq21_id/I", buffersize);

	tree->Branch("mc_truth_tW2_id", &mc_truth_tW2_id, "mc_truth_tW2_id/I", buffersize);
	tree->Branch("mc_truth_tWnu2_id", &mc_truth_tWnu2_id, "mc_truth_tWnu2_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_id", &mc_truth_tWnutau2_id, "mc_truth_tWnutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWl2_id", &mc_truth_tWl2_id, "mc_truth_tWl2_id/I", buffersize);
	tree->Branch("mc_truth_tWtau2_id", &mc_truth_tWtau2_id, "mc_truth_tWtau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_id", &mc_truth_tWtaunu2_id, "mc_truth_tWtaunu2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_id", &mc_truth_tWtaunutau2_id, "mc_truth_tWtaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_id", &mc_truth_tWtaul2_id, "mc_truth_tWtaul2_id/I", buffersize);
	tree->Branch("mc_truth_tWq12_id", &mc_truth_tWq12_id, "mc_truth_tWq12_id/I", buffersize);
	tree->Branch("mc_truth_tWq22_id", &mc_truth_tWq22_id, "mc_truth_tWq22_id/I", buffersize);

	tree->Branch("mc_truth_j1_id", &mc_truth_j1_id, "mc_truth_j1_id/I", buffersize);
	tree->Branch("mc_truth_j2_id", &mc_truth_j2_id, "mc_truth_j2_id/I", buffersize);
	tree->Branch("mc_truth_j3_id", &mc_truth_j3_id, "mc_truth_j3_id/I", buffersize);

	tree->Branch("mc_truth_W_status", &mc_truth_W_status, "mc_truth_W_status/I", buffersize);
	tree->Branch("mc_truth_Wnu_status", &mc_truth_Wnu_status, "mc_truth_Wnu_status/I", buffersize);
	tree->Branch("mc_truth_Wnutau_status", &mc_truth_Wnutau_status, "mc_truth_Wnutau_status/I", buffersize);
	tree->Branch("mc_truth_Wl_status", &mc_truth_Wl_status, "mc_truth_Wl_status/I", buffersize);
	tree->Branch("mc_truth_Wtau_status", &mc_truth_Wtau_status, "mc_truth_Wtau_status/I", buffersize);
	tree->Branch("mc_truth_Wtaunu_status", &mc_truth_Wtaunu_status, "mc_truth_Wtaunu_status/I", buffersize);
	tree->Branch("mc_truth_Wtaunutau_status", &mc_truth_Wtaunutau_status, "mc_truth_Wtaunutau_status/I", buffersize);
	tree->Branch("mc_truth_Wtaul_status", &mc_truth_Wtaul_status, "mc_truth_Wtaul_status/I", buffersize);
	tree->Branch("mc_truth_Wq1_status", &mc_truth_Wq1_status, "mc_truth_Wq1_status/I", buffersize);
	tree->Branch("mc_truth_Wq2_status", &mc_truth_Wq2_status, "mc_truth_Wq2_status/I", buffersize);

	tree->Branch("mc_truth_gammal1_status", &mc_truth_gammal1_status, "mc_truth_gammal1_status/I", buffersize);
	tree->Branch("mc_truth_gammal2_status", &mc_truth_gammal2_status, "mc_truth_gammal2_status/I", buffersize);
	tree->Branch("mc_truth_gammatau1_status", &mc_truth_gammatau1_status, "mc_truth_gammatau1_status/I", buffersize);
	tree->Branch("mc_truth_gammatau2_status", &mc_truth_gammatau2_status, "mc_truth_gammatau2_status/I", buffersize);
	tree->Branch("mc_truth_gammataul1_status", &mc_truth_gammataul1_status, "mc_truth_gammataul1_status/I", buffersize);
	tree->Branch("mc_truth_gammataul2_status", &mc_truth_gammataul2_status, "mc_truth_gammataul2_status/I", buffersize);
	tree->Branch("mc_truth_gammataunu1_status", &mc_truth_gammataunu1_status, "mc_truth_gammataunu1_status/I", buffersize);
	tree->Branch("mc_truth_gammataunu2_status", &mc_truth_gammataunu2_status, "mc_truth_gammataunu2_status/I", buffersize);
	tree->Branch("mc_truth_gammataunutau1_status", &mc_truth_gammataunutau1_status, "mc_truth_gammataunutau1_status/I", buffersize);
	tree->Branch("mc_truth_gammataunutau2_status", &mc_truth_gammataunutau2_status, "mc_truth_gammataunutau2_status/I", buffersize);
	
	tree->Branch("mc_truth_t1_status", &mc_truth_t1_status, "mc_truth_t1_status/I", buffersize);
	tree->Branch("mc_truth_t2_status", &mc_truth_t2_status, "mc_truth_t2_status/I", buffersize);
	tree->Branch("mc_truth_tb1_status", &mc_truth_tb1_status, "mc_truth_tb1_status/I", buffersize);
	tree->Branch("mc_truth_tb2_status", &mc_truth_tb2_status, "mc_truth_tb2_status/I", buffersize);

	tree->Branch("mc_truth_tW1_status", &mc_truth_tW1_status, "mc_truth_tW1_status/I", buffersize);
	tree->Branch("mc_truth_tWnu1_status", &mc_truth_tWnu1_status, "mc_truth_tWnu1_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_status", &mc_truth_tWnutau1_status, "mc_truth_tWnutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWl1_status", &mc_truth_tWl1_status, "mc_truth_tWl1_status/I", buffersize);
	tree->Branch("mc_truth_tWtau1_status", &mc_truth_tWtau1_status, "mc_truth_tWtau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_status", &mc_truth_tWtaunu1_status, "mc_truth_tWtaunu1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_status", &mc_truth_tWtaunutau1_status, "mc_truth_tWtaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_status", &mc_truth_tWtaul1_status, "mc_truth_tWtaul1_status/I", buffersize);
	tree->Branch("mc_truth_tWq11_status", &mc_truth_tWq11_status, "mc_truth_tWq11_status/I", buffersize);
	tree->Branch("mc_truth_tWq21_status", &mc_truth_tWq21_status, "mc_truth_tWq21_status/I", buffersize);

	tree->Branch("mc_truth_tW2_status", &mc_truth_tW2_status, "mc_truth_tW2_status/I", buffersize);
	tree->Branch("mc_truth_tWnu2_status", &mc_truth_tWnu2_status, "mc_truth_tWnu2_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_status", &mc_truth_tWnutau2_status, "mc_truth_tWnutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWl2_status", &mc_truth_tWl2_status, "mc_truth_tWl2_status/I", buffersize);
	tree->Branch("mc_truth_tWtau2_status", &mc_truth_tWtau2_status, "mc_truth_tWtau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_status", &mc_truth_tWtaunu2_status, "mc_truth_tWtaunu2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_status", &mc_truth_tWtaunutau2_status, "mc_truth_tWtaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_status", &mc_truth_tWtaul2_status, "mc_truth_tWtaul2_status/I", buffersize);
	tree->Branch("mc_truth_tWq12_status", &mc_truth_tWq12_status, "mc_truth_tWq12_status/I", buffersize);
	tree->Branch("mc_truth_tWq22_status", &mc_truth_tWq22_status, "mc_truth_tWq22_status/I", buffersize);

	tree->Branch("mc_truth_j1_status", &mc_truth_j1_status, "mc_truth_j1_status/I", buffersize);
	tree->Branch("mc_truth_j2_status", &mc_truth_j2_status, "mc_truth_j2_status/I", buffersize);
	tree->Branch("mc_truth_j3_status", &mc_truth_j3_status, "mc_truth_j3_status/I", buffersize);
     }
   
   if( doWrite("mc_truth_tzq") )
     {
	tree->Branch("mc_truth_tzq_channel", &mc_truth_tzq_channel, "mc_truth_tzq_channel/I", buffersize);

	if( doWrite("mc_truth_p4") )
	  {	
	     tree->Branch("mc_truth_Z_p4", "TLorentzVector", &mc_truth_Z_p4, buffersize);
	     tree->Branch("mc_truth_Zl1_p4", "TLorentzVector", &mc_truth_Zl1_p4, buffersize);
	     tree->Branch("mc_truth_Zl2_p4", "TLorentzVector", &mc_truth_Zl2_p4, buffersize);
	     tree->Branch("mc_truth_Ztau1_p4", "TLorentzVector", &mc_truth_Ztau1_p4, buffersize);
	     tree->Branch("mc_truth_Ztau2_p4", "TLorentzVector", &mc_truth_Ztau2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaul1_p4", "TLorentzVector", &mc_truth_Ztaul1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaul2_p4", "TLorentzVector", &mc_truth_Ztaul2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunu1_p4", "TLorentzVector", &mc_truth_Ztaunu1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunu2_p4", "TLorentzVector", &mc_truth_Ztaunu2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunutau1_p4", "TLorentzVector", &mc_truth_Ztaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunutau2_p4", "TLorentzVector", &mc_truth_Ztaunutau2_p4, buffersize);
	     
	     tree->Branch("mc_truth_t_p4", "TLorentzVector", &mc_truth_t_p4, buffersize);
	     tree->Branch("mc_truth_tb_p4", "TLorentzVector", &mc_truth_tb_p4, buffersize);
	     tree->Branch("mc_truth_tW_p4", "TLorentzVector", &mc_truth_tW_p4, buffersize);
	     tree->Branch("mc_truth_tWnu_p4", "TLorentzVector", &mc_truth_tWnu_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau_p4", "TLorentzVector", &mc_truth_tWnutau_p4, buffersize);
	     tree->Branch("mc_truth_tWl_p4", "TLorentzVector", &mc_truth_tWl_p4, buffersize);
	     tree->Branch("mc_truth_tWtau_p4", "TLorentzVector", &mc_truth_tWtau_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu_p4", "TLorentzVector", &mc_truth_tWtaunu_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau_p4", "TLorentzVector", &mc_truth_tWtaunutau_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul_p4", "TLorentzVector", &mc_truth_tWtaul_p4, buffersize);
	     tree->Branch("mc_truth_tWq1_p4", "TLorentzVector", &mc_truth_tWq1_p4, buffersize);
	     tree->Branch("mc_truth_tWq2_p4", "TLorentzVector", &mc_truth_tWq2_p4, buffersize);
	     
	     tree->Branch("mc_truth_j1_p4", "TLorentzVector", &mc_truth_j1_p4, buffersize);
	     tree->Branch("mc_truth_j2_p4", "TLorentzVector", &mc_truth_j2_p4, buffersize);
	     tree->Branch("mc_truth_j3_p4", "TLorentzVector", &mc_truth_j3_p4, buffersize);
	  }

	tree->Branch("mc_truth_Z_pt", &mc_truth_Z_pt, "mc_truth_Z_pt/F", buffersize);
	tree->Branch("mc_truth_Zl1_pt", &mc_truth_Zl1_pt, "mc_truth_Zl1_pt/F", buffersize);
	tree->Branch("mc_truth_Zl2_pt", &mc_truth_Zl2_pt, "mc_truth_Zl2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztau1_pt", &mc_truth_Ztau1_pt, "mc_truth_Ztau1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztau2_pt", &mc_truth_Ztau2_pt, "mc_truth_Ztau2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_pt", &mc_truth_Ztaul1_pt, "mc_truth_Ztaul1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_pt", &mc_truth_Ztaul2_pt, "mc_truth_Ztaul2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_pt", &mc_truth_Ztaunu1_pt, "mc_truth_Ztaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_pt", &mc_truth_Ztaunu2_pt, "mc_truth_Ztaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_pt", &mc_truth_Ztaunutau1_pt, "mc_truth_Ztaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_pt", &mc_truth_Ztaunutau2_pt, "mc_truth_Ztaunutau2_pt/F", buffersize);

	tree->Branch("mc_truth_t_pt", &mc_truth_t_pt, "mc_truth_t_pt/F", buffersize);
	tree->Branch("mc_truth_tb_pt", &mc_truth_tb_pt, "mc_truth_tb_pt/F", buffersize);
	tree->Branch("mc_truth_tW_pt", &mc_truth_tW_pt, "mc_truth_tW_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu_pt", &mc_truth_tWnu_pt, "mc_truth_tWnu_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau_pt", &mc_truth_tWnutau_pt, "mc_truth_tWnutau_pt/F", buffersize);
	tree->Branch("mc_truth_tWl_pt", &mc_truth_tWl_pt, "mc_truth_tWl_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau_pt", &mc_truth_tWtau_pt, "mc_truth_tWtau_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_pt", &mc_truth_tWtaunu_pt, "mc_truth_tWtaunu_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_pt", &mc_truth_tWtaunutau_pt, "mc_truth_tWtaunutau_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul_pt", &mc_truth_tWtaul_pt, "mc_truth_tWtaul_pt/F", buffersize);
	tree->Branch("mc_truth_tWq1_pt", &mc_truth_tWq1_pt, "mc_truth_tWq1_pt/F", buffersize);
	tree->Branch("mc_truth_tWq2_pt", &mc_truth_tWq2_pt, "mc_truth_tWq2_pt/F", buffersize);

	tree->Branch("mc_truth_j1_pt", &mc_truth_j1_pt, "mc_truth_j1_pt/F", buffersize);
	tree->Branch("mc_truth_j2_pt", &mc_truth_j2_pt, "mc_truth_j2_pt/F", buffersize);
	tree->Branch("mc_truth_j3_pt", &mc_truth_j3_pt, "mc_truth_j3_pt/F", buffersize);	
	
	tree->Branch("mc_truth_Z_eta", &mc_truth_Z_eta, "mc_truth_Z_eta/F", buffersize);
	tree->Branch("mc_truth_Zl1_eta", &mc_truth_Zl1_eta, "mc_truth_Zl1_eta/F", buffersize);
	tree->Branch("mc_truth_Zl2_eta", &mc_truth_Zl2_eta, "mc_truth_Zl2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztau1_eta", &mc_truth_Ztau1_eta, "mc_truth_Ztau1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztau2_eta", &mc_truth_Ztau2_eta, "mc_truth_Ztau2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_eta", &mc_truth_Ztaul1_eta, "mc_truth_Ztaul1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_eta", &mc_truth_Ztaul2_eta, "mc_truth_Ztaul2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_eta", &mc_truth_Ztaunu1_eta, "mc_truth_Ztaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_eta", &mc_truth_Ztaunu2_eta, "mc_truth_Ztaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_eta", &mc_truth_Ztaunutau1_eta, "mc_truth_Ztaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_eta", &mc_truth_Ztaunutau2_eta, "mc_truth_Ztaunutau2_eta/F", buffersize);

	tree->Branch("mc_truth_t_eta", &mc_truth_t_eta, "mc_truth_t_eta/F", buffersize);
	tree->Branch("mc_truth_tb_eta", &mc_truth_tb_eta, "mc_truth_tb_eta/F", buffersize);
	tree->Branch("mc_truth_tW_eta", &mc_truth_tW_eta, "mc_truth_tW_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu_eta", &mc_truth_tWnu_eta, "mc_truth_tWnu_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau_eta", &mc_truth_tWnutau_eta, "mc_truth_tWnutau_eta/F", buffersize);
	tree->Branch("mc_truth_tWl_eta", &mc_truth_tWl_eta, "mc_truth_tWl_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau_eta", &mc_truth_tWtau_eta, "mc_truth_tWtau_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_eta", &mc_truth_tWtaunu_eta, "mc_truth_tWtaunu_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_eta", &mc_truth_tWtaunutau_eta, "mc_truth_tWtaunutau_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul_eta", &mc_truth_tWtaul_eta, "mc_truth_tWtaul_eta/F", buffersize);
	tree->Branch("mc_truth_tWq1_eta", &mc_truth_tWq1_eta, "mc_truth_tWq1_eta/F", buffersize);
	tree->Branch("mc_truth_tWq2_eta", &mc_truth_tWq2_eta, "mc_truth_tWq2_eta/F", buffersize);

	tree->Branch("mc_truth_j1_eta", &mc_truth_j1_eta, "mc_truth_j1_eta/F", buffersize);
	tree->Branch("mc_truth_j2_eta", &mc_truth_j2_eta, "mc_truth_j2_eta/F", buffersize);
	tree->Branch("mc_truth_j3_eta", &mc_truth_j3_eta, "mc_truth_j3_eta/F", buffersize);
	
	tree->Branch("mc_truth_Z_phi", &mc_truth_Z_phi, "mc_truth_Z_phi/F", buffersize);
	tree->Branch("mc_truth_Zl1_phi", &mc_truth_Zl1_phi, "mc_truth_Zl1_phi/F", buffersize);
	tree->Branch("mc_truth_Zl2_phi", &mc_truth_Zl2_phi, "mc_truth_Zl2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztau1_phi", &mc_truth_Ztau1_phi, "mc_truth_Ztau1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztau2_phi", &mc_truth_Ztau2_phi, "mc_truth_Ztau2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_phi", &mc_truth_Ztaul1_phi, "mc_truth_Ztaul1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_phi", &mc_truth_Ztaul2_phi, "mc_truth_Ztaul2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_phi", &mc_truth_Ztaunu1_phi, "mc_truth_Ztaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_phi", &mc_truth_Ztaunu2_phi, "mc_truth_Ztaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_phi", &mc_truth_Ztaunutau1_phi, "mc_truth_Ztaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_phi", &mc_truth_Ztaunutau2_phi, "mc_truth_Ztaunutau2_phi/F", buffersize);

	tree->Branch("mc_truth_t_phi", &mc_truth_t_phi, "mc_truth_t_phi/F", buffersize);
	tree->Branch("mc_truth_tb_phi", &mc_truth_tb_phi, "mc_truth_tb_phi/F", buffersize);
	tree->Branch("mc_truth_tW_phi", &mc_truth_tW_phi, "mc_truth_tW_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu_phi", &mc_truth_tWnu_phi, "mc_truth_tWnu_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau_phi", &mc_truth_tWnutau_phi, "mc_truth_tWnutau_phi/F", buffersize);
	tree->Branch("mc_truth_tWl_phi", &mc_truth_tWl_phi, "mc_truth_tWl_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau_phi", &mc_truth_tWtau_phi, "mc_truth_tWtau_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_phi", &mc_truth_tWtaunu_phi, "mc_truth_tWtaunu_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_phi", &mc_truth_tWtaunutau_phi, "mc_truth_tWtaunutau_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul_phi", &mc_truth_tWtaul_phi, "mc_truth_tWtaul_phi/F", buffersize);
	tree->Branch("mc_truth_tWq1_phi", &mc_truth_tWq1_phi, "mc_truth_tWq1_phi/F", buffersize);
	tree->Branch("mc_truth_tWq2_phi", &mc_truth_tWq2_phi, "mc_truth_tWq2_phi/F", buffersize);

	tree->Branch("mc_truth_j1_phi", &mc_truth_j1_phi, "mc_truth_j1_phi/F", buffersize);
	tree->Branch("mc_truth_j2_phi", &mc_truth_j2_phi, "mc_truth_j2_phi/F", buffersize);
	tree->Branch("mc_truth_j3_phi", &mc_truth_j3_phi, "mc_truth_j3_phi/F", buffersize);
	
	tree->Branch("mc_truth_Z_E", &mc_truth_Z_E, "mc_truth_Z_E/F", buffersize);
	tree->Branch("mc_truth_Zl1_E", &mc_truth_Zl1_E, "mc_truth_Zl1_E/F", buffersize);
	tree->Branch("mc_truth_Zl2_E", &mc_truth_Zl2_E, "mc_truth_Zl2_E/F", buffersize);
	tree->Branch("mc_truth_Ztau1_E", &mc_truth_Ztau1_E, "mc_truth_Ztau1_E/F", buffersize);
	tree->Branch("mc_truth_Ztau2_E", &mc_truth_Ztau2_E, "mc_truth_Ztau2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_E", &mc_truth_Ztaul1_E, "mc_truth_Ztaul1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_E", &mc_truth_Ztaul2_E, "mc_truth_Ztaul2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_E", &mc_truth_Ztaunu1_E, "mc_truth_Ztaunu1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_E", &mc_truth_Ztaunu2_E, "mc_truth_Ztaunu2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_E", &mc_truth_Ztaunutau1_E, "mc_truth_Ztaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_E", &mc_truth_Ztaunutau2_E, "mc_truth_Ztaunutau2_E/F", buffersize);

	tree->Branch("mc_truth_t_E", &mc_truth_t_E, "mc_truth_t_E/F", buffersize);
	tree->Branch("mc_truth_tb_E", &mc_truth_tb_E, "mc_truth_tb_E/F", buffersize);
	tree->Branch("mc_truth_tW_E", &mc_truth_tW_E, "mc_truth_tW_E/F", buffersize);
	tree->Branch("mc_truth_tWnu_E", &mc_truth_tWnu_E, "mc_truth_tWnu_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau_E", &mc_truth_tWnutau_E, "mc_truth_tWnutau_E/F", buffersize);
	tree->Branch("mc_truth_tWl_E", &mc_truth_tWl_E, "mc_truth_tWl_E/F", buffersize);
	tree->Branch("mc_truth_tWtau_E", &mc_truth_tWtau_E, "mc_truth_tWtau_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_E", &mc_truth_tWtaunu_E, "mc_truth_tWtaunu_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_E", &mc_truth_tWtaunutau_E, "mc_truth_tWtaunutau_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul_E", &mc_truth_tWtaul_E, "mc_truth_tWtaul_E/F", buffersize);
	tree->Branch("mc_truth_tWq1_E", &mc_truth_tWq1_E, "mc_truth_tWq1_E/F", buffersize);
	tree->Branch("mc_truth_tWq2_E", &mc_truth_tWq2_E, "mc_truth_tWq2_E/F", buffersize);

	tree->Branch("mc_truth_j1_E", &mc_truth_j1_E, "mc_truth_j1_E/F", buffersize);
	tree->Branch("mc_truth_j2_E", &mc_truth_j2_E, "mc_truth_j2_E/F", buffersize);
	tree->Branch("mc_truth_j3_E", &mc_truth_j3_E, "mc_truth_j3_E/F", buffersize);
		
	tree->Branch("mc_truth_Z_id", &mc_truth_Z_id, "mc_truth_Z_id/I", buffersize);
	tree->Branch("mc_truth_Zl1_id", &mc_truth_Zl1_id, "mc_truth_Zl1_id/I", buffersize);
	tree->Branch("mc_truth_Zl2_id", &mc_truth_Zl2_id, "mc_truth_Zl2_id/I", buffersize);
	tree->Branch("mc_truth_Ztau1_id", &mc_truth_Ztau1_id, "mc_truth_Ztau1_id/I", buffersize);
	tree->Branch("mc_truth_Ztau2_id", &mc_truth_Ztau2_id, "mc_truth_Ztau2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaul1_id", &mc_truth_Ztaul1_id, "mc_truth_Ztaul1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaul2_id", &mc_truth_Ztaul2_id, "mc_truth_Ztaul2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunu1_id", &mc_truth_Ztaunu1_id, "mc_truth_Ztaunu1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunu2_id", &mc_truth_Ztaunu2_id, "mc_truth_Ztaunu2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_id", &mc_truth_Ztaunutau1_id, "mc_truth_Ztaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_id", &mc_truth_Ztaunutau2_id, "mc_truth_Ztaunutau2_id/I", buffersize);

	tree->Branch("mc_truth_t_id", &mc_truth_t_id, "mc_truth_t_id/I", buffersize);
	tree->Branch("mc_truth_tb_id", &mc_truth_tb_id, "mc_truth_tb_id/I", buffersize);
	tree->Branch("mc_truth_tW_id", &mc_truth_tW_id, "mc_truth_tW_id/I", buffersize);
	tree->Branch("mc_truth_tWnu_id", &mc_truth_tWnu_id, "mc_truth_tWnu_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau_id", &mc_truth_tWnutau_id, "mc_truth_tWnutau_id/I", buffersize);
	tree->Branch("mc_truth_tWl_id", &mc_truth_tWl_id, "mc_truth_tWl_id/I", buffersize);
	tree->Branch("mc_truth_tWtau_id", &mc_truth_tWtau_id, "mc_truth_tWtau_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu_id", &mc_truth_tWtaunu_id, "mc_truth_tWtaunu_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau_id", &mc_truth_tWtaunutau_id, "mc_truth_tWtaunutau_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul_id", &mc_truth_tWtaul_id, "mc_truth_tWtaul_id/I", buffersize);
	tree->Branch("mc_truth_tWq1_id", &mc_truth_tWq1_id, "mc_truth_tWq1_id/I", buffersize);
	tree->Branch("mc_truth_tWq2_id", &mc_truth_tWq2_id, "mc_truth_tWq2_id/I", buffersize);

	tree->Branch("mc_truth_j1_id", &mc_truth_j1_id, "mc_truth_j1_id/I", buffersize);
	tree->Branch("mc_truth_j2_id", &mc_truth_j2_id, "mc_truth_j2_id/I", buffersize);
	tree->Branch("mc_truth_j3_id", &mc_truth_j3_id, "mc_truth_j3_id/I", buffersize);

	tree->Branch("mc_truth_Z_status", &mc_truth_Z_status, "mc_truth_Z_status/I", buffersize);
	tree->Branch("mc_truth_Zl1_status", &mc_truth_Zl1_status, "mc_truth_Zl1_status/I", buffersize);
	tree->Branch("mc_truth_Zl2_status", &mc_truth_Zl2_status, "mc_truth_Zl2_status/I", buffersize);
	tree->Branch("mc_truth_Ztau1_status", &mc_truth_Ztau1_status, "mc_truth_Ztau1_status/I", buffersize);
	tree->Branch("mc_truth_Ztau2_status", &mc_truth_Ztau2_status, "mc_truth_Ztau2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaul1_status", &mc_truth_Ztaul1_status, "mc_truth_Ztaul1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaul2_status", &mc_truth_Ztaul2_status, "mc_truth_Ztaul2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunu1_status", &mc_truth_Ztaunu1_status, "mc_truth_Ztaunu1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunu2_status", &mc_truth_Ztaunu2_status, "mc_truth_Ztaunu2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_status", &mc_truth_Ztaunutau1_status, "mc_truth_Ztaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_status", &mc_truth_Ztaunutau2_status, "mc_truth_Ztaunutau2_status/I", buffersize);

	tree->Branch("mc_truth_t_status", &mc_truth_t_status, "mc_truth_t_status/I", buffersize);
	tree->Branch("mc_truth_tb_status", &mc_truth_tb_status, "mc_truth_tb_status/I", buffersize);
	tree->Branch("mc_truth_tW_status", &mc_truth_tW_status, "mc_truth_tW_status/I", buffersize);
	tree->Branch("mc_truth_tWnu_status", &mc_truth_tWnu_status, "mc_truth_tWnu_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau_status", &mc_truth_tWnutau_status, "mc_truth_tWnutau_status/I", buffersize);
	tree->Branch("mc_truth_tWl_status", &mc_truth_tWl_status, "mc_truth_tWl_status/I", buffersize);
	tree->Branch("mc_truth_tWtau_status", &mc_truth_tWtau_status, "mc_truth_tWtau_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu_status", &mc_truth_tWtaunu_status, "mc_truth_tWtaunu_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau_status", &mc_truth_tWtaunutau_status, "mc_truth_tWtaunutau_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul_status", &mc_truth_tWtaul_status, "mc_truth_tWtaul_status/I", buffersize);
	tree->Branch("mc_truth_tWq1_status", &mc_truth_tWq1_status, "mc_truth_tWq1_status/I", buffersize);
	tree->Branch("mc_truth_tWq2_status", &mc_truth_tWq2_status, "mc_truth_tWq2_status/I", buffersize);

	tree->Branch("mc_truth_j1_status", &mc_truth_j1_status, "mc_truth_j1_status/I", buffersize);
	tree->Branch("mc_truth_j2_status", &mc_truth_j2_status, "mc_truth_j2_status/I", buffersize);
	tree->Branch("mc_truth_j3_status", &mc_truth_j3_status, "mc_truth_j3_status/I", buffersize);
     }

   if( doWrite("mc_truth_thq") )
     {
	tree->Branch("mc_truth_thq_channel", &mc_truth_thq_channel, "mc_truth_thq_channel/I", buffersize);

	if( doWrite("mc_truth_p4") )
	  {	
	     tree->Branch("mc_truth_h0_p4", "TLorentzVector", &mc_truth_h0_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0W1_p4", "TLorentzVector", &mc_truth_h0W1_p4, buffersize);
	     tree->Branch("mc_truth_h0W2_p4", "TLorentzVector", &mc_truth_h0W2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wl1_p4", "TLorentzVector", &mc_truth_h0Wl1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnu1_p4", "TLorentzVector", &mc_truth_h0Wnu1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtau1_p4", "TLorentzVector", &mc_truth_h0Wtau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnutau1_p4", "TLorentzVector", &mc_truth_h0Wnutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaul1_p4", "TLorentzVector", &mc_truth_h0Wtaul1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunu1_p4", "TLorentzVector", &mc_truth_h0Wtaunu1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunutau1_p4", "TLorentzVector", &mc_truth_h0Wtaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wl2_p4", "TLorentzVector", &mc_truth_h0Wl2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnu2_p4", "TLorentzVector", &mc_truth_h0Wnu2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtau2_p4", "TLorentzVector", &mc_truth_h0Wtau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnutau2_p4", "TLorentzVector", &mc_truth_h0Wnutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaul2_p4", "TLorentzVector", &mc_truth_h0Wtaul2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunu2_p4", "TLorentzVector", &mc_truth_h0Wtaunu2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunutau2_p4", "TLorentzVector", &mc_truth_h0Wtaunutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq11_p4", "TLorentzVector", &mc_truth_h0Wq11_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq21_p4", "TLorentzVector", &mc_truth_h0Wq21_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq12_p4", "TLorentzVector", &mc_truth_h0Wq12_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq22_p4", "TLorentzVector", &mc_truth_h0Wq22_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0Z1_p4", "TLorentzVector", &mc_truth_h0Z1_p4, buffersize);
	     tree->Branch("mc_truth_h0Z2_p4", "TLorentzVector", &mc_truth_h0Z2_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl11_p4", "TLorentzVector", &mc_truth_h0Zl11_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl21_p4", "TLorentzVector", &mc_truth_h0Zl21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl12_p4", "TLorentzVector", &mc_truth_h0Zl12_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl22_p4", "TLorentzVector", &mc_truth_h0Zl22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau11_p4", "TLorentzVector", &mc_truth_h0Ztau11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau21_p4", "TLorentzVector", &mc_truth_h0Ztau21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul11_p4", "TLorentzVector", &mc_truth_h0Ztaul11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul21_p4", "TLorentzVector", &mc_truth_h0Ztaul21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu11_p4", "TLorentzVector", &mc_truth_h0Ztaunu11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu21_p4", "TLorentzVector", &mc_truth_h0Ztaunu21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau11_p4", "TLorentzVector", &mc_truth_h0Ztaunutau11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau21_p4", "TLorentzVector", &mc_truth_h0Ztaunutau21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq11_p4", "TLorentzVector", &mc_truth_h0Zq11_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq21_p4", "TLorentzVector", &mc_truth_h0Zq21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq12_p4", "TLorentzVector", &mc_truth_h0Zq12_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq22_p4", "TLorentzVector", &mc_truth_h0Zq22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau12_p4", "TLorentzVector", &mc_truth_h0Ztau12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau22_p4", "TLorentzVector", &mc_truth_h0Ztau22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul12_p4", "TLorentzVector", &mc_truth_h0Ztaul12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul22_p4", "TLorentzVector", &mc_truth_h0Ztaul22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu12_p4", "TLorentzVector", &mc_truth_h0Ztaunu12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu22_p4", "TLorentzVector", &mc_truth_h0Ztaunu22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau12_p4", "TLorentzVector", &mc_truth_h0Ztaunutau12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau22_p4", "TLorentzVector", &mc_truth_h0Ztaunutau22_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq12_p4", "TLorentzVector", &mc_truth_h0Zq12_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq22_p4", "TLorentzVector", &mc_truth_h0Zq22_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu11_p4", "TLorentzVector", &mc_truth_h0Znu11_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu21_p4", "TLorentzVector", &mc_truth_h0Znu21_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu12_p4", "TLorentzVector", &mc_truth_h0Znu12_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu22_p4", "TLorentzVector", &mc_truth_h0Znu22_p4, buffersize);

	     tree->Branch("mc_truth_h0tau1_p4", "TLorentzVector", &mc_truth_h0tau1_p4, buffersize);
	     tree->Branch("mc_truth_h0tau2_p4", "TLorentzVector", &mc_truth_h0tau2_p4, buffersize);
	     tree->Branch("mc_truth_h0taul1_p4", "TLorentzVector", &mc_truth_h0taul1_p4, buffersize);
	     tree->Branch("mc_truth_h0taunutau1_p4", "TLorentzVector", &mc_truth_h0taunutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0taunu1_p4", "TLorentzVector", &mc_truth_h0taunu1_p4, buffersize);
	     tree->Branch("mc_truth_h0taul2_p4", "TLorentzVector", &mc_truth_h0taul2_p4, buffersize);
	     tree->Branch("mc_truth_h0taunutau2_p4", "TLorentzVector", &mc_truth_h0taunutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0taunu2_p4", "TLorentzVector", &mc_truth_h0taunu2_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0b1_p4", "TLorentzVector", &mc_truth_h0b1_p4, buffersize);
	     tree->Branch("mc_truth_h0b2_p4", "TLorentzVector", &mc_truth_h0b2_p4, buffersize);
	     tree->Branch("mc_truth_h0b1_IS_p4", "TLorentzVector", &mc_truth_h0b1_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0b2_IS_p4", "TLorentzVector", &mc_truth_h0b2_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_t_p4", "TLorentzVector", &mc_truth_t_p4, buffersize);
	     tree->Branch("mc_truth_tb_p4", "TLorentzVector", &mc_truth_tb_p4, buffersize);
	     tree->Branch("mc_truth_tb_IS_p4", "TLorentzVector", &mc_truth_tb_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW_p4", "TLorentzVector", &mc_truth_tW_p4, buffersize);
	     tree->Branch("mc_truth_tWnu_p4", "TLorentzVector", &mc_truth_tWnu_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau_p4", "TLorentzVector", &mc_truth_tWnutau_p4, buffersize);
	     tree->Branch("mc_truth_tWl_p4", "TLorentzVector", &mc_truth_tWl_p4, buffersize);
	     tree->Branch("mc_truth_tWtau_p4", "TLorentzVector", &mc_truth_tWtau_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu_p4", "TLorentzVector", &mc_truth_tWtaunu_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau_p4", "TLorentzVector", &mc_truth_tWtaunutau_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul_p4", "TLorentzVector", &mc_truth_tWtaul_p4, buffersize);
	     tree->Branch("mc_truth_tWq1_p4", "TLorentzVector", &mc_truth_tWq1_p4, buffersize);
	     tree->Branch("mc_truth_tWq2_p4", "TLorentzVector", &mc_truth_tWq2_p4, buffersize);
	     
	     tree->Branch("mc_truth_j1_p4", "TLorentzVector", &mc_truth_j1_p4, buffersize);
	     tree->Branch("mc_truth_j2_p4", "TLorentzVector", &mc_truth_j2_p4, buffersize);
	     tree->Branch("mc_truth_j3_p4", "TLorentzVector", &mc_truth_j3_p4, buffersize);
	  }

	tree->Branch("mc_truth_h0_pt", &mc_truth_h0_pt, "mc_truth_h0_pt/F", buffersize);

	tree->Branch("mc_truth_h0W1_pt", &mc_truth_h0W1_pt, "mc_truth_h0W1_pt/F", buffersize);
	tree->Branch("mc_truth_h0W2_pt", &mc_truth_h0W2_pt, "mc_truth_h0W2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_pt", &mc_truth_h0Wl1_pt, "mc_truth_h0Wl1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_pt", &mc_truth_h0Wnu1_pt, "mc_truth_h0Wnu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_pt", &mc_truth_h0Wtau1_pt, "mc_truth_h0Wtau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_pt", &mc_truth_h0Wnutau1_pt, "mc_truth_h0Wnutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_pt", &mc_truth_h0Wtaul1_pt, "mc_truth_h0Wtaul1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_pt", &mc_truth_h0Wtaunu1_pt, "mc_truth_h0Wtaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_pt", &mc_truth_h0Wtaunutau1_pt, "mc_truth_h0Wtaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_pt", &mc_truth_h0Wl2_pt, "mc_truth_h0Wl2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_pt", &mc_truth_h0Wnu2_pt, "mc_truth_h0Wnu2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_pt", &mc_truth_h0Wtau2_pt, "mc_truth_h0Wtau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_pt", &mc_truth_h0Wnutau2_pt, "mc_truth_h0Wnutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_pt", &mc_truth_h0Wtaul2_pt, "mc_truth_h0Wtaul2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_pt", &mc_truth_h0Wtaunu2_pt, "mc_truth_h0Wtaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_pt", &mc_truth_h0Wtaunutau2_pt, "mc_truth_h0Wtaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_pt", &mc_truth_h0Wq11_pt, "mc_truth_h0Wq11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_pt", &mc_truth_h0Wq21_pt, "mc_truth_h0Wq21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_pt", &mc_truth_h0Wq12_pt, "mc_truth_h0Wq12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_pt", &mc_truth_h0Wq22_pt, "mc_truth_h0Wq22_pt/F", buffersize);

	tree->Branch("mc_truth_h0Z1_pt", &mc_truth_h0Z1_pt, "mc_truth_h0Z1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Z2_pt", &mc_truth_h0Z2_pt, "mc_truth_h0Z2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_pt", &mc_truth_h0Zl11_pt, "mc_truth_h0Zl11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_pt", &mc_truth_h0Zl21_pt, "mc_truth_h0Zl21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_pt", &mc_truth_h0Ztau11_pt, "mc_truth_h0Ztau11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_pt", &mc_truth_h0Ztau21_pt, "mc_truth_h0Ztau21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_pt", &mc_truth_h0Ztaul11_pt, "mc_truth_h0Ztaul11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_pt", &mc_truth_h0Ztaul21_pt, "mc_truth_h0Ztaul21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_pt", &mc_truth_h0Ztaunu11_pt, "mc_truth_h0Ztaunu11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_pt", &mc_truth_h0Ztaunu21_pt, "mc_truth_h0Ztaunu21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_pt", &mc_truth_h0Ztaunutau11_pt, "mc_truth_h0Ztaunutau11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_pt", &mc_truth_h0Ztaunutau21_pt, "mc_truth_h0Ztaunutau21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_pt", &mc_truth_h0Zq11_pt, "mc_truth_h0Zq11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_pt", &mc_truth_h0Zq21_pt, "mc_truth_h0Zq21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_pt", &mc_truth_h0Zl12_pt, "mc_truth_h0Zl12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_pt", &mc_truth_h0Zl22_pt, "mc_truth_h0Zl22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_pt", &mc_truth_h0Ztau12_pt, "mc_truth_h0Ztau12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_pt", &mc_truth_h0Ztau22_pt, "mc_truth_h0Ztau22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_pt", &mc_truth_h0Ztaul12_pt, "mc_truth_h0Ztaul12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_pt", &mc_truth_h0Ztaul22_pt, "mc_truth_h0Ztaul22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_pt", &mc_truth_h0Ztaunu12_pt, "mc_truth_h0Ztaunu12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_pt", &mc_truth_h0Ztaunu22_pt, "mc_truth_h0Ztaunu22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_pt", &mc_truth_h0Ztaunutau12_pt, "mc_truth_h0Ztaunutau12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_pt", &mc_truth_h0Ztaunutau22_pt, "mc_truth_h0Ztaunutau22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_pt", &mc_truth_h0Zq12_pt, "mc_truth_h0Zq12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_pt", &mc_truth_h0Zq22_pt, "mc_truth_h0Zq22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_pt", &mc_truth_h0Znu11_pt, "mc_truth_h0Znu11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_pt", &mc_truth_h0Znu21_pt, "mc_truth_h0Znu21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_pt", &mc_truth_h0Znu12_pt, "mc_truth_h0Znu12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_pt", &mc_truth_h0Znu22_pt, "mc_truth_h0Znu22_pt/F", buffersize);

	tree->Branch("mc_truth_h0tau1_pt", &mc_truth_h0tau1_pt, "mc_truth_h0tau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0tau2_pt", &mc_truth_h0tau2_pt, "mc_truth_h0tau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taul1_pt", &mc_truth_h0taul1_pt, "mc_truth_h0taul1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_pt", &mc_truth_h0taunutau1_pt, "mc_truth_h0taunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_pt", &mc_truth_h0taunu1_pt, "mc_truth_h0taunu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taul2_pt", &mc_truth_h0taul2_pt, "mc_truth_h0taul2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_pt", &mc_truth_h0taunutau2_pt, "mc_truth_h0taunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_pt", &mc_truth_h0taunu2_pt, "mc_truth_h0taunu2_pt/F", buffersize);

	tree->Branch("mc_truth_h0b1_pt", &mc_truth_h0b1_pt, "mc_truth_h0b1_pt/F", buffersize);
	tree->Branch("mc_truth_h0b2_pt", &mc_truth_h0b2_pt, "mc_truth_h0b2_pt/F", buffersize);
	tree->Branch("mc_truth_h0b1_IS_pt", &mc_truth_h0b1_IS_pt, "mc_truth_h0b1_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0b2_IS_pt", &mc_truth_h0b2_IS_pt, "mc_truth_h0b2_IS_pt/F", buffersize);
	
	tree->Branch("mc_truth_t_pt", &mc_truth_t_pt, "mc_truth_t_pt/F", buffersize);
	tree->Branch("mc_truth_tb_pt", &mc_truth_tb_pt, "mc_truth_tb_pt/F", buffersize);
	tree->Branch("mc_truth_tb_IS_pt", &mc_truth_tb_IS_pt, "mc_truth_tb_IS_pt/F", buffersize);

	tree->Branch("mc_truth_tW_pt", &mc_truth_tW_pt, "mc_truth_tW_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu_pt", &mc_truth_tWnu_pt, "mc_truth_tWnu_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau_pt", &mc_truth_tWnutau_pt, "mc_truth_tWnutau_pt/F", buffersize);
	tree->Branch("mc_truth_tWl_pt", &mc_truth_tWl_pt, "mc_truth_tWl_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau_pt", &mc_truth_tWtau_pt, "mc_truth_tWtau_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_pt", &mc_truth_tWtaunu_pt, "mc_truth_tWtaunu_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_pt", &mc_truth_tWtaunutau_pt, "mc_truth_tWtaunutau_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul_pt", &mc_truth_tWtaul_pt, "mc_truth_tWtaul_pt/F", buffersize);
	tree->Branch("mc_truth_tWq1_pt", &mc_truth_tWq1_pt, "mc_truth_tWq1_pt/F", buffersize);
	tree->Branch("mc_truth_tWq2_pt", &mc_truth_tWq2_pt, "mc_truth_tWq2_pt/F", buffersize);

	tree->Branch("mc_truth_j1_pt", &mc_truth_j1_pt, "mc_truth_j1_pt/F", buffersize);
	tree->Branch("mc_truth_j2_pt", &mc_truth_j2_pt, "mc_truth_j2_pt/F", buffersize);
	tree->Branch("mc_truth_j3_pt", &mc_truth_j3_pt, "mc_truth_j3_pt/F", buffersize);
		
	tree->Branch("mc_truth_h0_eta", &mc_truth_h0_eta, "mc_truth_h0_eta/F", buffersize);

	tree->Branch("mc_truth_h0W1_eta", &mc_truth_h0W1_eta, "mc_truth_h0W1_eta/F", buffersize);
	tree->Branch("mc_truth_h0W2_eta", &mc_truth_h0W2_eta, "mc_truth_h0W2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_eta", &mc_truth_h0Wl1_eta, "mc_truth_h0Wl1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_eta", &mc_truth_h0Wnu1_eta, "mc_truth_h0Wnu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_eta", &mc_truth_h0Wtau1_eta, "mc_truth_h0Wtau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_eta", &mc_truth_h0Wnutau1_eta, "mc_truth_h0Wnutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_eta", &mc_truth_h0Wtaul1_eta, "mc_truth_h0Wtaul1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_eta", &mc_truth_h0Wtaunu1_eta, "mc_truth_h0Wtaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_eta", &mc_truth_h0Wtaunutau1_eta, "mc_truth_h0Wtaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_eta", &mc_truth_h0Wl2_eta, "mc_truth_h0Wl2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_eta", &mc_truth_h0Wnu2_eta, "mc_truth_h0Wnu2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_eta", &mc_truth_h0Wtau2_eta, "mc_truth_h0Wtau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_eta", &mc_truth_h0Wnutau2_eta, "mc_truth_h0Wnutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_eta", &mc_truth_h0Wtaul2_eta, "mc_truth_h0Wtaul2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_eta", &mc_truth_h0Wtaunu2_eta, "mc_truth_h0Wtaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_eta", &mc_truth_h0Wtaunutau2_eta, "mc_truth_h0Wtaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_eta", &mc_truth_h0Wq11_eta, "mc_truth_h0Wq11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_eta", &mc_truth_h0Wq21_eta, "mc_truth_h0Wq21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_eta", &mc_truth_h0Wq12_eta, "mc_truth_h0Wq12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_eta", &mc_truth_h0Wq22_eta, "mc_truth_h0Wq22_eta/F", buffersize);

	tree->Branch("mc_truth_h0Z1_eta", &mc_truth_h0Z1_eta, "mc_truth_h0Z1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Z2_eta", &mc_truth_h0Z2_eta, "mc_truth_h0Z2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_eta", &mc_truth_h0Zl11_eta, "mc_truth_h0Zl11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_eta", &mc_truth_h0Zl21_eta, "mc_truth_h0Zl21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_eta", &mc_truth_h0Ztau11_eta, "mc_truth_h0Ztau11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_eta", &mc_truth_h0Ztau21_eta, "mc_truth_h0Ztau21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_eta", &mc_truth_h0Ztaul11_eta, "mc_truth_h0Ztaul11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_eta", &mc_truth_h0Ztaul21_eta, "mc_truth_h0Ztaul21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_eta", &mc_truth_h0Ztaunu11_eta, "mc_truth_h0Ztaunu11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_eta", &mc_truth_h0Ztaunu21_eta, "mc_truth_h0Ztaunu21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_eta", &mc_truth_h0Ztaunutau11_eta, "mc_truth_h0Ztaunutau11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_eta", &mc_truth_h0Ztaunutau21_eta, "mc_truth_h0Ztaunutau21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_eta", &mc_truth_h0Zq11_eta, "mc_truth_h0Zq11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_eta", &mc_truth_h0Zq21_eta, "mc_truth_h0Zq21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_eta", &mc_truth_h0Zl12_eta, "mc_truth_h0Zl12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_eta", &mc_truth_h0Zl22_eta, "mc_truth_h0Zl22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_eta", &mc_truth_h0Ztau12_eta, "mc_truth_h0Ztau12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_eta", &mc_truth_h0Ztau22_eta, "mc_truth_h0Ztau22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_eta", &mc_truth_h0Ztaul12_eta, "mc_truth_h0Ztaul12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_eta", &mc_truth_h0Ztaul22_eta, "mc_truth_h0Ztaul22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_eta", &mc_truth_h0Ztaunu12_eta, "mc_truth_h0Ztaunu12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_eta", &mc_truth_h0Ztaunu22_eta, "mc_truth_h0Ztaunu22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_eta", &mc_truth_h0Ztaunutau12_eta, "mc_truth_h0Ztaunutau12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_eta", &mc_truth_h0Ztaunutau22_eta, "mc_truth_h0Ztaunutau22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_eta", &mc_truth_h0Zq12_eta, "mc_truth_h0Zq12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_eta", &mc_truth_h0Zq22_eta, "mc_truth_h0Zq22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_eta", &mc_truth_h0Znu11_eta, "mc_truth_h0Znu11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_eta", &mc_truth_h0Znu21_eta, "mc_truth_h0Znu21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_eta", &mc_truth_h0Znu12_eta, "mc_truth_h0Znu12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_eta", &mc_truth_h0Znu22_eta, "mc_truth_h0Znu22_eta/F", buffersize);

	tree->Branch("mc_truth_h0tau1_eta", &mc_truth_h0tau1_eta, "mc_truth_h0tau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0tau2_eta", &mc_truth_h0tau2_eta, "mc_truth_h0tau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taul1_eta", &mc_truth_h0taul1_eta, "mc_truth_h0taul1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_eta", &mc_truth_h0taunutau1_eta, "mc_truth_h0taunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_eta", &mc_truth_h0taunu1_eta, "mc_truth_h0taunu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taul2_eta", &mc_truth_h0taul2_eta, "mc_truth_h0taul2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_eta", &mc_truth_h0taunutau2_eta, "mc_truth_h0taunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_eta", &mc_truth_h0taunu2_eta, "mc_truth_h0taunu2_eta/F", buffersize);

	tree->Branch("mc_truth_h0b1_eta", &mc_truth_h0b1_eta, "mc_truth_h0b1_eta/F", buffersize);
	tree->Branch("mc_truth_h0b2_eta", &mc_truth_h0b2_eta, "mc_truth_h0b2_eta/F", buffersize);
	tree->Branch("mc_truth_h0b1_IS_eta", &mc_truth_h0b1_IS_eta, "mc_truth_h0b1_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0b2_IS_eta", &mc_truth_h0b2_IS_eta, "mc_truth_h0b2_IS_eta/F", buffersize);
	
	tree->Branch("mc_truth_t_eta", &mc_truth_t_eta, "mc_truth_t_eta/F", buffersize);
	tree->Branch("mc_truth_tb_eta", &mc_truth_tb_eta, "mc_truth_tb_eta/F", buffersize);
	tree->Branch("mc_truth_tb_IS_eta", &mc_truth_tb_IS_eta, "mc_truth_tb_IS_eta/F", buffersize);

	tree->Branch("mc_truth_tW_eta", &mc_truth_tW_eta, "mc_truth_tW_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu_eta", &mc_truth_tWnu_eta, "mc_truth_tWnu_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau_eta", &mc_truth_tWnutau_eta, "mc_truth_tWnutau_eta/F", buffersize);
	tree->Branch("mc_truth_tWl_eta", &mc_truth_tWl_eta, "mc_truth_tWl_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau_eta", &mc_truth_tWtau_eta, "mc_truth_tWtau_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_eta", &mc_truth_tWtaunu_eta, "mc_truth_tWtaunu_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_eta", &mc_truth_tWtaunutau_eta, "mc_truth_tWtaunutau_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul_eta", &mc_truth_tWtaul_eta, "mc_truth_tWtaul_eta/F", buffersize);
	tree->Branch("mc_truth_tWq1_eta", &mc_truth_tWq1_eta, "mc_truth_tWq1_eta/F", buffersize);
	tree->Branch("mc_truth_tWq2_eta", &mc_truth_tWq2_eta, "mc_truth_tWq2_eta/F", buffersize);

	tree->Branch("mc_truth_j1_eta", &mc_truth_j1_eta, "mc_truth_j1_eta/F", buffersize);
	tree->Branch("mc_truth_j2_eta", &mc_truth_j2_eta, "mc_truth_j2_eta/F", buffersize);
	tree->Branch("mc_truth_j3_eta", &mc_truth_j3_eta, "mc_truth_j3_eta/F", buffersize);
		
	tree->Branch("mc_truth_h0_phi", &mc_truth_h0_phi, "mc_truth_h0_phi/F", buffersize);

	tree->Branch("mc_truth_h0W1_phi", &mc_truth_h0W1_phi, "mc_truth_h0W1_phi/F", buffersize);
	tree->Branch("mc_truth_h0W2_phi", &mc_truth_h0W2_phi, "mc_truth_h0W2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_phi", &mc_truth_h0Wl1_phi, "mc_truth_h0Wl1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_phi", &mc_truth_h0Wnu1_phi, "mc_truth_h0Wnu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_phi", &mc_truth_h0Wtau1_phi, "mc_truth_h0Wtau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_phi", &mc_truth_h0Wnutau1_phi, "mc_truth_h0Wnutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_phi", &mc_truth_h0Wtaul1_phi, "mc_truth_h0Wtaul1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_phi", &mc_truth_h0Wtaunu1_phi, "mc_truth_h0Wtaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_phi", &mc_truth_h0Wtaunutau1_phi, "mc_truth_h0Wtaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_phi", &mc_truth_h0Wl2_phi, "mc_truth_h0Wl2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_phi", &mc_truth_h0Wnu2_phi, "mc_truth_h0Wnu2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_phi", &mc_truth_h0Wtau2_phi, "mc_truth_h0Wtau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_phi", &mc_truth_h0Wnutau2_phi, "mc_truth_h0Wnutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_phi", &mc_truth_h0Wtaul2_phi, "mc_truth_h0Wtaul2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_phi", &mc_truth_h0Wtaunu2_phi, "mc_truth_h0Wtaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_phi", &mc_truth_h0Wtaunutau2_phi, "mc_truth_h0Wtaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_phi", &mc_truth_h0Wq11_phi, "mc_truth_h0Wq11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_phi", &mc_truth_h0Wq21_phi, "mc_truth_h0Wq21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_phi", &mc_truth_h0Wq12_phi, "mc_truth_h0Wq12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_phi", &mc_truth_h0Wq22_phi, "mc_truth_h0Wq22_phi/F", buffersize);

	tree->Branch("mc_truth_h0Z1_phi", &mc_truth_h0Z1_phi, "mc_truth_h0Z1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Z2_phi", &mc_truth_h0Z2_phi, "mc_truth_h0Z2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_phi", &mc_truth_h0Zl11_phi, "mc_truth_h0Zl11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_phi", &mc_truth_h0Zl21_phi, "mc_truth_h0Zl21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_phi", &mc_truth_h0Ztau11_phi, "mc_truth_h0Ztau11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_phi", &mc_truth_h0Ztau21_phi, "mc_truth_h0Ztau21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_phi", &mc_truth_h0Ztaul11_phi, "mc_truth_h0Ztaul11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_phi", &mc_truth_h0Ztaul21_phi, "mc_truth_h0Ztaul21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_phi", &mc_truth_h0Ztaunu11_phi, "mc_truth_h0Ztaunu11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_phi", &mc_truth_h0Ztaunu21_phi, "mc_truth_h0Ztaunu21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_phi", &mc_truth_h0Ztaunutau11_phi, "mc_truth_h0Ztaunutau11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_phi", &mc_truth_h0Ztaunutau21_phi, "mc_truth_h0Ztaunutau21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_phi", &mc_truth_h0Zq11_phi, "mc_truth_h0Zq11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_phi", &mc_truth_h0Zq21_phi, "mc_truth_h0Zq21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_phi", &mc_truth_h0Zl12_phi, "mc_truth_h0Zl12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_phi", &mc_truth_h0Zl22_phi, "mc_truth_h0Zl22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_phi", &mc_truth_h0Ztau12_phi, "mc_truth_h0Ztau12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_phi", &mc_truth_h0Ztau22_phi, "mc_truth_h0Ztau22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_phi", &mc_truth_h0Ztaul12_phi, "mc_truth_h0Ztaul12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_phi", &mc_truth_h0Ztaul22_phi, "mc_truth_h0Ztaul22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_phi", &mc_truth_h0Ztaunu12_phi, "mc_truth_h0Ztaunu12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_phi", &mc_truth_h0Ztaunu22_phi, "mc_truth_h0Ztaunu22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_phi", &mc_truth_h0Ztaunutau12_phi, "mc_truth_h0Ztaunutau12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_phi", &mc_truth_h0Ztaunutau22_phi, "mc_truth_h0Ztaunutau22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_phi", &mc_truth_h0Zq12_phi, "mc_truth_h0Zq12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_phi", &mc_truth_h0Zq22_phi, "mc_truth_h0Zq22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_phi", &mc_truth_h0Znu11_phi, "mc_truth_h0Znu11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_phi", &mc_truth_h0Znu21_phi, "mc_truth_h0Znu21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_phi", &mc_truth_h0Znu12_phi, "mc_truth_h0Znu12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_phi", &mc_truth_h0Znu22_phi, "mc_truth_h0Znu22_phi/F", buffersize);

	tree->Branch("mc_truth_h0tau1_phi", &mc_truth_h0tau1_phi, "mc_truth_h0tau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0tau2_phi", &mc_truth_h0tau2_phi, "mc_truth_h0tau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taul1_phi", &mc_truth_h0taul1_phi, "mc_truth_h0taul1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_phi", &mc_truth_h0taunutau1_phi, "mc_truth_h0taunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_phi", &mc_truth_h0taunu1_phi, "mc_truth_h0taunu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taul2_phi", &mc_truth_h0taul2_phi, "mc_truth_h0taul2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_phi", &mc_truth_h0taunutau2_phi, "mc_truth_h0taunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_phi", &mc_truth_h0taunu2_phi, "mc_truth_h0taunu2_phi/F", buffersize);

	tree->Branch("mc_truth_h0b1_phi", &mc_truth_h0b1_phi, "mc_truth_h0b1_phi/F", buffersize);
	tree->Branch("mc_truth_h0b2_phi", &mc_truth_h0b2_phi, "mc_truth_h0b2_phi/F", buffersize);
	tree->Branch("mc_truth_h0b1_IS_phi", &mc_truth_h0b1_IS_phi, "mc_truth_h0b1_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0b2_IS_phi", &mc_truth_h0b2_IS_phi, "mc_truth_h0b2_IS_phi/F", buffersize);
	
	tree->Branch("mc_truth_t_phi", &mc_truth_t_phi, "mc_truth_t_phi/F", buffersize);
	tree->Branch("mc_truth_tb_phi", &mc_truth_tb_phi, "mc_truth_tb_phi/F", buffersize);
	tree->Branch("mc_truth_tb_IS_phi", &mc_truth_tb_IS_phi, "mc_truth_tb_IS_phi/F", buffersize);

	tree->Branch("mc_truth_tW_phi", &mc_truth_tW_phi, "mc_truth_tW_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu_phi", &mc_truth_tWnu_phi, "mc_truth_tWnu_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau_phi", &mc_truth_tWnutau_phi, "mc_truth_tWnutau_phi/F", buffersize);
	tree->Branch("mc_truth_tWl_phi", &mc_truth_tWl_phi, "mc_truth_tWl_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau_phi", &mc_truth_tWtau_phi, "mc_truth_tWtau_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_phi", &mc_truth_tWtaunu_phi, "mc_truth_tWtaunu_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_phi", &mc_truth_tWtaunutau_phi, "mc_truth_tWtaunutau_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul_phi", &mc_truth_tWtaul_phi, "mc_truth_tWtaul_phi/F", buffersize);
	tree->Branch("mc_truth_tWq1_phi", &mc_truth_tWq1_phi, "mc_truth_tWq1_phi/F", buffersize);
	tree->Branch("mc_truth_tWq2_phi", &mc_truth_tWq2_phi, "mc_truth_tWq2_phi/F", buffersize);

	tree->Branch("mc_truth_j1_phi", &mc_truth_j1_phi, "mc_truth_j1_phi/F", buffersize);
	tree->Branch("mc_truth_j2_phi", &mc_truth_j2_phi, "mc_truth_j2_phi/F", buffersize);
	tree->Branch("mc_truth_j3_phi", &mc_truth_j3_phi, "mc_truth_j3_phi/F", buffersize);
	
	tree->Branch("mc_truth_h0_E", &mc_truth_h0_E, "mc_truth_h0_E/F", buffersize);

	tree->Branch("mc_truth_h0W1_E", &mc_truth_h0W1_E, "mc_truth_h0W1_E/F", buffersize);
	tree->Branch("mc_truth_h0W2_E", &mc_truth_h0W2_E, "mc_truth_h0W2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_E", &mc_truth_h0Wl1_E, "mc_truth_h0Wl1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_E", &mc_truth_h0Wnu1_E, "mc_truth_h0Wnu1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_E", &mc_truth_h0Wtau1_E, "mc_truth_h0Wtau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_E", &mc_truth_h0Wnutau1_E, "mc_truth_h0Wnutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_E", &mc_truth_h0Wtaul1_E, "mc_truth_h0Wtaul1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_E", &mc_truth_h0Wtaunu1_E, "mc_truth_h0Wtaunu1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_E", &mc_truth_h0Wtaunutau1_E, "mc_truth_h0Wtaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_E", &mc_truth_h0Wl2_E, "mc_truth_h0Wl2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_E", &mc_truth_h0Wnu2_E, "mc_truth_h0Wnu2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_E", &mc_truth_h0Wtau2_E, "mc_truth_h0Wtau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_E", &mc_truth_h0Wnutau2_E, "mc_truth_h0Wnutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_E", &mc_truth_h0Wtaul2_E, "mc_truth_h0Wtaul2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_E", &mc_truth_h0Wtaunu2_E, "mc_truth_h0Wtaunu2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_E", &mc_truth_h0Wtaunutau2_E, "mc_truth_h0Wtaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_E", &mc_truth_h0Wq11_E, "mc_truth_h0Wq11_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_E", &mc_truth_h0Wq21_E, "mc_truth_h0Wq21_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_E", &mc_truth_h0Wq12_E, "mc_truth_h0Wq12_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_E", &mc_truth_h0Wq22_E, "mc_truth_h0Wq22_E/F", buffersize);

	tree->Branch("mc_truth_h0Z1_E", &mc_truth_h0Z1_E, "mc_truth_h0Z1_E/F", buffersize);
	tree->Branch("mc_truth_h0Z2_E", &mc_truth_h0Z2_E, "mc_truth_h0Z2_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_E", &mc_truth_h0Zl11_E, "mc_truth_h0Zl11_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_E", &mc_truth_h0Zl21_E, "mc_truth_h0Zl21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_E", &mc_truth_h0Ztau11_E, "mc_truth_h0Ztau11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_E", &mc_truth_h0Ztau21_E, "mc_truth_h0Ztau21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_E", &mc_truth_h0Ztaul11_E, "mc_truth_h0Ztaul11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_E", &mc_truth_h0Ztaul21_E, "mc_truth_h0Ztaul21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_E", &mc_truth_h0Ztaunu11_E, "mc_truth_h0Ztaunu11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_E", &mc_truth_h0Ztaunu21_E, "mc_truth_h0Ztaunu21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_E", &mc_truth_h0Ztaunutau11_E, "mc_truth_h0Ztaunutau11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_E", &mc_truth_h0Ztaunutau21_E, "mc_truth_h0Ztaunutau21_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_E", &mc_truth_h0Zq11_E, "mc_truth_h0Zq11_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_E", &mc_truth_h0Zq21_E, "mc_truth_h0Zq21_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_E", &mc_truth_h0Zl12_E, "mc_truth_h0Zl12_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_E", &mc_truth_h0Zl22_E, "mc_truth_h0Zl22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_E", &mc_truth_h0Ztau12_E, "mc_truth_h0Ztau12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_E", &mc_truth_h0Ztau22_E, "mc_truth_h0Ztau22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_E", &mc_truth_h0Ztaul12_E, "mc_truth_h0Ztaul12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_E", &mc_truth_h0Ztaul22_E, "mc_truth_h0Ztaul22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_E", &mc_truth_h0Ztaunu12_E, "mc_truth_h0Ztaunu12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_E", &mc_truth_h0Ztaunu22_E, "mc_truth_h0Ztaunu22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_E", &mc_truth_h0Ztaunutau12_E, "mc_truth_h0Ztaunutau12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_E", &mc_truth_h0Ztaunutau22_E, "mc_truth_h0Ztaunutau22_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_E", &mc_truth_h0Zq12_E, "mc_truth_h0Zq12_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_E", &mc_truth_h0Zq22_E, "mc_truth_h0Zq22_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_E", &mc_truth_h0Znu11_E, "mc_truth_h0Znu11_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_E", &mc_truth_h0Znu21_E, "mc_truth_h0Znu21_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_E", &mc_truth_h0Znu12_E, "mc_truth_h0Znu12_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_E", &mc_truth_h0Znu22_E, "mc_truth_h0Znu22_E/F", buffersize);

	tree->Branch("mc_truth_h0tau1_E", &mc_truth_h0tau1_E, "mc_truth_h0tau1_E/F", buffersize);
	tree->Branch("mc_truth_h0tau2_E", &mc_truth_h0tau2_E, "mc_truth_h0tau2_E/F", buffersize);
	tree->Branch("mc_truth_h0taul1_E", &mc_truth_h0taul1_E, "mc_truth_h0taul1_E/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_E", &mc_truth_h0taunutau1_E, "mc_truth_h0taunutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_E", &mc_truth_h0taunu1_E, "mc_truth_h0taunu1_E/F", buffersize);
	tree->Branch("mc_truth_h0taul2_E", &mc_truth_h0taul2_E, "mc_truth_h0taul2_E/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_E", &mc_truth_h0taunutau2_E, "mc_truth_h0taunutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_E", &mc_truth_h0taunu2_E, "mc_truth_h0taunu2_E/F", buffersize);

	tree->Branch("mc_truth_h0b1_E", &mc_truth_h0b1_E, "mc_truth_h0b1_E/F", buffersize);
	tree->Branch("mc_truth_h0b2_E", &mc_truth_h0b2_E, "mc_truth_h0b2_E/F", buffersize);
	tree->Branch("mc_truth_h0b1_IS_E", &mc_truth_h0b1_IS_E, "mc_truth_h0b1_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0b2_IS_E", &mc_truth_h0b2_IS_E, "mc_truth_h0b2_IS_E/F", buffersize);
	
	tree->Branch("mc_truth_t_E", &mc_truth_t_E, "mc_truth_t_E/F", buffersize);
	tree->Branch("mc_truth_tb_E", &mc_truth_tb_E, "mc_truth_tb_E/F", buffersize);
	tree->Branch("mc_truth_tb_IS_E", &mc_truth_tb_IS_E, "mc_truth_tb_IS_E/F", buffersize);

	tree->Branch("mc_truth_tW_E", &mc_truth_tW_E, "mc_truth_tW_E/F", buffersize);
	tree->Branch("mc_truth_tWnu_E", &mc_truth_tWnu_E, "mc_truth_tWnu_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau_E", &mc_truth_tWnutau_E, "mc_truth_tWnutau_E/F", buffersize);
	tree->Branch("mc_truth_tWl_E", &mc_truth_tWl_E, "mc_truth_tWl_E/F", buffersize);
	tree->Branch("mc_truth_tWtau_E", &mc_truth_tWtau_E, "mc_truth_tWtau_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_E", &mc_truth_tWtaunu_E, "mc_truth_tWtaunu_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_E", &mc_truth_tWtaunutau_E, "mc_truth_tWtaunutau_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul_E", &mc_truth_tWtaul_E, "mc_truth_tWtaul_E/F", buffersize);
	tree->Branch("mc_truth_tWq1_E", &mc_truth_tWq1_E, "mc_truth_tWq1_E/F", buffersize);
	tree->Branch("mc_truth_tWq2_E", &mc_truth_tWq2_E, "mc_truth_tWq2_E/F", buffersize);

	tree->Branch("mc_truth_j1_E", &mc_truth_j1_E, "mc_truth_j1_E/F", buffersize);
	tree->Branch("mc_truth_j2_E", &mc_truth_j2_E, "mc_truth_j2_E/F", buffersize);
	tree->Branch("mc_truth_j3_E", &mc_truth_j3_E, "mc_truth_j3_E/F", buffersize);
				
	tree->Branch("mc_truth_h0_id", &mc_truth_h0_id, "mc_truth_h0_id/I", buffersize);

	tree->Branch("mc_truth_h0W1_id", &mc_truth_h0W1_id, "mc_truth_h0W1_id/I", buffersize);
	tree->Branch("mc_truth_h0W2_id", &mc_truth_h0W2_id, "mc_truth_h0W2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wl1_id", &mc_truth_h0Wl1_id, "mc_truth_h0Wl1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnu1_id", &mc_truth_h0Wnu1_id, "mc_truth_h0Wnu1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtau1_id", &mc_truth_h0Wtau1_id, "mc_truth_h0Wtau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_id", &mc_truth_h0Wnutau1_id, "mc_truth_h0Wnutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_id", &mc_truth_h0Wtaul1_id, "mc_truth_h0Wtaul1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_id", &mc_truth_h0Wtaunu1_id, "mc_truth_h0Wtaunu1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_id", &mc_truth_h0Wtaunutau1_id, "mc_truth_h0Wtaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wl2_id", &mc_truth_h0Wl2_id, "mc_truth_h0Wl2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnu2_id", &mc_truth_h0Wnu2_id, "mc_truth_h0Wnu2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtau2_id", &mc_truth_h0Wtau2_id, "mc_truth_h0Wtau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_id", &mc_truth_h0Wnutau2_id, "mc_truth_h0Wnutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_id", &mc_truth_h0Wtaul2_id, "mc_truth_h0Wtaul2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_id", &mc_truth_h0Wtaunu2_id, "mc_truth_h0Wtaunu2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_id", &mc_truth_h0Wtaunutau2_id, "mc_truth_h0Wtaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_id", &mc_truth_h0Wq11_id, "mc_truth_h0Wq11_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_id", &mc_truth_h0Wq21_id, "mc_truth_h0Wq21_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_id", &mc_truth_h0Wq12_id, "mc_truth_h0Wq12_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_id", &mc_truth_h0Wq22_id, "mc_truth_h0Wq22_id/I", buffersize);

	tree->Branch("mc_truth_h0Z1_id", &mc_truth_h0Z1_id, "mc_truth_h0Z1_id/I", buffersize);
	tree->Branch("mc_truth_h0Z2_id", &mc_truth_h0Z2_id, "mc_truth_h0Z2_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl11_id", &mc_truth_h0Zl11_id, "mc_truth_h0Zl11_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl21_id", &mc_truth_h0Zl21_id, "mc_truth_h0Zl21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau11_id", &mc_truth_h0Ztau11_id, "mc_truth_h0Ztau11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau21_id", &mc_truth_h0Ztau21_id, "mc_truth_h0Ztau21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_id", &mc_truth_h0Ztaul11_id, "mc_truth_h0Ztaul11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_id", &mc_truth_h0Ztaul21_id, "mc_truth_h0Ztaul21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_id", &mc_truth_h0Ztaunu11_id, "mc_truth_h0Ztaunu11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_id", &mc_truth_h0Ztaunu21_id, "mc_truth_h0Ztaunu21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_id", &mc_truth_h0Ztaunutau11_id, "mc_truth_h0Ztaunutau11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_id", &mc_truth_h0Ztaunutau21_id, "mc_truth_h0Ztaunutau21_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_id", &mc_truth_h0Zq11_id, "mc_truth_h0Zq11_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_id", &mc_truth_h0Zq21_id, "mc_truth_h0Zq21_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl12_id", &mc_truth_h0Zl12_id, "mc_truth_h0Zl12_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl22_id", &mc_truth_h0Zl22_id, "mc_truth_h0Zl22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau12_id", &mc_truth_h0Ztau12_id, "mc_truth_h0Ztau12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau22_id", &mc_truth_h0Ztau22_id, "mc_truth_h0Ztau22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_id", &mc_truth_h0Ztaul12_id, "mc_truth_h0Ztaul12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_id", &mc_truth_h0Ztaul22_id, "mc_truth_h0Ztaul22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_id", &mc_truth_h0Ztaunu12_id, "mc_truth_h0Ztaunu12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_id", &mc_truth_h0Ztaunu22_id, "mc_truth_h0Ztaunu22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_id", &mc_truth_h0Ztaunutau12_id, "mc_truth_h0Ztaunutau12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_id", &mc_truth_h0Ztaunutau22_id, "mc_truth_h0Ztaunutau22_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_id", &mc_truth_h0Zq12_id, "mc_truth_h0Zq12_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_id", &mc_truth_h0Zq22_id, "mc_truth_h0Zq22_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu11_id", &mc_truth_h0Znu11_id, "mc_truth_h0Znu11_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu21_id", &mc_truth_h0Znu21_id, "mc_truth_h0Znu21_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu12_id", &mc_truth_h0Znu12_id, "mc_truth_h0Znu12_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu22_id", &mc_truth_h0Znu22_id, "mc_truth_h0Znu22_id/I", buffersize);

	tree->Branch("mc_truth_h0tau1_id", &mc_truth_h0tau1_id, "mc_truth_h0tau1_id/I", buffersize);
	tree->Branch("mc_truth_h0tau2_id", &mc_truth_h0tau2_id, "mc_truth_h0tau2_id/I", buffersize);
	tree->Branch("mc_truth_h0taul1_id", &mc_truth_h0taul1_id, "mc_truth_h0taul1_id/I", buffersize);
	tree->Branch("mc_truth_h0taunutau1_id", &mc_truth_h0taunutau1_id, "mc_truth_h0taunutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0taunu1_id", &mc_truth_h0taunu1_id, "mc_truth_h0taunu1_id/I", buffersize);
	tree->Branch("mc_truth_h0taul2_id", &mc_truth_h0taul2_id, "mc_truth_h0taul2_id/I", buffersize);
	tree->Branch("mc_truth_h0taunutau2_id", &mc_truth_h0taunutau2_id, "mc_truth_h0taunutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0taunu2_id", &mc_truth_h0taunu2_id, "mc_truth_h0taunu2_id/I", buffersize);

	tree->Branch("mc_truth_h0b1_id", &mc_truth_h0b1_id, "mc_truth_h0b1_id/I", buffersize);
	tree->Branch("mc_truth_h0b2_id", &mc_truth_h0b2_id, "mc_truth_h0b2_id/I", buffersize);
	tree->Branch("mc_truth_h0b1_IS_id", &mc_truth_h0b1_IS_id, "mc_truth_h0b1_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0b2_IS_id", &mc_truth_h0b2_IS_id, "mc_truth_h0b2_IS_id/I", buffersize);
	
	tree->Branch("mc_truth_t_id", &mc_truth_t_id, "mc_truth_t_id/I", buffersize);
	tree->Branch("mc_truth_tb_id", &mc_truth_tb_id, "mc_truth_tb_id/I", buffersize);
	tree->Branch("mc_truth_tb_IS_id", &mc_truth_tb_IS_id, "mc_truth_tb_IS_id/I", buffersize);

	tree->Branch("mc_truth_tW_id", &mc_truth_tW_id, "mc_truth_tW_id/I", buffersize);
	tree->Branch("mc_truth_tWnu_id", &mc_truth_tWnu_id, "mc_truth_tWnu_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau_id", &mc_truth_tWnutau_id, "mc_truth_tWnutau_id/I", buffersize);
	tree->Branch("mc_truth_tWl_id", &mc_truth_tWl_id, "mc_truth_tWl_id/I", buffersize);
	tree->Branch("mc_truth_tWtau_id", &mc_truth_tWtau_id, "mc_truth_tWtau_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu_id", &mc_truth_tWtaunu_id, "mc_truth_tWtaunu_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau_id", &mc_truth_tWtaunutau_id, "mc_truth_tWtaunutau_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul_id", &mc_truth_tWtaul_id, "mc_truth_tWtaul_id/I", buffersize);
	tree->Branch("mc_truth_tWq1_id", &mc_truth_tWq1_id, "mc_truth_tWq1_id/I", buffersize);
	tree->Branch("mc_truth_tWq2_id", &mc_truth_tWq2_id, "mc_truth_tWq2_id/I", buffersize);

	tree->Branch("mc_truth_j1_id", &mc_truth_j1_id, "mc_truth_j1_id/I", buffersize);
	tree->Branch("mc_truth_j2_id", &mc_truth_j2_id, "mc_truth_j2_id/I", buffersize);
	tree->Branch("mc_truth_j3_id", &mc_truth_j3_id, "mc_truth_j3_id/I", buffersize);

	tree->Branch("mc_truth_h0_status", &mc_truth_h0_status, "mc_truth_h0_status/I", buffersize);

	tree->Branch("mc_truth_h0W1_status", &mc_truth_h0W1_status, "mc_truth_h0W1_status/I", buffersize);
	tree->Branch("mc_truth_h0W2_status", &mc_truth_h0W2_status, "mc_truth_h0W2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wl1_status", &mc_truth_h0Wl1_status, "mc_truth_h0Wl1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnu1_status", &mc_truth_h0Wnu1_status, "mc_truth_h0Wnu1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtau1_status", &mc_truth_h0Wtau1_status, "mc_truth_h0Wtau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_status", &mc_truth_h0Wnutau1_status, "mc_truth_h0Wnutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_status", &mc_truth_h0Wtaul1_status, "mc_truth_h0Wtaul1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_status", &mc_truth_h0Wtaunu1_status, "mc_truth_h0Wtaunu1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_status", &mc_truth_h0Wtaunutau1_status, "mc_truth_h0Wtaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wl2_status", &mc_truth_h0Wl2_status, "mc_truth_h0Wl2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnu2_status", &mc_truth_h0Wnu2_status, "mc_truth_h0Wnu2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtau2_status", &mc_truth_h0Wtau2_status, "mc_truth_h0Wtau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_status", &mc_truth_h0Wnutau2_status, "mc_truth_h0Wnutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_status", &mc_truth_h0Wtaul2_status, "mc_truth_h0Wtaul2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_status", &mc_truth_h0Wtaunu2_status, "mc_truth_h0Wtaunu2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_status", &mc_truth_h0Wtaunutau2_status, "mc_truth_h0Wtaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_status", &mc_truth_h0Wq11_status, "mc_truth_h0Wq11_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_status", &mc_truth_h0Wq21_status, "mc_truth_h0Wq21_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_status", &mc_truth_h0Wq12_status, "mc_truth_h0Wq12_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_status", &mc_truth_h0Wq22_status, "mc_truth_h0Wq22_status/I", buffersize);

	tree->Branch("mc_truth_h0Z1_status", &mc_truth_h0Z1_status, "mc_truth_h0Z1_status/I", buffersize);
	tree->Branch("mc_truth_h0Z2_status", &mc_truth_h0Z2_status, "mc_truth_h0Z2_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl11_status", &mc_truth_h0Zl11_status, "mc_truth_h0Zl11_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl21_status", &mc_truth_h0Zl21_status, "mc_truth_h0Zl21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau11_status", &mc_truth_h0Ztau11_status, "mc_truth_h0Ztau11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau21_status", &mc_truth_h0Ztau21_status, "mc_truth_h0Ztau21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_status", &mc_truth_h0Ztaul11_status, "mc_truth_h0Ztaul11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_status", &mc_truth_h0Ztaul21_status, "mc_truth_h0Ztaul21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_status", &mc_truth_h0Ztaunu11_status, "mc_truth_h0Ztaunu11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_status", &mc_truth_h0Ztaunu21_status, "mc_truth_h0Ztaunu21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_status", &mc_truth_h0Ztaunutau11_status, "mc_truth_h0Ztaunutau11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_status", &mc_truth_h0Ztaunutau21_status, "mc_truth_h0Ztaunutau21_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_status", &mc_truth_h0Zq11_status, "mc_truth_h0Zq11_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_status", &mc_truth_h0Zq21_status, "mc_truth_h0Zq21_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl12_status", &mc_truth_h0Zl12_status, "mc_truth_h0Zl12_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl22_status", &mc_truth_h0Zl22_status, "mc_truth_h0Zl22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau12_status", &mc_truth_h0Ztau12_status, "mc_truth_h0Ztau12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau22_status", &mc_truth_h0Ztau22_status, "mc_truth_h0Ztau22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_status", &mc_truth_h0Ztaul12_status, "mc_truth_h0Ztaul12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_status", &mc_truth_h0Ztaul22_status, "mc_truth_h0Ztaul22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_status", &mc_truth_h0Ztaunu12_status, "mc_truth_h0Ztaunu12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_status", &mc_truth_h0Ztaunu22_status, "mc_truth_h0Ztaunu22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_status", &mc_truth_h0Ztaunutau12_status, "mc_truth_h0Ztaunutau12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_status", &mc_truth_h0Ztaunutau22_status, "mc_truth_h0Ztaunutau22_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_status", &mc_truth_h0Zq12_status, "mc_truth_h0Zq12_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_status", &mc_truth_h0Zq22_status, "mc_truth_h0Zq22_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu11_status", &mc_truth_h0Znu11_status, "mc_truth_h0Znu11_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu21_status", &mc_truth_h0Znu21_status, "mc_truth_h0Znu21_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu12_status", &mc_truth_h0Znu12_status, "mc_truth_h0Znu12_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu22_status", &mc_truth_h0Znu22_status, "mc_truth_h0Znu22_status/I", buffersize);

	tree->Branch("mc_truth_h0tau1_status", &mc_truth_h0tau1_status, "mc_truth_h0tau1_status/I", buffersize);
	tree->Branch("mc_truth_h0tau2_status", &mc_truth_h0tau2_status, "mc_truth_h0tau2_status/I", buffersize);
	tree->Branch("mc_truth_h0taul1_status", &mc_truth_h0taul1_status, "mc_truth_h0taul1_status/I", buffersize);
	tree->Branch("mc_truth_h0taunutau1_status", &mc_truth_h0taunutau1_status, "mc_truth_h0taunutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0taunu1_status", &mc_truth_h0taunu1_status, "mc_truth_h0taunu1_status/I", buffersize);
	tree->Branch("mc_truth_h0taul2_status", &mc_truth_h0taul2_status, "mc_truth_h0taul2_status/I", buffersize);
	tree->Branch("mc_truth_h0taunutau2_status", &mc_truth_h0taunutau2_status, "mc_truth_h0taunutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0taunu2_status", &mc_truth_h0taunu2_status, "mc_truth_h0taunu2_status/I", buffersize);

	tree->Branch("mc_truth_h0b1_status", &mc_truth_h0b1_status, "mc_truth_h0b1_status/I", buffersize);
	tree->Branch("mc_truth_h0b2_status", &mc_truth_h0b2_status, "mc_truth_h0b2_status/I", buffersize);
	tree->Branch("mc_truth_h0b1_IS_status", &mc_truth_h0b1_IS_status, "mc_truth_h0b1_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0b2_IS_status", &mc_truth_h0b2_IS_status, "mc_truth_h0b2_IS_status/I", buffersize);
	
	tree->Branch("mc_truth_t_status", &mc_truth_t_status, "mc_truth_t_status/I", buffersize);
	tree->Branch("mc_truth_tb_status", &mc_truth_tb_status, "mc_truth_tb_status/I", buffersize);
	tree->Branch("mc_truth_tb_IS_status", &mc_truth_tb_IS_status, "mc_truth_tb_IS_status/I", buffersize);

	tree->Branch("mc_truth_tW_status", &mc_truth_tW_status, "mc_truth_tW_status/I", buffersize);
	tree->Branch("mc_truth_tWnu_status", &mc_truth_tWnu_status, "mc_truth_tWnu_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau_status", &mc_truth_tWnutau_status, "mc_truth_tWnutau_status/I", buffersize);
	tree->Branch("mc_truth_tWl_status", &mc_truth_tWl_status, "mc_truth_tWl_status/I", buffersize);
	tree->Branch("mc_truth_tWtau_status", &mc_truth_tWtau_status, "mc_truth_tWtau_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu_status", &mc_truth_tWtaunu_status, "mc_truth_tWtaunu_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau_status", &mc_truth_tWtaunutau_status, "mc_truth_tWtaunutau_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul_status", &mc_truth_tWtaul_status, "mc_truth_tWtaul_status/I", buffersize);
	tree->Branch("mc_truth_tWq1_status", &mc_truth_tWq1_status, "mc_truth_tWq1_status/I", buffersize);
	tree->Branch("mc_truth_tWq2_status", &mc_truth_tWq2_status, "mc_truth_tWq2_status/I", buffersize);

	tree->Branch("mc_truth_j1_status", &mc_truth_j1_status, "mc_truth_j1_status/I", buffersize);
	tree->Branch("mc_truth_j2_status", &mc_truth_j2_status, "mc_truth_j2_status/I", buffersize);
	tree->Branch("mc_truth_j3_status", &mc_truth_j3_status, "mc_truth_j3_status/I", buffersize);
     }
   
   if( doWrite("gen_all") )
     {
	tree->Branch("gen_n", &gen_n, "gen_n/I", buffersize);
	tree->Branch("gen_pt", "std::vector<float>", &gen_pt, buffersize);
	tree->Branch("gen_eta", "std::vector<float>", &gen_eta, buffersize);
	tree->Branch("gen_phi", "std::vector<float>", &gen_phi, buffersize);
	tree->Branch("gen_m", "std::vector<float>", &gen_m, buffersize);
	tree->Branch("gen_E", "std::vector<float>", &gen_E, buffersize);
	tree->Branch("gen_status", "std::vector<int>", &gen_status, buffersize);
	tree->Branch("gen_id", "std::vector<int>", &gen_id, buffersize);
	tree->Branch("gen_charge", "std::vector<int>", &gen_charge, buffersize);
	tree->Branch("gen_index", "std::vector<int>", &gen_index, buffersize);
	tree->Branch("gen_mother_index", "std::vector<int>", &gen_mother_index, buffersize);
	tree->Branch("gen_daughter_n", "std::vector<int>", &gen_daughter_n, buffersize);
	tree->Branch("gen_daughter_index", "std::vector<std::vector<int> >", &gen_daughter_index, buffersize);
     }

  tree->Branch("n_presel_jets",     &n_presel_jets,     "n_presel_jets/I",      buffersize);
  tree->Branch("n_presel_btag",     &n_presel_btag,     "n_presel_btag/I",      buffersize);
  tree->Branch("n_presel_electron", &n_presel_electron, "n_presel_electron/I",  buffersize);
  tree->Branch("n_presel_muon",     &n_presel_muon,     "n_presel_muon/I",      buffersize);
  tree->Branch("n_presel_tau",      &n_presel_tau,      "n_presel_tau/I",       buffersize);

}

bool FlatTree::doWrite(const std::string& name)
{
  std::map<std::string,bool>::iterator it = conf.find(name);
  if( it != conf.end() )
    return it->second;
  return 0;
}
