#ifndef FLATTREE_H
#define FLATTREE_H

#include <TTree.h>
#include <TLorentzVector.h>
#include <string>
#include <iostream>
#include <vector>

#include <boost/any.hpp>
using boost::any_cast;
typedef std::map<std::string, std::map<std::string, boost::any> >  t_map;

#define DEFVAL -666

#include "IPHCFlatTree/FlatTreeProducer/interface/Helper.hh"

class FlatTree
{
 public:

   FlatTree(TTree* _tree) {tree = _tree;};
   TTree* tree;

   int n_presel_jets;
   int n_presel_btag;
   int n_presel_electron;
   int n_presel_muon;
   int n_presel_tau;

   std::map<std::string, bool> conf;
   t_map keep_conf;

   void Init();
   void CreateBranches(int buffersize);
   bool doWrite(const std::string& name);

   int ev_run;
   int ev_id;
   int ev_lumi;
   float ev_rho;

   float met_pt;
   float met_phi;
   float met_sumet;
   
   float metPuppi_pt;
   float metPuppi_phi;
   float metPuppi_sumet;
   
   float pv_x;
   float pv_y;
   float pv_z;
   
   int pv_ndof;
   float pv_rho;
   int pv_isFake;

   float mc_weight;
   int mc_id;
   int mc_f1;
   int mc_f2;
   float mc_x1;
   float mc_x2;
   float mc_scale;
   float mc_ptHat;

   int mc_pu_intime_NumInt;
   int mc_pu_trueNumInt;
   int mc_pu_before_npu;
   int mc_pu_after_npu;

   int mc_pu_Npvi;
   std::vector<int> mc_pu_Nzpositions;
   std::vector<int> mc_pu_BunchCrossing;
   std::vector<std::vector<float> > mc_pu_zpositions;
   std::vector<std::vector<float> > mc_pu_sumpT_lowpT;
   std::vector<std::vector<float> > mc_pu_sumpT_highpT;
   std::vector<std::vector<int> > mc_pu_ntrks_lowpT;
   std::vector<std::vector<int> > mc_pu_ntrks_highpT;

   // Trigger
   
   int                       trigger_n;
   std::vector<int>          trigger;
   std::vector<std::string>  trigger_name;
   std::vector<bool>         trigger_pass;
   std::vector<int>          trigger_prescale;

     // trigger object general informations
   int                       triggerobject_n;
   std::vector<float>        triggerobject_pt;
   std::vector<float>        triggerobject_eta;
   std::vector<float>        triggerobject_phi;

   std::vector<std::string>  triggerobject_collection;

     // filter Ids...
   std::vector<int>     triggerobject_filterIds_n;
   std::vector<int>     triggerobject_filterIds;

   std::vector<bool>    triggerobject_isTriggerL1Mu;            //-81
   std::vector<bool>    triggerobject_isTriggerL1NoIsoEG;       //-82
   std::vector<bool>    triggerobject_isTriggerL1IsoEG;
   std::vector<bool>    triggerobject_isTriggerL1CenJet;
   std::vector<bool>    triggerobject_isTriggerL1ForJet;
   std::vector<bool>    triggerobject_isTriggerL1TauJet;
   std::vector<bool>    triggerobject_isTriggerL1ETM;
   std::vector<bool>    triggerobject_isTriggerL1ETT;
   std::vector<bool>    triggerobject_isTriggerL1HTT;
   std::vector<bool>    triggerobject_isTriggerL1HTM;
   std::vector<bool>    triggerobject_isTriggerL1JetCounts;
   std::vector<bool>    triggerobject_isTriggerL1HfBitCounts;
   std::vector<bool>    triggerobject_isTriggerL1HfRingEtSums;
   std::vector<bool>    triggerobject_isTriggerL1TechTrig;
   std::vector<bool>    triggerobject_isTriggerL1Castor;
   std::vector<bool>    triggerobject_isTriggerL1BPTX;
   std::vector<bool>    triggerobject_isTriggerL1GtExternal;

   std::vector<bool>    triggerobject_isHLT_TriggerPhoton;      //+81
   std::vector<bool>    triggerobject_isHLT_TriggerElectron;    //+82
   std::vector<bool>    triggerobject_isHLT_TriggerMuon;
   std::vector<bool>    triggerobject_isHLT_TriggerTau;
   std::vector<bool>    triggerobject_isHLT_TriggerJet;
   std::vector<bool>    triggerobject_isHLT_TriggerBJet;
   std::vector<bool>    triggerobject_isHLT_TriggerMET;
   std::vector<bool>    triggerobject_isHLT_TriggerTET;
   std::vector<bool>    triggerobject_isHLT_TriggerTHT;
   std::vector<bool>    triggerobject_isHLT_TriggerMHT;
   std::vector<bool>    triggerobject_isHLT_TriggerTrack;
   std::vector<bool>    triggerobject_isHLT_TriggerCluster;
   std::vector<bool>    triggerobject_isHLT_TriggerMETSig;
   std::vector<bool>    triggerobject_isHLT_TriggerELongit;
   std::vector<bool>    triggerobject_isHLT_TriggerMHTSig;
   std::vector<bool>    triggerobject_isHLT_TriggerHLongit;

     // filters label...
   std::vector<int>           triggerobject_filterLabels_n;
   std::vector<std::string>   triggerobject_filterLabels;

     // paths names and status
   std::vector<int>           triggerobject_pathNamesAll_n;
   std::vector<std::string>   triggerobject_pathNamesAll;
   std::vector<bool>          triggerobject_pathNamesAll_isBoth;
   std::vector<bool>          triggerobject_pathNamesAll_isL3;
   std::vector<bool>          triggerobject_pathNamesAll_isLF;
   std::vector<bool>          triggerobject_pathNamesAll_isNone;


   // Electrons

   int el_n;
   std::vector<float> el_pt;
   std::vector<float> el_eta;
   std::vector<float> el_phi;
   std::vector<float> el_m;
   std::vector<float> el_E;
   std::vector<int> el_id;
   std::vector<int> el_charge;

   std::vector<float> el_scleta;
   std::vector<int> el_passConversionVeto;
   std::vector<int> el_isGsfCtfScPixChargeConsistent;
   std::vector<float> el_dB3D;
   std::vector<float> el_edB3D;
   std::vector<float> el_dB;
   std::vector<float> el_edB;

   std::vector<float> el_miniIso;

   std::vector<float> el_neutralHadronIso;
   std::vector<float> el_chargedHadronIso;
   std::vector<float> el_puChargedHadronIso;
   std::vector<float> el_ecalIso;
   std::vector<float> el_hcalIso;
   std::vector<float> el_particleIso;
   std::vector<float> el_photonIso;
   std::vector<float> el_trackIso;

   std::vector<float> el_pfIso_sumChargedHadronPt;
   std::vector<float> el_pfIso_sumNeutralHadronEt;
   std::vector<float> el_pfIso_sumPhotonEt;
   std::vector<float> el_pfIso_sumPUPt;

   std::vector<int> el_isLoose;
   std::vector<int> el_isTight;
   std::vector<int> el_isRobustLoose;
   std::vector<int> el_isRobustTight;
   std::vector<int> el_isRobustHighEnergy;

   std::vector<float> el_vx;
   std::vector<float> el_vy;
   std::vector<float> el_vz;

   std::vector<bool> el_isGsf;
   std::vector<float> el_dxy;
   std::vector<float> el_dz;
   std::vector<float> el_dxyError;
   std::vector<float> el_dzError;

   std::vector<int> el_numberOfHits;

   std::vector<float> el_sigmaIetaIeta;
   std::vector<float> el_sigmaIphiIphi;
   std::vector<float> el_hadronicOverEm;
   std::vector<float> el_dr03TkSumPt;
   std::vector<float> el_dr03EcalRecHitSumEt;
   std::vector<float> el_dr03HcalTowerSumEt;
   std::vector<int> el_numberOfLostHits;

   std::vector<float> el_fbrem;
   std::vector<float> el_kf_normalizedChi2;
   std::vector<int> el_trackerLayersWithMeasurement;
   std::vector<float> el_gsf_normalizedChi2;
   std::vector<float> el_deltaEtaSuperClusterTrackAtVtx;
   std::vector<float> el_deltaPhiSuperClusterTrackAtVtx;
   std::vector<float> el_deltaEtaSeedClusterTrackAtCalo;
   std::vector<float> el_see;
   std::vector<float> el_spp;
   std::vector<float> el_superClusterEtaWidth;
   std::vector<float> el_superClusterPhiWidth;
   std::vector<float> el_full5x5_OneMinusE1x5E5x5;
   std::vector<float> el_OneMinusE1x5E5x5;
   std::vector<float> el_full5x5_r9;
   std::vector<float> el_r9;
   std::vector<float> el_eSuperClusterOverP;
   std::vector<float> el_IoEmIoP;
   std::vector<float> el_ooEmooP;
   std::vector<float> el_eleEoPout;
   std::vector<float> el_PreShowerOverRaw;
   std::vector<float> el_ecalEnergy;

   std::vector<float> el_mvaNonTrigV0;
   std::vector<float> el_mvaNonTrigCat;
   std::vector<bool> el_mvaPassMedium;
   std::vector<bool> el_mvaPassTight;

   std::vector<float> el_lepMVA;

   std::vector<float> el_lepMVA_neuRelIso;
   std::vector<float> el_lepMVA_chRelIso;
   std::vector<float> el_lepMVA_jetDR;
   std::vector<float> el_lepMVA_jetPtRatio;
   std::vector<float> el_lepMVA_jetBTagCSV;
   std::vector<float> el_lepMVA_sip3d;
   std::vector<float> el_lepMVA_dxy;
   std::vector<float> el_lepMVA_dz;
   std::vector<float> el_lepMVA_mvaId;

   std::vector<int> el_hasMCMatch;
   std::vector<float> el_gen_pt;
   std::vector<float> el_gen_eta;
   std::vector<float> el_gen_phi;
   std::vector<float> el_gen_m;
   std::vector<float> el_gen_E;
   std::vector<int> el_gen_status;
   std::vector<int> el_gen_id;
   std::vector<int> el_gen_charge;
   std::vector<float> el_gen_dr;

   std::vector<bool> el_hasMatchedConversion;
   std::vector<int> el_expectedMissingInnerHits;

   // Muons

   int mu_n;
   std::vector<float> mu_pt;
   std::vector<float> mu_eta;
   std::vector<float> mu_phi;
   std::vector<float> mu_m;
   std::vector<float> mu_E;
   std::vector<int> mu_id;
   std::vector<int> mu_charge;

   std::vector<float> mu_dB3D;
   std::vector<float> mu_edB3D;

   std::vector<float> mu_dB;
   std::vector<float> mu_edB;
   
   std::vector<float> mu_muonBest_dz;
   
   std::vector<float> mu_neutralHadronIso;
   std::vector<float> mu_chargedHadronIso;
   std::vector<float> mu_puChargedHadronIso;
   std::vector<float> mu_ecalIso;
   std::vector<float> mu_hcalIso;
   std::vector<float> mu_photonIso;
   std::vector<float> mu_trackIso;

   std::vector<float> mu_pfIso03_sumChargedHadronPt;
   std::vector<float> mu_pfIso03_sumNeutralHadronEt;
   std::vector<float> mu_pfIso03_sumPhotonEt;
   std::vector<float> mu_pfIso03_sumPUPt;

   std::vector<float> mu_pfIso04_sumChargedHadronPt;
   std::vector<float> mu_pfIso04_sumNeutralHadronEt;
   std::vector<float> mu_pfIso04_sumPhotonEt;
   std::vector<float> mu_pfIso04_sumPUPt;
   
   std::vector<float> mu_miniIso;

   std::vector<int> mu_isGlobalMuon;
   std::vector<int> mu_isTrackerMuon;
   std::vector<int> mu_isStandAloneMuon;
   std::vector<int> mu_isCaloMuon;
   std::vector<int> mu_isPFMuon;

   std::vector<float> mu_vx;
   std::vector<float> mu_vy;
   std::vector<float> mu_vz;
   
   std::vector<float> mu_segmentCompatibility;

   std::vector<bool> mu_isTightMuon;

   std::vector<int> mu_hasTrack;
   std::vector<int> mu_track_trackerLayersWithMeasurement;
   
   std::vector<int> mu_hasGlobalTrack;
   std::vector<float> mu_globalTrack_dxy;
   std::vector<float> mu_globalTrack_dz;
   std::vector<float> mu_globalTrack_dxyError;
   std::vector<float> mu_globalTrack_dzError;   
   std::vector<float> mu_globalTrack_normalizedChi2;
   
   std::vector<float> mu_combinedQuality_chi2LocalPosition;
   std::vector<float> mu_combinedQuality_trkKink;

   std::vector<int> mu_hasInnerTrack;
   std::vector<float> mu_innerTrack_dxy;
   std::vector<float> mu_innerTrack_dz;
   std::vector<float> mu_innerTrack_dxyError;
   std::vector<float> mu_innerTrack_dzError;
   std::vector<float> mu_innerTrack_normalizedChi2;
   
   std::vector<float> mu_innerTrack_validFraction;

   std::vector<int> mu_bestTrackType;
   std::vector<float> mu_bestTrack_dxy;
   std::vector<float> mu_bestTrack_dz;
   std::vector<float> mu_bestTrack_dxyError;
   std::vector<float> mu_bestTrack_dzError;
   std::vector<float> mu_bestTrack_normalizedChi2;

   std::vector<float> mu_innerTrack_pt;
   std::vector<float> mu_innerTrack_ptError;
   std::vector<int> mu_innerTrack_numberOfValidPixelHits;

   std::vector<int> mu_numberOfMatches;
   std::vector<int> mu_numberOfMatchedStations;
   std::vector<int> mu_numberOfValidMuonHits;

   std::vector<float> mu_lepMVA;

   std::vector<float> mu_lepMVA_neuRelIso;
   std::vector<float> mu_lepMVA_chRelIso;
   std::vector<float> mu_lepMVA_jetDR;
   std::vector<float> mu_lepMVA_jetPtRatio;
   std::vector<float> mu_lepMVA_jetBTagCSV;
   std::vector<float> mu_lepMVA_sip3d;
   std::vector<float> mu_lepMVA_dxy;
   std::vector<float> mu_lepMVA_dz;
   std::vector<float> mu_lepMVA_mvaId;

   std::vector<int> mu_hasMCMatch;
   std::vector<float> mu_gen_pt;
   std::vector<float> mu_gen_eta;
   std::vector<float> mu_gen_phi;
   std::vector<float> mu_gen_m;
   std::vector<int> mu_gen_status;
   std::vector<int> mu_gen_id;
   std::vector<int> mu_gen_charge;
   std::vector<float> mu_gen_dr;

   // Taus
   
   int tau_n;
   std::vector<float> tau_pt;
   std::vector<float> tau_eta;
   std::vector<float> tau_phi;
   std::vector<float> tau_m;
   std::vector<float> tau_E;
   std::vector<int> tau_id;
   std::vector<int> tau_charge;
   
   std::vector<bool> tau_hasLeadChargedHadrCand;
   std::vector<float> tau_leadingTrackPt;
   std::vector<int> tau_leadingTrackCharge;
   
   std::vector<int> tau_decayMode;
//   std::vector<float> tau_decayModeFindingOldDMs;
   std::vector<float> tau_decayModeFindingNewDMs;
   
   std::vector<float> tau_puCorrPtSum;
   std::vector<float> tau_neutralIsoPtSum;
   std::vector<float> tau_chargedIsoPtSum;
   std::vector<float> tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;
   
   std::vector<float> tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   std::vector<float> tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   std::vector<float> tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
   std::vector<float> tau_byMediumIsolationMVA3newDMwLT;

   std::vector<float> tau_againstMuonLoose3;
   std::vector<float> tau_againstMuonTight3;

   std::vector<float> tau_againstElectronVLooseMVA5;
   std::vector<float> tau_againstElectronLooseMVA5;
   std::vector<float> tau_againstElectronMediumMVA5;
   
   std::vector<float> tau_pfEssential_jet_pt;
   std::vector<float> tau_pfEssential_jet_eta;
   std::vector<float> tau_pfEssential_jet_phi;
   std::vector<float> tau_pfEssential_jet_m;

   std::vector<float> tau_pfEssential_jetCorr_pt;
   std::vector<float> tau_pfEssential_jetCorr_eta;
   std::vector<float> tau_pfEssential_jetCorr_phi;
   std::vector<float> tau_pfEssential_jetCorr_m;
   
   std::vector<bool> tau_pfEssential_hasSV;
   std::vector<float> tau_pfEssential_sv_x;
   std::vector<float> tau_pfEssential_sv_y;
   std::vector<float> tau_pfEssential_sv_z;
   
   std::vector<float> tau_pfEssential_flightLengthSig;
   std::vector<float> tau_pfEssential_dxy;
   std::vector<float> tau_pfEssential_dxy_error;
   std::vector<float> tau_pfEssential_dxy_Sig;
   
   // Jets

   int jet_n;
   std::vector<float> jet_pt;
   std::vector<float> jet_eta;
   std::vector<float> jet_phi;
   std::vector<float> jet_m;
   std::vector<float> jet_E;

   std::vector<int> jet_ntrk;

   std::vector<float> jet_JBP;
   std::vector<float> jet_JP;
   std::vector<float> jet_TCHP;
   std::vector<float> jet_TCHE;
   std::vector<float> jet_SSVHE;
   std::vector<float> jet_SSVHP;
   std::vector<float> jet_CMVA;
   std::vector<float> jet_CSV;
   std::vector<float> jet_CSVv2;
   std::vector<int> jet_partonFlavour;
   std::vector<int> jet_hadronFlavour;

   std::vector<float> jet_neutralHadronEnergy;
   std::vector<float> jet_neutralEmEnergy;
   std::vector<float> jet_chargedHadronEnergy;
   std::vector<float> jet_chargedEmEnergy;
   std::vector<float> jet_electronEnergy;
   std::vector<float> jet_muonEnergy;
   std::vector<float> jet_photonEnergy;

   std::vector<int> jet_chargedMultiplicity;
   std::vector<int> jet_neutralMultiplicity;
   std::vector<int> jet_chargedHadronMultiplicity;
   
   std::vector<float> jet_jecFactorUncorrected;
   std::vector<float> jet_jecFactorL1FastJet;
   std::vector<float> jet_jecFactorL2Relative;
   std::vector<float> jet_jecFactorL3Absolute;
   
   std::vector<float> jet_Unc;
   
   std::vector<float> jet_pileupJetId;
   std::vector<float> jet_looseJetID;
   std::vector<float> jet_tightJetID;

   std::vector<bool> jet_hasGenJet;   
   std::vector<float> jet_genJet_pt;
   std::vector<float> jet_genJet_eta;
   std::vector<float> jet_genJet_phi;
   std::vector<float> jet_genJet_m;
   std::vector<float> jet_genJet_E;
   std::vector<int> jet_genJet_status;
   std::vector<int> jet_genJet_id;

   std::vector<bool> jet_hasGenParton;
   std::vector<float> jet_genParton_pt;
   std::vector<float> jet_genParton_eta;
   std::vector<float> jet_genParton_phi;
   std::vector<float> jet_genParton_m;
   std::vector<float> jet_genParton_E;
   std::vector<int> jet_genParton_status;
   std::vector<int> jet_genParton_id;
   
   // Puppi Jets

   int jetPuppi_n;
   std::vector<float> jetPuppi_pt;
   std::vector<float> jetPuppi_eta;
   std::vector<float> jetPuppi_phi;
   std::vector<float> jetPuppi_m;
   std::vector<float> jetPuppi_E;

   std::vector<int> jetPuppi_ntrk;

   std::vector<float> jetPuppi_JBP;
   std::vector<float> jetPuppi_JP;
   std::vector<float> jetPuppi_TCHP;
   std::vector<float> jetPuppi_TCHE;
   std::vector<float> jetPuppi_SSVHE;
   std::vector<float> jetPuppi_SSVHP;
   std::vector<float> jetPuppi_CMVA;
   std::vector<float> jetPuppi_CSV;
   std::vector<float> jetPuppi_CSVv2;
   std::vector<int> jetPuppi_partonFlavour;
   std::vector<int> jetPuppi_hadronFlavour;

   std::vector<float> jetPuppi_neutralHadronEnergy;
   std::vector<float> jetPuppi_neutralEmEnergy;
   std::vector<float> jetPuppi_chargedHadronEnergy;
   std::vector<float> jetPuppi_chargedEmEnergy;
   std::vector<float> jetPuppi_electronEnergy;
   std::vector<float> jetPuppi_muonEnergy;
   std::vector<float> jetPuppi_photonEnergy;

   std::vector<int> jetPuppi_chargedMultiplicity;
   std::vector<int> jetPuppi_neutralMultiplicity;
   std::vector<int> jetPuppi_chargedHadronMultiplicity;
   
   std::vector<float> jetPuppi_jecFactorUncorrected;
   std::vector<float> jetPuppi_jecFactorL1FastJet;
   std::vector<float> jetPuppi_jecFactorL2Relative;
   std::vector<float> jetPuppi_jecFactorL3Absolute;
   
   std::vector<float> jetPuppi_pileupJetId;

   std::vector<bool> jetPuppi_hasGenJet;   
   std::vector<float> jetPuppi_genJet_pt;
   std::vector<float> jetPuppi_genJet_eta;
   std::vector<float> jetPuppi_genJet_phi;
   std::vector<float> jetPuppi_genJet_m;
   std::vector<float> jetPuppi_genJet_E;
   std::vector<int> jetPuppi_genJet_status;
   std::vector<int> jetPuppi_genJet_id;

   std::vector<bool> jetPuppi_hasGenParton;   
   std::vector<float> jetPuppi_genParton_pt;
   std::vector<float> jetPuppi_genParton_eta;
   std::vector<float> jetPuppi_genParton_phi;
   std::vector<float> jetPuppi_genParton_m;
   std::vector<float> jetPuppi_genParton_E;
   std::vector<int> jetPuppi_genParton_status;
   std::vector<int> jetPuppi_genParton_id;
   
   // GenJets

   int genJet_n;
   std::vector<float> genJet_pt;
   std::vector<float> genJet_eta;
   std::vector<float> genJet_phi;
   std::vector<float> genJet_m;
   std::vector<float> genJet_E;
   std::vector<float> genJet_emEnergy;
   std::vector<float> genJet_hadEnergy;
   std::vector<float> genJet_invisibleEnergy;
   std::vector<float> genJet_auxiliaryEnergy;
   std::vector<int>   genJet_flavour;
  
  
   //PFcands
   
   int pfcand_n;
   std::vector<float> pfcand_pt;
   std::vector<float> pfcand_eta;
   std::vector<float> pfcand_phi;
   std::vector<float> pfcand_E;
   std::vector<float> pfcand_charge;
   std::vector<int> pfcand_id;
   std::vector<float> pfcand_dz;


   // ttH
   int mc_truth_tth_channel;

   int mc_truth_thq_channel;
   
   TLorentzVector mc_truth_h0_p4;

   TLorentzVector mc_truth_h0W1_p4;
   TLorentzVector mc_truth_h0W2_p4;
   TLorentzVector mc_truth_h0Wl1_p4;
   TLorentzVector mc_truth_h0Wnu1_p4;
   TLorentzVector mc_truth_h0Wtau1_p4;
   TLorentzVector mc_truth_h0Wnutau1_p4;
   TLorentzVector mc_truth_h0Wtaul1_p4;
   TLorentzVector mc_truth_h0Wtaunu1_p4;
   TLorentzVector mc_truth_h0Wtaunutau1_p4;
   TLorentzVector mc_truth_h0Wl2_p4;
   TLorentzVector mc_truth_h0Wnu2_p4;
   TLorentzVector mc_truth_h0Wtau2_p4;
   TLorentzVector mc_truth_h0Wnutau2_p4;
   TLorentzVector mc_truth_h0Wtaul2_p4;
   TLorentzVector mc_truth_h0Wtaunu2_p4;
   TLorentzVector mc_truth_h0Wtaunutau2_p4;
   TLorentzVector mc_truth_h0Wq11_p4;
   TLorentzVector mc_truth_h0Wq21_p4;
   TLorentzVector mc_truth_h0Wq12_p4;
   TLorentzVector mc_truth_h0Wq22_p4;

   TLorentzVector mc_truth_h0Z1_p4;
   TLorentzVector mc_truth_h0Z2_p4;
   TLorentzVector mc_truth_h0Zl11_p4;
   TLorentzVector mc_truth_h0Zl21_p4;
   TLorentzVector mc_truth_h0Ztau11_p4;
   TLorentzVector mc_truth_h0Ztau21_p4;
   TLorentzVector mc_truth_h0Ztaul11_p4;
   TLorentzVector mc_truth_h0Ztaul21_p4;
   TLorentzVector mc_truth_h0Ztaunu11_p4;
   TLorentzVector mc_truth_h0Ztaunu21_p4;
   TLorentzVector mc_truth_h0Ztaunutau11_p4;
   TLorentzVector mc_truth_h0Ztaunutau21_p4;
   TLorentzVector mc_truth_h0Zq11_p4;
   TLorentzVector mc_truth_h0Zq21_p4;
   TLorentzVector mc_truth_h0Zl12_p4;
   TLorentzVector mc_truth_h0Zl22_p4;
   TLorentzVector mc_truth_h0Ztau12_p4;
   TLorentzVector mc_truth_h0Ztau22_p4;
   TLorentzVector mc_truth_h0Ztaul12_p4;
   TLorentzVector mc_truth_h0Ztaul22_p4;
   TLorentzVector mc_truth_h0Ztaunu12_p4;
   TLorentzVector mc_truth_h0Ztaunu22_p4;
   TLorentzVector mc_truth_h0Ztaunutau12_p4;
   TLorentzVector mc_truth_h0Ztaunutau22_p4;
   TLorentzVector mc_truth_h0Zq12_p4;
   TLorentzVector mc_truth_h0Zq22_p4;
   TLorentzVector mc_truth_h0Znu11_p4;
   TLorentzVector mc_truth_h0Znu21_p4;
   TLorentzVector mc_truth_h0Znu12_p4;
   TLorentzVector mc_truth_h0Znu22_p4;

   TLorentzVector mc_truth_h0tau1_p4;
   TLorentzVector mc_truth_h0tau2_p4;
   TLorentzVector mc_truth_h0taul1_p4;
   TLorentzVector mc_truth_h0taunutau1_p4;
   TLorentzVector mc_truth_h0taunu1_p4;
   TLorentzVector mc_truth_h0taul2_p4;
   TLorentzVector mc_truth_h0taunutau2_p4;
   TLorentzVector mc_truth_h0taunu2_p4;

   TLorentzVector mc_truth_h0b1_p4;
   TLorentzVector mc_truth_h0b2_p4;
   TLorentzVector mc_truth_h0b1_IS_p4;
   TLorentzVector mc_truth_h0b2_IS_p4;
   
   TLorentzVector mc_truth_t1_p4;
   TLorentzVector mc_truth_t2_p4;
   TLorentzVector mc_truth_tb1_p4;
   TLorentzVector mc_truth_tb2_p4;

   TLorentzVector mc_truth_tW1_p4;
   TLorentzVector mc_truth_tWnu1_p4;
   TLorentzVector mc_truth_tWnutau1_p4;
   TLorentzVector mc_truth_tWl1_p4;
   TLorentzVector mc_truth_tWtau1_p4;
   TLorentzVector mc_truth_tWtaunu1_p4;
   TLorentzVector mc_truth_tWtaunutau1_p4;
   TLorentzVector mc_truth_tWtaul1_p4;
   TLorentzVector mc_truth_tWq11_p4;
   TLorentzVector mc_truth_tWq21_p4;

   TLorentzVector mc_truth_tW2_p4;
   TLorentzVector mc_truth_tWnu2_p4;
   TLorentzVector mc_truth_tWnutau2_p4;
   TLorentzVector mc_truth_tWl2_p4;
   TLorentzVector mc_truth_tWtau2_p4;
   TLorentzVector mc_truth_tWtaunu2_p4;
   TLorentzVector mc_truth_tWtaunutau2_p4;
   TLorentzVector mc_truth_tWtaul2_p4;
   TLorentzVector mc_truth_tWq12_p4;
   TLorentzVector mc_truth_tWq22_p4;

   TLorentzVector mc_truth_j1_p4;
   TLorentzVector mc_truth_j2_p4;
   TLorentzVector mc_truth_j3_p4;

   // pt

   float mc_truth_h0_pt;

   float mc_truth_h0W1_pt;
   float mc_truth_h0W2_pt;
   float mc_truth_h0Wl1_pt;
   float mc_truth_h0Wnu1_pt;
   float mc_truth_h0Wtau1_pt;
   float mc_truth_h0Wnutau1_pt;
   float mc_truth_h0Wtaul1_pt;
   float mc_truth_h0Wtaunu1_pt;
   float mc_truth_h0Wtaunutau1_pt;
   float mc_truth_h0Wl2_pt;
   float mc_truth_h0Wnu2_pt;
   float mc_truth_h0Wtau2_pt;
   float mc_truth_h0Wnutau2_pt;
   float mc_truth_h0Wtaul2_pt;
   float mc_truth_h0Wtaunu2_pt;
   float mc_truth_h0Wtaunutau2_pt;
   float mc_truth_h0Wq11_pt;
   float mc_truth_h0Wq21_pt;
   float mc_truth_h0Wq12_pt;
   float mc_truth_h0Wq22_pt;

   float mc_truth_h0Z1_pt;
   float mc_truth_h0Z2_pt;
   float mc_truth_h0Zl11_pt;
   float mc_truth_h0Zl21_pt;
   float mc_truth_h0Ztau11_pt;
   float mc_truth_h0Ztau21_pt;
   float mc_truth_h0Ztaul11_pt;
   float mc_truth_h0Ztaul21_pt;
   float mc_truth_h0Ztaunu11_pt;
   float mc_truth_h0Ztaunu21_pt;
   float mc_truth_h0Ztaunutau11_pt;
   float mc_truth_h0Ztaunutau21_pt;
   float mc_truth_h0Zq11_pt;
   float mc_truth_h0Zq21_pt;
   float mc_truth_h0Zl12_pt;
   float mc_truth_h0Zl22_pt;
   float mc_truth_h0Ztau12_pt;
   float mc_truth_h0Ztau22_pt;
   float mc_truth_h0Ztaul12_pt;
   float mc_truth_h0Ztaul22_pt;
   float mc_truth_h0Ztaunu12_pt;
   float mc_truth_h0Ztaunu22_pt;
   float mc_truth_h0Ztaunutau12_pt;
   float mc_truth_h0Ztaunutau22_pt;
   float mc_truth_h0Zq12_pt;
   float mc_truth_h0Zq22_pt;
   float mc_truth_h0Znu11_pt;
   float mc_truth_h0Znu21_pt;
   float mc_truth_h0Znu12_pt;
   float mc_truth_h0Znu22_pt;

   float mc_truth_h0tau1_pt;
   float mc_truth_h0tau2_pt;
   float mc_truth_h0taul1_pt;
   float mc_truth_h0taunutau1_pt;
   float mc_truth_h0taunu1_pt;
   float mc_truth_h0taul2_pt;
   float mc_truth_h0taunutau2_pt;
   float mc_truth_h0taunu2_pt;

   float mc_truth_h0b1_pt;
   float mc_truth_h0b2_pt;
   float mc_truth_h0b1_IS_pt;
   float mc_truth_h0b2_IS_pt;
   
   float mc_truth_t1_pt;
   float mc_truth_t2_pt;
   float mc_truth_tb1_pt;
   float mc_truth_tb2_pt;

   float mc_truth_tW1_pt;
   float mc_truth_tWnu1_pt;
   float mc_truth_tWnutau1_pt;
   float mc_truth_tWl1_pt;
   float mc_truth_tWtau1_pt;
   float mc_truth_tWtaunu1_pt;
   float mc_truth_tWtaunutau1_pt;
   float mc_truth_tWtaul1_pt;
   float mc_truth_tWq11_pt;
   float mc_truth_tWq21_pt;

   float mc_truth_tW2_pt;
   float mc_truth_tWnu2_pt;
   float mc_truth_tWnutau2_pt;
   float mc_truth_tWl2_pt;
   float mc_truth_tWtau2_pt;
   float mc_truth_tWtaunu2_pt;
   float mc_truth_tWtaunutau2_pt;
   float mc_truth_tWtaul2_pt;
   float mc_truth_tWq12_pt;
   float mc_truth_tWq22_pt;

   float mc_truth_j1_pt;
   float mc_truth_j2_pt;
   float mc_truth_j3_pt;

   // eta

   float mc_truth_h0_eta;

   float mc_truth_h0W1_eta;
   float mc_truth_h0W2_eta;
   float mc_truth_h0Wl1_eta;
   float mc_truth_h0Wnu1_eta;
   float mc_truth_h0Wtau1_eta;
   float mc_truth_h0Wnutau1_eta;
   float mc_truth_h0Wtaul1_eta;
   float mc_truth_h0Wtaunu1_eta;
   float mc_truth_h0Wtaunutau1_eta;
   float mc_truth_h0Wl2_eta;
   float mc_truth_h0Wnu2_eta;
   float mc_truth_h0Wtau2_eta;
   float mc_truth_h0Wnutau2_eta;
   float mc_truth_h0Wtaul2_eta;
   float mc_truth_h0Wtaunu2_eta;
   float mc_truth_h0Wtaunutau2_eta;
   float mc_truth_h0Wq11_eta;
   float mc_truth_h0Wq21_eta;
   float mc_truth_h0Wq12_eta;
   float mc_truth_h0Wq22_eta;

   float mc_truth_h0Z1_eta;
   float mc_truth_h0Z2_eta;
   float mc_truth_h0Zl11_eta;
   float mc_truth_h0Zl21_eta;
   float mc_truth_h0Ztau11_eta;
   float mc_truth_h0Ztau21_eta;
   float mc_truth_h0Ztaul11_eta;
   float mc_truth_h0Ztaul21_eta;
   float mc_truth_h0Ztaunu11_eta;
   float mc_truth_h0Ztaunu21_eta;
   float mc_truth_h0Ztaunutau11_eta;
   float mc_truth_h0Ztaunutau21_eta;
   float mc_truth_h0Zq11_eta;
   float mc_truth_h0Zq21_eta;
   float mc_truth_h0Zl12_eta;
   float mc_truth_h0Zl22_eta;
   float mc_truth_h0Ztau12_eta;
   float mc_truth_h0Ztau22_eta;
   float mc_truth_h0Ztaul12_eta;
   float mc_truth_h0Ztaul22_eta;
   float mc_truth_h0Ztaunu12_eta;
   float mc_truth_h0Ztaunu22_eta;
   float mc_truth_h0Ztaunutau12_eta;
   float mc_truth_h0Ztaunutau22_eta;
   float mc_truth_h0Zq12_eta;
   float mc_truth_h0Zq22_eta;
   float mc_truth_h0Znu11_eta;
   float mc_truth_h0Znu21_eta;
   float mc_truth_h0Znu12_eta;
   float mc_truth_h0Znu22_eta;

   float mc_truth_h0tau1_eta;
   float mc_truth_h0tau2_eta;
   float mc_truth_h0taul1_eta;
   float mc_truth_h0taunutau1_eta;
   float mc_truth_h0taunu1_eta;
   float mc_truth_h0taul2_eta;
   float mc_truth_h0taunutau2_eta;
   float mc_truth_h0taunu2_eta;

   float mc_truth_h0b1_eta;
   float mc_truth_h0b2_eta;
   float mc_truth_h0b1_IS_eta;
   float mc_truth_h0b2_IS_eta;
   
   float mc_truth_t1_eta;
   float mc_truth_t2_eta;
   float mc_truth_tb1_eta;
   float mc_truth_tb2_eta;

   float mc_truth_tW1_eta;
   float mc_truth_tWnu1_eta;
   float mc_truth_tWnutau1_eta;
   float mc_truth_tWl1_eta;
   float mc_truth_tWtau1_eta;
   float mc_truth_tWtaunu1_eta;
   float mc_truth_tWtaunutau1_eta;
   float mc_truth_tWtaul1_eta;
   float mc_truth_tWq11_eta;
   float mc_truth_tWq21_eta;

   float mc_truth_tW2_eta;
   float mc_truth_tWnu2_eta;
   float mc_truth_tWnutau2_eta;
   float mc_truth_tWl2_eta;
   float mc_truth_tWtau2_eta;
   float mc_truth_tWtaunu2_eta;
   float mc_truth_tWtaunutau2_eta;
   float mc_truth_tWtaul2_eta;
   float mc_truth_tWq12_eta;
   float mc_truth_tWq22_eta;

   float mc_truth_j1_eta;
   float mc_truth_j2_eta;
   float mc_truth_j3_eta;
   
   // phi

   float mc_truth_h0_phi;

   float mc_truth_h0W1_phi;
   float mc_truth_h0W2_phi;
   float mc_truth_h0Wl1_phi;
   float mc_truth_h0Wnu1_phi;
   float mc_truth_h0Wtau1_phi;
   float mc_truth_h0Wnutau1_phi;
   float mc_truth_h0Wtaul1_phi;
   float mc_truth_h0Wtaunu1_phi;
   float mc_truth_h0Wtaunutau1_phi;
   float mc_truth_h0Wl2_phi;
   float mc_truth_h0Wnu2_phi;
   float mc_truth_h0Wtau2_phi;
   float mc_truth_h0Wnutau2_phi;
   float mc_truth_h0Wtaul2_phi;
   float mc_truth_h0Wtaunu2_phi;
   float mc_truth_h0Wtaunutau2_phi;
   float mc_truth_h0Wq11_phi;
   float mc_truth_h0Wq21_phi;
   float mc_truth_h0Wq12_phi;
   float mc_truth_h0Wq22_phi;

   float mc_truth_h0Z1_phi;
   float mc_truth_h0Z2_phi;
   float mc_truth_h0Zl11_phi;
   float mc_truth_h0Zl21_phi;
   float mc_truth_h0Ztau11_phi;
   float mc_truth_h0Ztau21_phi;
   float mc_truth_h0Ztaul11_phi;
   float mc_truth_h0Ztaul21_phi;
   float mc_truth_h0Ztaunu11_phi;
   float mc_truth_h0Ztaunu21_phi;
   float mc_truth_h0Ztaunutau11_phi;
   float mc_truth_h0Ztaunutau21_phi;
   float mc_truth_h0Zq11_phi;
   float mc_truth_h0Zq21_phi;
   float mc_truth_h0Zl12_phi;
   float mc_truth_h0Zl22_phi;
   float mc_truth_h0Ztau12_phi;
   float mc_truth_h0Ztau22_phi;
   float mc_truth_h0Ztaul12_phi;
   float mc_truth_h0Ztaul22_phi;
   float mc_truth_h0Ztaunu12_phi;
   float mc_truth_h0Ztaunu22_phi;
   float mc_truth_h0Ztaunutau12_phi;
   float mc_truth_h0Ztaunutau22_phi;
   float mc_truth_h0Zq12_phi;
   float mc_truth_h0Zq22_phi;
   float mc_truth_h0Znu11_phi;
   float mc_truth_h0Znu21_phi;
   float mc_truth_h0Znu12_phi;
   float mc_truth_h0Znu22_phi;

   float mc_truth_h0tau1_phi;
   float mc_truth_h0tau2_phi;
   float mc_truth_h0taul1_phi;
   float mc_truth_h0taunutau1_phi;
   float mc_truth_h0taunu1_phi;
   float mc_truth_h0taul2_phi;
   float mc_truth_h0taunutau2_phi;
   float mc_truth_h0taunu2_phi;

   float mc_truth_h0b1_phi;
   float mc_truth_h0b2_phi;
   float mc_truth_h0b1_IS_phi;
   float mc_truth_h0b2_IS_phi;
   
   float mc_truth_t1_phi;
   float mc_truth_t2_phi;
   float mc_truth_tb1_phi;
   float mc_truth_tb2_phi;

   float mc_truth_tW1_phi;
   float mc_truth_tWnu1_phi;
   float mc_truth_tWnutau1_phi;
   float mc_truth_tWl1_phi;
   float mc_truth_tWtau1_phi;
   float mc_truth_tWtaunu1_phi;
   float mc_truth_tWtaunutau1_phi;
   float mc_truth_tWtaul1_phi;
   float mc_truth_tWq11_phi;
   float mc_truth_tWq21_phi;

   float mc_truth_tW2_phi;
   float mc_truth_tWnu2_phi;
   float mc_truth_tWnutau2_phi;
   float mc_truth_tWl2_phi;
   float mc_truth_tWtau2_phi;
   float mc_truth_tWtaunu2_phi;
   float mc_truth_tWtaunutau2_phi;
   float mc_truth_tWtaul2_phi;
   float mc_truth_tWq12_phi;
   float mc_truth_tWq22_phi;

   float mc_truth_j1_phi;
   float mc_truth_j2_phi;
   float mc_truth_j3_phi;
   
   // E

   float mc_truth_h0_E;

   float mc_truth_h0W1_E;
   float mc_truth_h0W2_E;
   float mc_truth_h0Wl1_E;
   float mc_truth_h0Wnu1_E;
   float mc_truth_h0Wtau1_E;
   float mc_truth_h0Wnutau1_E;
   float mc_truth_h0Wtaul1_E;
   float mc_truth_h0Wtaunu1_E;
   float mc_truth_h0Wtaunutau1_E;
   float mc_truth_h0Wl2_E;
   float mc_truth_h0Wnu2_E;
   float mc_truth_h0Wtau2_E;
   float mc_truth_h0Wnutau2_E;
   float mc_truth_h0Wtaul2_E;
   float mc_truth_h0Wtaunu2_E;
   float mc_truth_h0Wtaunutau2_E;
   float mc_truth_h0Wq11_E;
   float mc_truth_h0Wq21_E;
   float mc_truth_h0Wq12_E;
   float mc_truth_h0Wq22_E;

   float mc_truth_h0Z1_E;
   float mc_truth_h0Z2_E;
   float mc_truth_h0Zl11_E;
   float mc_truth_h0Zl21_E;
   float mc_truth_h0Ztau11_E;
   float mc_truth_h0Ztau21_E;
   float mc_truth_h0Ztaul11_E;
   float mc_truth_h0Ztaul21_E;
   float mc_truth_h0Ztaunu11_E;
   float mc_truth_h0Ztaunu21_E;
   float mc_truth_h0Ztaunutau11_E;
   float mc_truth_h0Ztaunutau21_E;
   float mc_truth_h0Zq11_E;
   float mc_truth_h0Zq21_E;
   float mc_truth_h0Zl12_E;
   float mc_truth_h0Zl22_E;
   float mc_truth_h0Ztau12_E;
   float mc_truth_h0Ztau22_E;
   float mc_truth_h0Ztaul12_E;
   float mc_truth_h0Ztaul22_E;
   float mc_truth_h0Ztaunu12_E;
   float mc_truth_h0Ztaunu22_E;
   float mc_truth_h0Ztaunutau12_E;
   float mc_truth_h0Ztaunutau22_E;
   float mc_truth_h0Zq12_E;
   float mc_truth_h0Zq22_E;
   float mc_truth_h0Znu11_E;
   float mc_truth_h0Znu21_E;
   float mc_truth_h0Znu12_E;
   float mc_truth_h0Znu22_E;

   float mc_truth_h0tau1_E;
   float mc_truth_h0tau2_E;
   float mc_truth_h0taul1_E;
   float mc_truth_h0taunutau1_E;
   float mc_truth_h0taunu1_E;
   float mc_truth_h0taul2_E;
   float mc_truth_h0taunutau2_E;
   float mc_truth_h0taunu2_E;

   float mc_truth_h0b1_E;
   float mc_truth_h0b2_E;
   float mc_truth_h0b1_IS_E;
   float mc_truth_h0b2_IS_E;
   
   float mc_truth_t1_E;
   float mc_truth_t2_E;
   float mc_truth_tb1_E;
   float mc_truth_tb2_E;

   float mc_truth_tW1_E;
   float mc_truth_tWnu1_E;
   float mc_truth_tWnutau1_E;
   float mc_truth_tWl1_E;
   float mc_truth_tWtau1_E;
   float mc_truth_tWtaunu1_E;
   float mc_truth_tWtaunutau1_E;
   float mc_truth_tWtaul1_E;
   float mc_truth_tWq11_E;
   float mc_truth_tWq21_E;

   float mc_truth_tW2_E;
   float mc_truth_tWnu2_E;
   float mc_truth_tWnutau2_E;
   float mc_truth_tWl2_E;
   float mc_truth_tWtau2_E;
   float mc_truth_tWtaunu2_E;
   float mc_truth_tWtaunutau2_E;
   float mc_truth_tWtaul2_E;
   float mc_truth_tWq12_E;
   float mc_truth_tWq22_E;

   float mc_truth_j1_E;
   float mc_truth_j2_E;
   float mc_truth_j3_E;
   
   // pdgId

   int mc_truth_h0_id;

   int mc_truth_h0W1_id;
   int mc_truth_h0W2_id;
   int mc_truth_h0Wl1_id;
   int mc_truth_h0Wnu1_id;
   int mc_truth_h0Wtau1_id;
   int mc_truth_h0Wnutau1_id;
   int mc_truth_h0Wtaul1_id;
   int mc_truth_h0Wtaunu1_id;
   int mc_truth_h0Wtaunutau1_id;
   int mc_truth_h0Wl2_id;
   int mc_truth_h0Wnu2_id;
   int mc_truth_h0Wtau2_id;
   int mc_truth_h0Wnutau2_id;
   int mc_truth_h0Wtaul2_id;
   int mc_truth_h0Wtaunu2_id;
   int mc_truth_h0Wtaunutau2_id;
   int mc_truth_h0Wq11_id;
   int mc_truth_h0Wq21_id;
   int mc_truth_h0Wq12_id;
   int mc_truth_h0Wq22_id;

   int mc_truth_h0Z1_id;
   int mc_truth_h0Z2_id;
   int mc_truth_h0Zl11_id;
   int mc_truth_h0Zl21_id;
   int mc_truth_h0Ztau11_id;
   int mc_truth_h0Ztau21_id;
   int mc_truth_h0Ztaul11_id;
   int mc_truth_h0Ztaul21_id;
   int mc_truth_h0Ztaunu11_id;
   int mc_truth_h0Ztaunu21_id;
   int mc_truth_h0Ztaunutau11_id;
   int mc_truth_h0Ztaunutau21_id;
   int mc_truth_h0Zq11_id;
   int mc_truth_h0Zq21_id;
   int mc_truth_h0Zl12_id;
   int mc_truth_h0Zl22_id;
   int mc_truth_h0Ztau12_id;
   int mc_truth_h0Ztau22_id;
   int mc_truth_h0Ztaul12_id;
   int mc_truth_h0Ztaul22_id;
   int mc_truth_h0Ztaunu12_id;
   int mc_truth_h0Ztaunu22_id;
   int mc_truth_h0Ztaunutau12_id;
   int mc_truth_h0Ztaunutau22_id;
   int mc_truth_h0Zq12_id;
   int mc_truth_h0Zq22_id;
   int mc_truth_h0Znu11_id;
   int mc_truth_h0Znu21_id;
   int mc_truth_h0Znu12_id;
   int mc_truth_h0Znu22_id;

   int mc_truth_h0tau1_id;
   int mc_truth_h0tau2_id;
   int mc_truth_h0taul1_id;
   int mc_truth_h0taunutau1_id;
   int mc_truth_h0taunu1_id;
   int mc_truth_h0taul2_id;
   int mc_truth_h0taunutau2_id;
   int mc_truth_h0taunu2_id;

   int mc_truth_h0b1_id;
   int mc_truth_h0b2_id;
   int mc_truth_h0b1_IS_id;
   int mc_truth_h0b2_IS_id;
   
   int mc_truth_t1_id;
   int mc_truth_t2_id;
   int mc_truth_tb1_id;
   int mc_truth_tb2_id;

   int mc_truth_tW1_id;
   int mc_truth_tWnu1_id;
   int mc_truth_tWnutau1_id;
   int mc_truth_tWl1_id;
   int mc_truth_tWtau1_id;
   int mc_truth_tWtaunu1_id;
   int mc_truth_tWtaunutau1_id;
   int mc_truth_tWtaul1_id;
   int mc_truth_tWq11_id;
   int mc_truth_tWq21_id;

   int mc_truth_tW2_id;
   int mc_truth_tWnu2_id;
   int mc_truth_tWnutau2_id;
   int mc_truth_tWl2_id;
   int mc_truth_tWtau2_id;
   int mc_truth_tWtaunu2_id;
   int mc_truth_tWtaunutau2_id;
   int mc_truth_tWtaul2_id;
   int mc_truth_tWq12_id;
   int mc_truth_tWq22_id;

   int mc_truth_j1_id;
   int mc_truth_j2_id;
   int mc_truth_j3_id;

   // status

   int mc_truth_h0_status;

   int mc_truth_h0W1_status;
   int mc_truth_h0W2_status;
   int mc_truth_h0Wl1_status;
   int mc_truth_h0Wnu1_status;
   int mc_truth_h0Wtau1_status;
   int mc_truth_h0Wnutau1_status;
   int mc_truth_h0Wtaul1_status;
   int mc_truth_h0Wtaunu1_status;
   int mc_truth_h0Wtaunutau1_status;
   int mc_truth_h0Wl2_status;
   int mc_truth_h0Wnu2_status;
   int mc_truth_h0Wtau2_status;
   int mc_truth_h0Wnutau2_status;
   int mc_truth_h0Wtaul2_status;
   int mc_truth_h0Wtaunu2_status;
   int mc_truth_h0Wtaunutau2_status;
   int mc_truth_h0Wq11_status;
   int mc_truth_h0Wq21_status;
   int mc_truth_h0Wq12_status;
   int mc_truth_h0Wq22_status;

   int mc_truth_h0Z1_status;
   int mc_truth_h0Z2_status;
   int mc_truth_h0Zl11_status;
   int mc_truth_h0Zl21_status;
   int mc_truth_h0Ztau11_status;
   int mc_truth_h0Ztau21_status;
   int mc_truth_h0Ztaul11_status;
   int mc_truth_h0Ztaul21_status;
   int mc_truth_h0Ztaunu11_status;
   int mc_truth_h0Ztaunu21_status;
   int mc_truth_h0Ztaunutau11_status;
   int mc_truth_h0Ztaunutau21_status;
   int mc_truth_h0Zq11_status;
   int mc_truth_h0Zq21_status;
   int mc_truth_h0Zl12_status;
   int mc_truth_h0Zl22_status;
   int mc_truth_h0Ztau12_status;
   int mc_truth_h0Ztau22_status;
   int mc_truth_h0Ztaul12_status;
   int mc_truth_h0Ztaul22_status;
   int mc_truth_h0Ztaunu12_status;
   int mc_truth_h0Ztaunu22_status;
   int mc_truth_h0Ztaunutau12_status;
   int mc_truth_h0Ztaunutau22_status;
   int mc_truth_h0Zq12_status;
   int mc_truth_h0Zq22_status;
   int mc_truth_h0Znu11_status;
   int mc_truth_h0Znu21_status;
   int mc_truth_h0Znu12_status;
   int mc_truth_h0Znu22_status;

   int mc_truth_h0tau1_status;
   int mc_truth_h0tau2_status;
   int mc_truth_h0taul1_status;
   int mc_truth_h0taunutau1_status;
   int mc_truth_h0taunu1_status;
   int mc_truth_h0taul2_status;
   int mc_truth_h0taunutau2_status;
   int mc_truth_h0taunu2_status;

   int mc_truth_h0b1_status;
   int mc_truth_h0b2_status;
   int mc_truth_h0b1_IS_status;
   int mc_truth_h0b2_IS_status;
   
   int mc_truth_t1_status;
   int mc_truth_t2_status;
   int mc_truth_tb1_status;
   int mc_truth_tb2_status;

   int mc_truth_tW1_status;
   int mc_truth_tWnu1_status;
   int mc_truth_tWnutau1_status;
   int mc_truth_tWl1_status;
   int mc_truth_tWtau1_status;
   int mc_truth_tWtaunu1_status;
   int mc_truth_tWtaunutau1_status;
   int mc_truth_tWtaul1_status;
   int mc_truth_tWq11_status;
   int mc_truth_tWq21_status;

   int mc_truth_tW2_status;
   int mc_truth_tWnu2_status;
   int mc_truth_tWnutau2_status;
   int mc_truth_tWl2_status;
   int mc_truth_tWtau2_status;
   int mc_truth_tWtaunu2_status;
   int mc_truth_tWtaunutau2_status;
   int mc_truth_tWtaul2_status;
   int mc_truth_tWq12_status;
   int mc_truth_tWq22_status;

   int mc_truth_j1_status;
   int mc_truth_j2_status;
   int mc_truth_j3_status;

   // tZq
   int mc_truth_tzq_channel;

   // TLV

   TLorentzVector mc_truth_gammal1_p4;
   TLorentzVector mc_truth_gammal2_p4;
   TLorentzVector mc_truth_gammatau1_p4;
   TLorentzVector mc_truth_gammatau2_p4;
   TLorentzVector mc_truth_gammataul1_p4;
   TLorentzVector mc_truth_gammataul2_p4;
   TLorentzVector mc_truth_gammataunu1_p4;
   TLorentzVector mc_truth_gammataunu2_p4;
   TLorentzVector mc_truth_gammataunutau1_p4;
   TLorentzVector mc_truth_gammataunutau2_p4;
   
   TLorentzVector mc_truth_Z_p4;
   TLorentzVector mc_truth_Zl1_p4;
   TLorentzVector mc_truth_Zl2_p4;
   TLorentzVector mc_truth_Ztau1_p4;
   TLorentzVector mc_truth_Ztau2_p4;
   TLorentzVector mc_truth_Ztaul1_p4;
   TLorentzVector mc_truth_Ztaul2_p4;
   TLorentzVector mc_truth_Ztaunu1_p4;
   TLorentzVector mc_truth_Ztaunu2_p4;
   TLorentzVector mc_truth_Ztaunutau1_p4;
   TLorentzVector mc_truth_Ztaunutau2_p4;
   TLorentzVector mc_truth_Zq1_p4;
   TLorentzVector mc_truth_Zq2_p4;
   TLorentzVector mc_truth_Znu1_p4;
   TLorentzVector mc_truth_Znu2_p4;

   TLorentzVector mc_truth_W_p4;
   TLorentzVector mc_truth_Wl_p4;
   TLorentzVector mc_truth_Wnu_p4;
   TLorentzVector mc_truth_Wtau_p4;
   TLorentzVector mc_truth_Wtaunu_p4;
   TLorentzVector mc_truth_Wtaunutau_p4;
   TLorentzVector mc_truth_Wtaul_p4;
   TLorentzVector mc_truth_Wnutau_p4;
   TLorentzVector mc_truth_Wq1_p4;
   TLorentzVector mc_truth_Wq2_p4;
   
   TLorentzVector mc_truth_t_p4;
   TLorentzVector mc_truth_tb_p4;
   TLorentzVector mc_truth_tb_IS_p4;
   TLorentzVector mc_truth_tW_p4;
   TLorentzVector mc_truth_tWnu_p4;
   TLorentzVector mc_truth_tWnutau_p4;
   TLorentzVector mc_truth_tWl_p4;
   TLorentzVector mc_truth_tWtau_p4;
   TLorentzVector mc_truth_tWtaunu_p4;
   TLorentzVector mc_truth_tWtaunutau_p4;
   TLorentzVector mc_truth_tWtaul_p4;
   TLorentzVector mc_truth_tWq1_p4;
   TLorentzVector mc_truth_tWq2_p4;

   // pdgId

   int mc_truth_gammal1_id;
   int mc_truth_gammal2_id;
   int mc_truth_gammatau1_id;
   int mc_truth_gammatau2_id;
   int mc_truth_gammataul1_id;
   int mc_truth_gammataul2_id;
   int mc_truth_gammataunu1_id;
   int mc_truth_gammataunu2_id;
   int mc_truth_gammataunutau1_id;
   int mc_truth_gammataunutau2_id;
   
   int mc_truth_Z_id;
   int mc_truth_Zl1_id;
   int mc_truth_Zl2_id;
   int mc_truth_Ztau1_id;
   int mc_truth_Ztau2_id;
   int mc_truth_Ztaul1_id;
   int mc_truth_Ztaul2_id;
   int mc_truth_Ztaunu1_id;
   int mc_truth_Ztaunu2_id;
   int mc_truth_Ztaunutau1_id;
   int mc_truth_Ztaunutau2_id;
   int mc_truth_Zq1_id;
   int mc_truth_Zq2_id;
   int mc_truth_Znu1_id;
   int mc_truth_Znu2_id;

   int mc_truth_W_id;
   int mc_truth_Wl_id;
   int mc_truth_Wnu_id;
   int mc_truth_Wtau_id;
   int mc_truth_Wtaunu_id;
   int mc_truth_Wtaunutau_id;
   int mc_truth_Wtaul_id;
   int mc_truth_Wnutau_id;
   int mc_truth_Wq1_id;
   int mc_truth_Wq2_id;
   
   int mc_truth_t_id;
   int mc_truth_tb_id;
   int mc_truth_tb_IS_id;
   int mc_truth_tW_id;
   int mc_truth_tWnu_id;
   int mc_truth_tWnutau_id;
   int mc_truth_tWl_id;
   int mc_truth_tWtau_id;
   int mc_truth_tWtaunu_id;
   int mc_truth_tWtaunutau_id;
   int mc_truth_tWtaul_id;
   int mc_truth_tWq1_id;
   int mc_truth_tWq2_id;

   // pt

   float mc_truth_gammal1_pt;
   float mc_truth_gammal2_pt;
   float mc_truth_gammatau1_pt;
   float mc_truth_gammatau2_pt;
   float mc_truth_gammataul1_pt;
   float mc_truth_gammataul2_pt;
   float mc_truth_gammataunu1_pt;
   float mc_truth_gammataunu2_pt;
   float mc_truth_gammataunutau1_pt;
   float mc_truth_gammataunutau2_pt;
   
   float mc_truth_Z_pt;
   float mc_truth_Zl1_pt;
   float mc_truth_Zl2_pt;
   float mc_truth_Ztau1_pt;
   float mc_truth_Ztau2_pt;
   float mc_truth_Ztaul1_pt;
   float mc_truth_Ztaul2_pt;
   float mc_truth_Ztaunu1_pt;
   float mc_truth_Ztaunu2_pt;
   float mc_truth_Ztaunutau1_pt;
   float mc_truth_Ztaunutau2_pt;
   float mc_truth_Zq1_pt;
   float mc_truth_Zq2_pt;
   float mc_truth_Znu1_pt;
   float mc_truth_Znu2_pt;

   float mc_truth_W_pt;
   float mc_truth_Wl_pt;
   float mc_truth_Wnu_pt;
   float mc_truth_Wtau_pt;
   float mc_truth_Wtaunu_pt;
   float mc_truth_Wtaunutau_pt;
   float mc_truth_Wtaul_pt;
   float mc_truth_Wnutau_pt;
   float mc_truth_Wq1_pt;
   float mc_truth_Wq2_pt;
   
   float mc_truth_t_pt;
   float mc_truth_tb_pt;
   float mc_truth_tb_IS_pt;
   float mc_truth_tW_pt;
   float mc_truth_tWnu_pt;
   float mc_truth_tWnutau_pt;
   float mc_truth_tWl_pt;
   float mc_truth_tWtau_pt;
   float mc_truth_tWtaunu_pt;
   float mc_truth_tWtaunutau_pt;
   float mc_truth_tWtaul_pt;
   float mc_truth_tWq1_pt;
   float mc_truth_tWq2_pt;

   // eta

   float mc_truth_gammal1_eta;
   float mc_truth_gammal2_eta;
   float mc_truth_gammatau1_eta;
   float mc_truth_gammatau2_eta;
   float mc_truth_gammataul1_eta;
   float mc_truth_gammataul2_eta;
   float mc_truth_gammataunu1_eta;
   float mc_truth_gammataunu2_eta;
   float mc_truth_gammataunutau1_eta;
   float mc_truth_gammataunutau2_eta;
   
   float mc_truth_Z_eta;
   float mc_truth_Zl1_eta;
   float mc_truth_Zl2_eta;
   float mc_truth_Ztau1_eta;
   float mc_truth_Ztau2_eta;
   float mc_truth_Ztaul1_eta;
   float mc_truth_Ztaul2_eta;
   float mc_truth_Ztaunu1_eta;
   float mc_truth_Ztaunu2_eta;
   float mc_truth_Ztaunutau1_eta;
   float mc_truth_Ztaunutau2_eta;
   float mc_truth_Zq1_eta;
   float mc_truth_Zq2_eta;
   float mc_truth_Znu1_eta;
   float mc_truth_Znu2_eta;

   float mc_truth_W_eta;
   float mc_truth_Wl_eta;
   float mc_truth_Wnu_eta;
   float mc_truth_Wtau_eta;
   float mc_truth_Wtaunu_eta;
   float mc_truth_Wtaunutau_eta;
   float mc_truth_Wtaul_eta;
   float mc_truth_Wnutau_eta;
   float mc_truth_Wq1_eta;
   float mc_truth_Wq2_eta;
   
   float mc_truth_t_eta;
   float mc_truth_tb_eta;
   float mc_truth_tb_IS_eta;
   float mc_truth_tW_eta;
   float mc_truth_tWnu_eta;
   float mc_truth_tWnutau_eta;
   float mc_truth_tWl_eta;
   float mc_truth_tWtau_eta;
   float mc_truth_tWtaunu_eta;
   float mc_truth_tWtaunutau_eta;
   float mc_truth_tWtaul_eta;
   float mc_truth_tWq1_eta;
   float mc_truth_tWq2_eta;
   
   // phi

   float mc_truth_gammal1_phi;
   float mc_truth_gammal2_phi;
   float mc_truth_gammatau1_phi;
   float mc_truth_gammatau2_phi;
   float mc_truth_gammataul1_phi;
   float mc_truth_gammataul2_phi;
   float mc_truth_gammataunu1_phi;
   float mc_truth_gammataunu2_phi;
   float mc_truth_gammataunutau1_phi;
   float mc_truth_gammataunutau2_phi;
   
   float mc_truth_Z_phi;
   float mc_truth_Zl1_phi;
   float mc_truth_Zl2_phi;
   float mc_truth_Ztau1_phi;
   float mc_truth_Ztau2_phi;
   float mc_truth_Ztaul1_phi;
   float mc_truth_Ztaul2_phi;
   float mc_truth_Ztaunu1_phi;
   float mc_truth_Ztaunu2_phi;
   float mc_truth_Ztaunutau1_phi;
   float mc_truth_Ztaunutau2_phi;
   float mc_truth_Zq1_phi;
   float mc_truth_Zq2_phi;
   float mc_truth_Znu1_phi;
   float mc_truth_Znu2_phi;

   float mc_truth_W_phi;
   float mc_truth_Wl_phi;
   float mc_truth_Wnu_phi;
   float mc_truth_Wtau_phi;
   float mc_truth_Wtaunu_phi;
   float mc_truth_Wtaunutau_phi;
   float mc_truth_Wtaul_phi;
   float mc_truth_Wnutau_phi;
   float mc_truth_Wq1_phi;
   float mc_truth_Wq2_phi;
   
   float mc_truth_t_phi;
   float mc_truth_tb_phi;
   float mc_truth_tb_IS_phi;
   float mc_truth_tW_phi;
   float mc_truth_tWnu_phi;
   float mc_truth_tWnutau_phi;
   float mc_truth_tWl_phi;
   float mc_truth_tWtau_phi;
   float mc_truth_tWtaunu_phi;
   float mc_truth_tWtaunutau_phi;
   float mc_truth_tWtaul_phi;
   float mc_truth_tWq1_phi;
   float mc_truth_tWq2_phi;
   
   // E

   float mc_truth_gammal1_E;
   float mc_truth_gammal2_E;
   float mc_truth_gammatau1_E;
   float mc_truth_gammatau2_E;
   float mc_truth_gammataul1_E;
   float mc_truth_gammataul2_E;
   float mc_truth_gammataunu1_E;
   float mc_truth_gammataunu2_E;
   float mc_truth_gammataunutau1_E;
   float mc_truth_gammataunutau2_E;
   
   float mc_truth_Z_E;
   float mc_truth_Zl1_E;
   float mc_truth_Zl2_E;
   float mc_truth_Ztau1_E;
   float mc_truth_Ztau2_E;
   float mc_truth_Ztaul1_E;
   float mc_truth_Ztaul2_E;
   float mc_truth_Ztaunu1_E;
   float mc_truth_Ztaunu2_E;
   float mc_truth_Ztaunutau1_E;
   float mc_truth_Ztaunutau2_E;
   float mc_truth_Zq1_E;
   float mc_truth_Zq2_E;
   float mc_truth_Znu1_E;
   float mc_truth_Znu2_E;

   float mc_truth_W_E;
   float mc_truth_Wl_E;
   float mc_truth_Wnu_E;
   float mc_truth_Wtau_E;
   float mc_truth_Wtaunu_E;
   float mc_truth_Wtaunutau_E;
   float mc_truth_Wtaul_E;
   float mc_truth_Wnutau_E;
   float mc_truth_Wq1_E;
   float mc_truth_Wq2_E;
   
   float mc_truth_t_E;
   float mc_truth_tb_E;
   float mc_truth_tb_IS_E;
   float mc_truth_tW_E;
   float mc_truth_tWnu_E;
   float mc_truth_tWnutau_E;
   float mc_truth_tWl_E;
   float mc_truth_tWtau_E;
   float mc_truth_tWtaunu_E;
   float mc_truth_tWtaunutau_E;
   float mc_truth_tWtaul_E;
   float mc_truth_tWq1_E;
   float mc_truth_tWq2_E;
   
   // status

   int mc_truth_gammal1_status;
   int mc_truth_gammal2_status;
   int mc_truth_gammatau1_status;
   int mc_truth_gammatau2_status;
   int mc_truth_gammataul1_status;
   int mc_truth_gammataul2_status;
   int mc_truth_gammataunu1_status;
   int mc_truth_gammataunu2_status;
   int mc_truth_gammataunutau1_status;
   int mc_truth_gammataunutau2_status;
   
   int mc_truth_Z_status;
   int mc_truth_Zl1_status;
   int mc_truth_Zl2_status;
   int mc_truth_Ztau1_status;
   int mc_truth_Ztau2_status;
   int mc_truth_Ztaul1_status;
   int mc_truth_Ztaul2_status;
   int mc_truth_Ztaunu1_status;
   int mc_truth_Ztaunu2_status;
   int mc_truth_Ztaunutau1_status;
   int mc_truth_Ztaunutau2_status;
   int mc_truth_Zq1_status;
   int mc_truth_Zq2_status;
   int mc_truth_Znu1_status;
   int mc_truth_Znu2_status;
   
   int mc_truth_W_status;
   int mc_truth_Wl_status;
   int mc_truth_Wnu_status;
   int mc_truth_Wtau_status;
   int mc_truth_Wtaunu_status;
   int mc_truth_Wtaunutau_status;
   int mc_truth_Wtaul_status;
   int mc_truth_Wnutau_status;
   int mc_truth_Wq1_status;
   int mc_truth_Wq2_status;
   
   int mc_truth_t_status;
   int mc_truth_tb_status;
   int mc_truth_tb_IS_status;
   int mc_truth_tW_status;
   int mc_truth_tWnu_status;
   int mc_truth_tWnutau_status;
   int mc_truth_tWl_status;
   int mc_truth_tWtau_status;
   int mc_truth_tWtaunu_status;
   int mc_truth_tWtaunutau_status;
   int mc_truth_tWtaul_status;
   int mc_truth_tWq1_status;
   int mc_truth_tWq2_status;

   // gen
   int gen_n;
   std::vector<float> gen_pt;
   std::vector<float> gen_eta;
   std::vector<float> gen_phi;
   std::vector<float> gen_m;
   std::vector<float> gen_E;
   std::vector<int> gen_id;
   std::vector<int> gen_charge;
   std::vector<int> gen_status;
   std::vector<int> gen_index;
   std::vector<int> gen_mother_index;
   std::vector<int> gen_daughter_n;
   std::vector<std::vector<int> > gen_daughter_index;
};

#endif
