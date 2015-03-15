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

   float pv_x;
   float pv_y;
   float pv_z;

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
   std::vector<float> el_eleEoPout;
   std::vector<float> el_PreShowerOverRaw;

   std::vector<float> el_mvaNonTrigV0;

   std::vector<float> el_lepMVA;

   std::vector<float> el_lepMVA_neuRelIso;
   std::vector<float> el_lepMVA_chRelIso;
   std::vector<float> el_lepMVA_jetDR;
   std::vector<float> el_lepMVA_jetPtRatio;
   std::vector<float> el_lepMVA_jetBTagCSV;
   std::vector<float> el_lepMVA_sip3d;
   std::vector<float> el_lepMVA_mvaId;
   std::vector<float> el_lepMVA_innerHits;

   std::vector<int> el_hasMCMatch;
   std::vector<float> el_gen_pt;
   std::vector<float> el_gen_eta;
   std::vector<float> el_gen_phi;
   std::vector<float> el_gen_m;
   std::vector<int> el_gen_status;
   std::vector<int> el_gen_id;
   std::vector<int> el_gen_charge;
   std::vector<float> el_gen_dr;

   std::vector<bool> el_hasMatchedConversion;

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

   std::vector<float> mu_miniIso;

   std::vector<int> mu_isGlobalMuon;
   std::vector<int> mu_isTrackerMuon;
   std::vector<int> mu_isStandAloneMuon;
   std::vector<int> mu_isCaloMuon;
   std::vector<int> mu_isPFMuon;

   std::vector<float> mu_vx;
   std::vector<float> mu_vy;
   std::vector<float> mu_vz;

   std::vector<bool> mu_isTightMuon;

   std::vector<int> mu_hasGlobalTrack;
   std::vector<float> mu_globalTrack_dxy;
   std::vector<float> mu_globalTrack_dz;
   std::vector<float> mu_globalTrack_dxyError;
   std::vector<float> mu_globalTrack_dzError;

   std::vector<int> mu_hasInnerTrack;
   std::vector<float> mu_innerTrack_dxy;
   std::vector<float> mu_innerTrack_dz;
   std::vector<float> mu_innerTrack_dxyError;
   std::vector<float> mu_innerTrack_dzError;

   std::vector<float> mu_bestTrack_dxy;
   std::vector<float> mu_bestTrack_dz;
   std::vector<float> mu_bestTrack_dxyError;
   std::vector<float> mu_bestTrack_dzError;

   std::vector<float> mu_innerTrack_pt;
   std::vector<float> mu_innerTrack_ptError;

   std::vector<int> mu_numberOfMatches;
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
   std::vector<int> jet_flavour;

   std::vector<float> jet_neutralHadronEnergy;
   std::vector<float> jet_neutralEmEnergy;
   std::vector<float> jet_chargedHadronEnergy;
   std::vector<float> jet_chargedEmEnergy;
   std::vector<float> jet_electronEnergy;
   std::vector<float> jet_muonEnergy;
   std::vector<float> jet_photonEnergy;

   std::vector<float> jet_pileupJetId;

   std::vector<float> jet_gen_pt;
   std::vector<float> jet_gen_eta;
   std::vector<float> jet_gen_phi;
   std::vector<float> jet_gen_m;
   std::vector<float> jet_gen_E;

   std::vector<int> jet_gen_status;
   std::vector<int> jet_gen_id;

   // ttH
   int mc_truth_tth_channel;

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
   std::vector<int> gen_id;
   std::vector<int> gen_charge;
   std::vector<int> gen_status;
   std::vector<int> gen_index;
   std::vector<int> gen_mother_index;
   std::vector<int> gen_daughter_n;
   std::vector<std::vector<int> > gen_daughter_index;
};

#endif
