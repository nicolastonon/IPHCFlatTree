#ifndef HELPER_H
#define HELPER_H 

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

namespace SelfVetoPolicy
{
   enum SelfVetoPolicy
     {
	selfVetoNone=0, selfVetoAll=1, selfVetoFirst=2
     };
}

float GetDeltaR(float,float,float,float);

float ptRelElec(const pat::Electron& elec,const pat::Jet& jet);
float ptRelMuon(const pat::Muon& muon,const pat::Jet& jet);

float ElecPfIsoCharged(const pat::Electron& elec,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR);
float ElecPfIsoNeutral(const pat::Electron& elec,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR);

float MuonPfIsoCharged(const pat::Muon& muon,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR);
float MuonPfIsoNeutral(const pat::Muon& muon,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR);

float isoSumRaw(const std::vector<const pat::PackedCandidate *> & cands, const reco::Candidate &cand, float dR, float innerR, float threshold, SelfVetoPolicy::SelfVetoPolicy selfVeto, int pdgId=-1);

double getPFIsolation(edm::Handle<pat::PackedCandidateCollection>,
		      const reco::Candidate*,
		      double,double,double,
		      bool,bool);

#endif
