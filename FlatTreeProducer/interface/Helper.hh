#ifndef HELPER_H
#define HELPER_H 

#include "DataFormats/PatCandidates/interface/Electron.h"

float GetDeltaR(float,float,float,float);

double getPFIsolation(edm::Handle<pat::PackedCandidateCollection>,
		      const reco::Candidate*,
		      double,double,double,
		      bool,bool);

#endif
