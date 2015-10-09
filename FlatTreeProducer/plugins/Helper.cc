#include "IPHCFlatTree/FlatTreeProducer/interface/Helper.hh"

#include "TMath.h"

namespace 
{   
   struct ByEta 
     {	
	bool operator()(const pat::PackedCandidate *c1, const pat::PackedCandidate *c2) const 
	  {	     
	     return c1->eta() < c2->eta();
	  }	
	bool operator()(float c1eta, const pat::PackedCandidate *c2) const 
	  {	     
	     return c1eta < c2->eta();
	  }	
	bool operator()(const pat::PackedCandidate *c1, float c2eta) const 
	  {	     
	     return c1->eta() < c2eta;
	  }	
     };
}

float GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
   float DeltaPhi = TMath::Abs(phi2 - phi1);
      if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}

// part of miniIso for ttH
float isoSumRaw(const std::vector<const pat::PackedCandidate *> & cands, const reco::Candidate &cand, float dR, float innerR, float threshold, SelfVetoPolicy::SelfVetoPolicy selfVeto, int pdgId)
{
   std::vector<const reco::Candidate *> vetos_;
   
   float dR2 = dR*dR, innerR2 = innerR*innerR;
   
   std::vector<const reco::Candidate *> vetos(vetos_);
   for( unsigned int i=0,n=cand.numberOfSourceCandidatePtrs();i<n;++i )
     {
	if(selfVeto == SelfVetoPolicy::selfVetoNone) break;
	const reco::CandidatePtr &cp = cand.sourceCandidatePtr(i);
	if( cp.isNonnull() && cp.isAvailable() )
	  {
	     vetos.push_back(&*cp);
	     if (selfVeto == SelfVetoPolicy::selfVetoFirst) break;
	  }
     }   
   
   typedef std::vector<const pat::PackedCandidate *>::const_iterator IT;
   IT candsbegin = std::lower_bound(cands.begin(), cands.end(), cand.eta() - dR, ByEta());
   IT candsend = std::upper_bound(candsbegin, cands.end(), cand.eta() + dR, ByEta());
   
   double isosum = 0;
   for( IT icharged=candsbegin;icharged<candsend;++icharged )
     {
	// pdgId
	if( pdgId > 0 && abs((*icharged)->pdgId()) != pdgId ) continue;
	// threshold
	if( threshold > 0 && (*icharged)->pt() < threshold ) continue;
	// cone
	float mydr2 = reco::deltaR2(**icharged, cand);
	if( mydr2 >= dR2 || mydr2 < innerR2 ) continue;
	// veto
	if( std::find(vetos.begin(), vetos.end(), *icharged) != vetos.end() )
	  {
	     continue;	     
	  }
	// add to sum
	isosum += (*icharged)->pt();
     }
   return isosum;
}

float ElecPfIsoCharged(const pat::Electron& elec,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR)
{
   std::vector<const pat::PackedCandidate *> charged;
   
   for( const pat::PackedCandidate &p : *pfcands )
     {
	if( p.charge() != 0 )
	  {
	     if( abs(p.pdgId()) == 211 )
	       {
		  if (p.fromPV() > 1 && fabs(p.dz()) < 9999. )
		    {
		       charged.push_back(&p);
		    }		  
	       }	     
	  }	
     }   
	
   return isoSumRaw(charged,elec,miniIsoR,0.0001,0.0,SelfVetoPolicy::selfVetoAll);
}

float MuonPfIsoCharged(const pat::Muon& muon,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR)
{
   std::vector<const pat::PackedCandidate *> charged;
   
   for( const pat::PackedCandidate &p : *pfcands )
     {
	if( p.charge() != 0 )
	  {
	     if( abs(p.pdgId()) == 211 )
	       {
		  if (p.fromPV() > 1 && fabs(p.dz()) < 9999. )
		    {
		       charged.push_back(&p);
		    }		  
	       }	     
	  }	
     }   
	
   return isoSumRaw(charged,muon,miniIsoR,0.0001,0.0,SelfVetoPolicy::selfVetoAll);
}

float ElecPfIsoNeutral(const pat::Electron& elec,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR)
{
   std::vector<const pat::PackedCandidate *> neutral;
   
   for( const pat::PackedCandidate &p : *pfcands )
     {
	if( p.charge() == 0 )
	  {
	     neutral.push_back(&p);
	  }
     }   
   
   return isoSumRaw(neutral,elec,miniIsoR,0.01,0.5,SelfVetoPolicy::selfVetoAll);
}

float MuonPfIsoNeutral(const pat::Muon& muon,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR)
{
   std::vector<const pat::PackedCandidate *> neutral;
   
   for( const pat::PackedCandidate &p : *pfcands )
     {
	if( p.charge() == 0 )
	  {
	     neutral.push_back(&p);
	  }
     }   
   
   return isoSumRaw(neutral,muon,miniIsoR,0.01,0.5,SelfVetoPolicy::selfVetoAll);
}

// https://twiki.cern.ch/twiki/bin/view/CMS/MiniIsolationSUSY
// https://github.com/manuelfs/CfANtupler/blob/master/minicfa/interface/miniAdHocNTupler.h
double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
		      const reco::Candidate* ptcl,
		      double r_iso_min, double r_iso_max, double kt_scale,
		      bool use_pfweight, bool charged_only) 
{
   if (ptcl->pt()<5.) return 99999.;
   double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
   if(ptcl->isElectron()) 
     {	
	if (fabs(ptcl->eta())>1.479) 
	  {
	     deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;
	  }	
     }
    else if(ptcl->isMuon()) 
     {	
	deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
     }
    else 
     {	
	//deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
     }
   
   double iso_nh(0.); double iso_ch(0.);
   double iso_ph(0.); double iso_pu(0.);
   double ptThresh(0.5);
   if(ptcl->isElectron()) ptThresh = 0;
   double r_iso = std::max(r_iso_min,std::min(r_iso_max, kt_scale/ptcl->pt()));
   for (const pat::PackedCandidate &pfc : *pfcands)
     {	
	if (abs(pfc.pdgId())<7) continue;
	double dr = deltaR(pfc, *ptcl);
	if (dr > r_iso) continue;
	////////////////// NEUTRALS /////////////////////////
	if (pfc.charge()==0)
	  {	     
	     if (pfc.pt()>ptThresh) 
	       {		  
		  double wpf(1.);
		  if (use_pfweight)
		    {		       
		       double wpv(0.), wpu(0.);
		       for (const pat::PackedCandidate &jpfc : *pfcands) 
			 {			    
			    double jdr = deltaR(pfc, jpfc);
			    if (pfc.charge()!=0 || jdr<0.00001) continue;
			    double jpt = jpfc.pt();
			    if (pfc.fromPV()>1) wpv *= jpt/jdr;
			    else wpu *= jpt/jdr;
			 }
		       
		       wpv = log(wpv);
		       wpu = log(wpu);
		       wpf = wpv/(wpv+wpu);
		    }
		  
		  /////////// PHOTONS ////////////
		  if (abs(pfc.pdgId())==22) 
		    {		       
		       if(dr < deadcone_ph) continue;
		       iso_ph += wpf*pfc.pt();
		       /////////// NEUTRAL HADRONS ////////////
		    }
		   else if (abs(pfc.pdgId())==130) 
		    {		       
		       if(dr < deadcone_nh) continue;
		       iso_nh += wpf*pfc.pt();
		    }		  
	       }	     
	     ////////////////// CHARGED from PV /////////////////////////
	  }
	 else if (pfc.fromPV()>1)
	  {	     
	     if (abs(pfc.pdgId())==211) 
	       {		  
		  if(dr < deadcone_ch) continue;
		  iso_ch += pfc.pt();
	       }
	     
	     ////////////////// CHARGED from PU /////////////////////////
	  }
	 else 
	  {	     
	     if (pfc.pt()>ptThresh)
	       {		  
		  if(dr < deadcone_pu) continue;
		  iso_pu += pfc.pt();
	       }	     
	  }	
     }
   
   double iso(0.);
   if (charged_only)
     {	
	iso = iso_ch;
     }
    else 
     {	
	iso = iso_ph + iso_nh;
	if (!use_pfweight) iso -= 0.5*iso_pu;
	if (iso>0) iso += iso_ch;
	else iso = iso_ch;
     }
   
   iso = iso/ptcl->pt();
   return iso;
}

float ptRelElec(const pat::Electron& elec,const pat::Jet& jet)
{
   float j_x = jet.px();
   float j_y = jet.py();
   float j_z = jet.pz();

   float l_x = elec.px();
   float l_y = elec.py();
   float l_z = elec.pz();
   
   float j2 = j_x*j_x+j_y*j_y+j_z*j_z;
   float l2 = l_x*l_x+l_y*l_y+l_z*l_z;
   
   float lXj = l_x*j_x+l_y*j_y+l_z*j_z;
   
   float pLrel2 = lXj*lXj/j2;
   
   float pTrel2 = l2-pLrel2;
   
   return (pTrel2 > 0) ? std::sqrt(pTrel2) : 0.0;
}

float ptRelMuon(const pat::Muon& muon,const pat::Jet& jet)
{
   float j_x = jet.px();
   float j_y = jet.py();
   float j_z = jet.pz();

   float l_x = muon.px();
   float l_y = muon.py();
   float l_z = muon.pz();
   
   float j2 = j_x*j_x+j_y*j_y+j_z*j_z;
   float l2 = l_x*l_x+l_y*l_y+l_z*l_z;
   
   float lXj = l_x*j_x+l_y*j_y+l_z*j_z;
   
   float pLrel2 = lXj*lXj/j2;
   
   float pTrel2 = l2-pLrel2;
   
   return (pTrel2 > 0) ? std::sqrt(pTrel2) : 0.0;
}
