#include <memory>
#include <iostream>
#include <sstream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "IPHCFlatTree/FlatTreeProducer/interface/tinyxml2.h"

#include "IPHCFlatTree/FlatTreeProducer/interface/FlatTree.hh"
#include "IPHCFlatTree/FlatTreeProducer/interface/MCTruth.hh"

#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimator.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TFile.h"
#include "TThreadSlots.h"
#include "TROOT.h"
#include "Compression.h"

using namespace tinyxml2;

class FlatTreeProducer : public edm::EDAnalyzer
{
 public:
   explicit FlatTreeProducer(const edm::ParameterSet&);
   ~FlatTreeProducer();
   
   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   
 private:
   virtual void beginJob() override;
   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
   virtual void endJob() override;

   TMVA::Reader* BookLeptonMVAReader(string basePath, string weightFileName, string type);

   void KeepEvent();
   bool isFloat(const std::string& s);
   bool isFloat(const boost::any & operand);
   bool isBool(const boost::any & operand);
   bool isInt(const boost::any & operand);
   void fillCutVector(const char* cut_type, std::string& cut_value, std::map<std::string, boost::any>& vec);
   void AddValue(const std::string& name);
   void ReadConfFile(const std::string& confFile);
   int CheckAlgo(const std::map<std::string, boost::any>& jet_algo, const char* name, std::string& algo);
   std::string CheckAlgos();

   template <typename T>
     boost::any comp(const std::string& op, T v1, T v2);
   template <typename T>
     void CheckVectorCut(const std::vector<T>& v, const std::string& name);   

   template <typename T>
     void CompareAlgo(const std::string& algo, T conf_algo_value);
   template <typename T>
     void CheckJet(const std::vector<T>& vJetPt, const std::string& jetpt, const std::vector<T>& vJetEta, const std::string& jeteta, const std::string& algo);
   template <typename T>
     void CheckElectron(const std::vector<T>& vElPt, const std::string& elpt, const std::vector<T>& vElEta, const std::string& eleta);
   template <typename T>
     void CheckMuon(const std::vector<T>& vMuonPt, const std::string& muonpt, const std::vector<T>& vMuonEta, const std::string& muoneta);
   
   FlatTree* ftree;
   const edm::Service<TFileService> fs;
   
   TH1D* hcount;
   TH1D* hweight;

   EGammaMvaEleEstimator* elecMVA;
   std::vector<std::string> elecMVACatWeights;

   TMVA::Reader* mu_reader_high_b;
   TMVA::Reader* mu_reader_high_e;
   TMVA::Reader* mu_reader_low_b;
   TMVA::Reader* mu_reader_low_e;
   TMVA::Reader* ele_reader_high_cb;
   TMVA::Reader* ele_reader_high_fb;
   TMVA::Reader* ele_reader_high_ec;
   TMVA::Reader* ele_reader_low_cb;
   TMVA::Reader* ele_reader_low_fb;
   TMVA::Reader* ele_reader_low_ec;

   float lepMVA_neuRelIso;
   float lepMVA_chRelIso;
   float lepMVA_jetDR;
   float lepMVA_jetPtRatio;
   float lepMVA_jetBTagCSV;
   float lepMVA_sip3d;
   
   float lepMVA_mvaId;
   float lepMVA_innerHits;

   float lepMVA_dxy;
   float lepMVA_dz;
   
   XMLDocument xmlconf;
   
   std::string dataFormat_;
   bool isData_;
   
   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
   edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
   
   edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
   edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
   edm::EDGetTokenT<pat::MuonCollection> muonToken_;
   edm::EDGetTokenT<pat::JetCollection> jetToken_;
   edm::EDGetTokenT<std::vector<pat::MET> > metTokenAOD_;
   edm::EDGetTokenT<pat::METCollection> metTokenMINIAOD_;
   edm::EDGetTokenT<double> rhoToken_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
};

bool FlatTreeProducer::isInt(const boost::any & operand)
{
   return operand.type() == typeid(int);
}
bool FlatTreeProducer::isFloat(const boost::any & operand)
{
   return operand.type() == typeid(float);
}
bool FlatTreeProducer::isBool(const boost::any& operand)
{
   return operand.type() == typeid(bool);
}

template<typename T>
  boost::any FlatTreeProducer::comp(const std::string& op, T v1, T v2)
{
   if( !op.compare("||") )
     {
	return v1 || v2;
     }
   else if( !op.compare("&&") )
     {
	return v1 && v2;
     }
   
   if( v1 == typeid(int) && v2 == typeid(int) )
     {
	try
	  {
	     if( !op.compare("|") )
	       {
		  return v1 | v2;
	       }
	     else if( !op.compare("&") )
	       {
		  return v1 & v2;
	       }
	  }
	catch (...)
	  {
	  }
     }
   return "";
}

void FlatTreeProducer::AddValue(const std::string& name)
{
   if( !name.compare("njets") )
     {
	ftree->njets += 1;
     }
   else if( !name.compare("nelectron") )
     {
	ftree->nelectron += 1;
     }
}

template <typename T>
  void FlatTreeProducer::CheckVectorCut(const std::vector<T>& v, const std::string& name)
{
   if( ftree->keep_conf.find(name) != ftree->keep_conf.end() )
     {
	std::map<std::string, std::map<std::string, boost::any> > keep_conf = ftree->keep_conf;
	std::map<std::string, boost::any> map_conf = keep_conf[name];
	for( unsigned int i=0;i<v.size();++i )
	  {
	     if( isInt(map_conf["cut_min"]) && isInt(map_conf["cut_max"]) )
	       {
		  if( v[i] < boost::any_cast<int>(map_conf["cut_min"]) || v[i] > boost::any_cast<int>(map_conf["cut_max"]) )
		    {
		       AddValue(name);
		    }
	       }
	     else if( isFloat(map_conf["cut_min"]) && isFloat(map_conf["cut_max"]) )
	       {
		  if( v[i] < boost::any_cast<float>(map_conf["cut_min"]) || v[i] > boost::any_cast<float>(map_conf["cut_max"]) )
		    {
		       AddValue(name);
		    }
	       }
	     else{std::cout << "Wrong types" << std::endl;}
	  }
     }
}

template <typename T>
  void FlatTreeProducer::CompareAlgo(const std::string& algo, T conf_algo_value)
{
   if( !algo.compare("jet_JBP") )
     {
	for( unsigned int i=0;i<ftree->jet_JBP.size();++i )
	  if( ftree->jet_JBP[i] > conf_algo_value )
	    AddValue("njets");
     }
   else if( !algo.compare("jet_JP") )
     {
	for( unsigned int i=0;i<ftree->jet_JP.size();++i )
	  if( ftree->jet_JP[i] > conf_algo_value )
	    AddValue("njets");
     }
   else if( !algo.compare("jet_TCHP") )
     {
	for( unsigned int i=0;i<ftree->jet_TCHP.size();++i )
	  if( ftree->jet_TCHP[i] > conf_algo_value )
	    AddValue("njets");
     }
   else if( !algo.compare("jet_TCHE") )
     {
	for( unsigned int i=0;i<ftree->jet_TCHE.size();++i )
	  if( ftree->jet_TCHE[i] > conf_algo_value )
	    AddValue("njets");
     }
   else if( !algo.compare("jet_SSVHP") )
     {
	for( unsigned int i=0;i<ftree->jet_SSVHP.size();++i )
	  if( ftree->jet_SSVHP[i] > conf_algo_value )
	    AddValue("njets");
     }
   else if( !algo.compare("jet_SSVHE") )
     {
	for( unsigned int i=0;i<ftree->jet_SSVHE.size();++i )
	  if( ftree->jet_SSVHE[i] > conf_algo_value )
	    AddValue("njets");
     }
   else if( !algo.compare("jet_CMVA") )
     {
	for( unsigned int i=0;i<ftree->jet_CMVA.size();++i )
	  if( ftree->jet_CMVA[i] > conf_algo_value )
	    AddValue("njets");
     }
   else if( !algo.compare("jet_CSV") )
     {
	for( unsigned int i=0;i<ftree->jet_CSV.size();++i )
	  if( ftree->jet_CSV[i] > conf_algo_value )
	    AddValue("njets");
     }
   else if( !algo.compare("jet_CSVv2") )
     {
	for( unsigned int i=0;i<ftree->jet_CSVv2.size();++i )
	  if( ftree->jet_CSVv2[i] > conf_algo_value )
	    AddValue("njets");
     }
   else if( !algo.compare("jet_flavour") )
     {
	for( unsigned int i=0;i<ftree->jet_flavour.size();++i )
	  if( ftree->jet_flavour[i] > conf_algo_value )
	    AddValue("njets");
     }
}

template <typename T>
void FlatTreeProducer::CheckJet(const std::vector<T>& vJetPt, const std::string& jetpt, const std::vector<T>& vJetEta, const std::string& jeteta, const std::string& algo)
{
   std::map<std::string, std::map<std::string, boost::any> > keep_conf = ftree->keep_conf;
   if( keep_conf.find(jetpt) != keep_conf.end() )
     {
	std::map<std::string, boost::any> conf_pt = keep_conf[jetpt];
	for( unsigned int i=0;i<vJetPt.size();++i )
	  {
	     if( isInt(conf_pt["cut_min"]) )
	       {
		  if( vJetPt[i] > boost::any_cast<int>(conf_pt["cut_min"]) )
		    {
		       if( keep_conf.find(jeteta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[jeteta];
			    for( unsigned int i=0;i<vJetEta.size();++i )
			      {
				 if( isInt(conf_eta["cut_max"]) )
				   {
				      if( vJetEta[i] < boost::any_cast<int>(conf_eta["cut_max"]) )
					{
					   if( !algo.empty() )
					     {
						std::map<std::string, boost::any> conf_algo = keep_conf[algo];
						if( isInt(conf_algo["cut_algo"]) )
						  {
						     CompareAlgo(algo, boost::any_cast<int>(conf_algo["cut_algo"]));
						  }
						else
						  {
						     std::cout << "'jet_pt' : 'cut_min' and 'cut_max' are type 'int', 'cut_algo' cannot be a 'float'." << std::endl;
						  }
					     }
					   else
					     {
						AddValue("njets");
					     }
					}
				   }
				 else
				   {
				      std::cout << "'jet_pt' : 'cut_min' is a type 'int', 'jet_eta' : 'cut_max' cannot be a 'float'." << std::endl;
				      break;
				   }
			      }
			 }
		       else
			 {
			    std::cout << "'jet_pt' : 'cut_min' set, but not 'jet_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else if( isFloat(conf_pt["cut_min"]) )
	       {
		  if( vJetPt[i] > boost::any_cast<float>(conf_pt["cut_min"]) )
		    {
		       if( keep_conf.find(jeteta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[jeteta];
			    for( unsigned int i=0;i<vJetEta.size();++i )
			      {
				 if( isFloat(conf_eta["cut_max"]) )
				   {
				      if( vJetEta[i] < boost::any_cast<float>(conf_eta["cut_max"]) )
					{
					   if( !algo.empty() )
					     {
						std::map<std::string, boost::any> conf_algo = keep_conf[algo];
						if( isFloat(conf_algo["cut_algo"]) )
						  {
						     CompareAlgo(algo, boost::any_cast<float>(conf_algo["cut_algo"]));
						  }
						else
						  {
						     std::cout << "'jet_pt' : 'cut_min' and 'cut_max' are type 'float', 'cut_algo' cannot be an 'int'." << std::endl;
						  }
					     }
					   else
					     {
						AddValue("njets");
					     }
					}
				   }
				 else
				   {
				      std::cout << "'jet_pt' : 'cut_min' is a type 'float', 'jet_eta' : 'cut_max' cannot be an 'int'." << std::endl;
				      break;
				   }
			      }
			 }
		       else
			 {
			    std::cout << "'jet_pt' : 'cut_min' set, but not 'jet_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else{std::cout << "'jet_pt' : wrong types." << std::endl;}
	  }
     }
}

template <typename T>
  void FlatTreeProducer::CheckElectron(const std::vector<T>& vElPt, const std::string& elpt, const std::vector<T>& vElEta, const std::string& eleta)
{
   std::map<std::string, std::map<std::string, boost::any> > keep_conf = ftree->keep_conf;
   if( keep_conf.find(elpt) != keep_conf.end() )
     {
	std::map<std::string, boost::any> conf_pt = keep_conf[elpt];
	for( unsigned int i=0;i<vElPt.size();++i )
	  {
	     if( isInt(conf_pt["cut_min"]) )
	       {
		  if( vElPt[i] > boost::any_cast<int>(conf_pt["cut_min"]) )
		    {
		       if( keep_conf.find(eleta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[eleta];
			    for( unsigned int i=0;i<vElEta.size();++i )
			      {
				 if( isInt(conf_eta["cut_max"]) )
				   {
				      if( vElEta[i] < boost::any_cast<int>(conf_eta["cut_max"]) )
					{
					   AddValue("nelectron");
					}
				   }
				 else
				   {
				      std::cout << "'el_pt' : 'cut_min' is a type 'int', 'el_eta' : 'cut_max' cannot be a 'float'." << std::endl;
				      break;
				   }
			      }
			 }
		       else
			 {
			    std::cout << "'el_pt' : 'cut_min' set, but not 'el_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else if( isFloat(conf_pt["cut_min"]) )
	       {
		  if( vElPt[i] > boost::any_cast<float>(conf_pt["cut_min"]) )
		    {
		       if( keep_conf.find(eleta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[eleta];
			    for( unsigned int i=0;i<vElEta.size();++i )
			      {
				 if( isFloat(conf_eta["cut_max"]) )
				   {
				      if( vElEta[i] < boost::any_cast<float>(conf_eta["cut_max"]) )
					{
					   AddValue("nelectron");
					}
				   }
				 else
				   {
				      std::cout << "'el_pt' : 'cut_min' is a type 'float', 'el_eta' : 'cut_max' cannot be an 'int'." << std::endl;
				      break;
				   }
			      }
			 }
		       else
			 {
			    std::cout << "'el_pt' : 'cut_min' set, but not 'el_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else{std::cout << "'el_pt' : wrong types." << std::endl;}
	  }
     }
}

template <typename T>
  void FlatTreeProducer::CheckMuon(const std::vector<T>& vMuonPt, const std::string& muonpt, const std::vector<T>& vMuonEta, const std::string& muoneta)
{
   std::map<std::string, std::map<std::string, boost::any> > keep_conf = ftree->keep_conf;
   if( keep_conf.find(muonpt) != keep_conf.end() )
     {
	std::map<std::string, boost::any> conf_pt = keep_conf[muonpt];
	for( unsigned int i=0;i<vMuonPt.size();++i )
	  {
	     if( isInt(conf_pt["cut_min"]) )
	       {
		  if( vMuonPt[i] > boost::any_cast<int>(conf_pt["cut_min"]) )
		    {
		       if( keep_conf.find(muoneta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[muoneta];
			    for( unsigned int i=0;i<vMuonEta.size();++i )
			      {
				 if( isInt(conf_eta["cut_max"]) )
				   {
				      if( vMuonEta[i] < boost::any_cast<int>(conf_eta["cut_max"]) )
					{
					   AddValue("nelectron");
					}
				   }
				 else
				   {
				      std::cout << "'mu_pt' : 'cut_min' is a type 'int', 'mu_eta' : 'cut_max' cannot be a 'float'." << std::endl;
				      break;
				   }
			      }
			 }
		       else
			 {
			    std::cout << "'mu_pt' : 'cut_min' set, but not 'mu_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else if( isFloat(conf_pt["cut_min"]) )
	       {
		  if( vMuonPt[i] > boost::any_cast<float>(conf_pt["cut_min"]) )
		    {
		       if( keep_conf.find(muoneta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[muoneta];
			    for( unsigned int i=0;i<vMuonEta.size();++i )
			      {
				 if( isFloat(conf_eta["cut_max"]) )
				   {
				      if( vMuonEta[i] < boost::any_cast<float>(conf_eta["cut_max"]) )
					{
					   AddValue("nelectron");
					}
				   }
				 else
				   {
				      std::cout << "'mu_pt' : 'cut_min' is a type 'float', 'mu_eta' : 'cut_max' cannot be an 'int'." << std::endl;
				      break;
				   }
			      }
			 }
		       else
			 {
			    std::cout << "'mu_pt' : 'cut_min' set, but not 'mu_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else{std::cout << "'mu_pt' : wrong types." << std::endl;}
	  }
     }
}

int FlatTreeProducer::CheckAlgo(const std::map<std::string, boost::any>& jet_algo, const char* name, std::string& algo)
{
   if( jet_algo.count("cut_algo") == 1 )
     {
	algo = name;
	return 1;
     }
   return 0;
}

std::string FlatTreeProducer::CheckAlgos()
{
   int nbAlgo = 0;
   std::string algo("");
   std::map<std::string, std::map<std::string, boost::any> > keep_conf = ftree->keep_conf;
   std::map<std::string, boost::any> jet_JBP = keep_conf["jet_JBP"];
   std::map<std::string, boost::any> jet_JP = keep_conf["jet_JP"];
   std::map<std::string, boost::any> jet_TCHP = keep_conf["jet_TCHP"];
   std::map<std::string, boost::any> jet_TCHE = keep_conf["jet_TCHE"];
   std::map<std::string, boost::any> jet_SSVHP = keep_conf["jet_SSVHP"];
   std::map<std::string, boost::any> jet_SSVHE = keep_conf["jet_SSVHE"];
   std::map<std::string, boost::any> jet_CMVA = keep_conf["jet_CMVA"];
   std::map<std::string, boost::any> jet_CSV = keep_conf["jet_CSV"];
   std::map<std::string, boost::any> jet_CSVv2 = keep_conf["jet_CSVv2"];
   std::map<std::string, boost::any> jet_flavour = keep_conf["jet_flavour"];
   
   nbAlgo += CheckAlgo(jet_JBP, "jet_JBP", algo);
   nbAlgo += CheckAlgo(jet_JP, "jet_JP", algo);
   nbAlgo += CheckAlgo(jet_TCHP, "jet_TCHP", algo);
   nbAlgo += CheckAlgo(jet_TCHE, "jet_TCHE", algo);
   nbAlgo += CheckAlgo(jet_SSVHP, "jet_SSVHP", algo);
   nbAlgo += CheckAlgo(jet_SSVHE, "jet_SSHVE", algo);
   nbAlgo += CheckAlgo(jet_CMVA, "jet_CMVA", algo);
   nbAlgo += CheckAlgo(jet_CSV, "jet_CSV", algo);
   nbAlgo += CheckAlgo(jet_CSVv2, "jet_CSVv2", algo);
   nbAlgo += CheckAlgo(jet_flavour, "jet_flavour", algo);
   
   if( nbAlgo > 1 )
     {
	std::cout << "Different algorithms are set, please choose only one. Algorithm not considered." << std::endl;
	algo.clear();
	return algo;
     }
   return algo;
}

void FlatTreeProducer::KeepEvent()
{
   // Jets
   std::string algo = this->CheckAlgos();
   // jet_pt > cut_min && jet_eta >
   this->CheckJet(ftree->jet_pt, "jet_pt", ftree->jet_eta, "jet_eta", algo);
   
   // Electron
   // el_pt> XX && fabs(el_eta)> YY && iso < ZZ
   // TODO ISO
   this->CheckElectron(ftree->el_pt, "el_pt", ftree->el_eta, "el_eta");
   
   // Muon idem
   this->CheckMuon(ftree->mu_pt, "mu_pt", ftree->mu_eta, "mu_eta");
}

bool FlatTreeProducer::isFloat(const std::string& s)
{
   return !std::all_of(s.begin(), s.end(), ::isdigit);
}

void FlatTreeProducer::fillCutVector(const char* cut_type, std::string& cut_value, std::map<std::string, boost::any> & vmap)
{
   if( !cut_value.empty() )
     {
	bool t_cut = isFloat(cut_value);
	if( t_cut )
	  {
	     float fcut = atof(cut_value.c_str());
	     vmap[cut_type] = fcut;
	  }
	else if( !t_cut )
	  {
	     int fcut = atoi(cut_value.c_str());
	     vmap[cut_type] = fcut;
	  }
	else
	  {
	     std::cout << "Warning ! Different types of cut (int or float). Cut Skipped." << std::endl;
	  }
     }
   cut_value = "";
}

void FlatTreeProducer::ReadConfFile(const std::string& confFile)
{
   xmlconf.LoadFile(confFile.c_str());
   XMLElement* tElement = xmlconf.FirstChildElement("var");
   
   std::string vcutmin("");
   std::string vcutmax("");
   std::string vcutiso("");
   std::string vcutalgo("");
   std::string vcutexpr("");
   for( XMLElement* child=tElement;child!=0;child=child->NextSiblingElement() )
     {
	const std::string& vname = child->ToElement()->Attribute("name");
	const std::string& vsave = child->ToElement()->Attribute("save");
	if (child->ToElement()->Attribute("cut_min"))
	  {
	     vcutmin = child->ToElement()->Attribute("cut_min");
	  }
	if (child->ToElement()->Attribute("cut_max"))
	  {
	     vcutmax = child->ToElement()->Attribute("cut_max");
	  }
	if (child->ToElement()->Attribute("cut_iso"))
	  {
	     vcutiso = child->ToElement()->Attribute("cut_iso");
	  }
	if (child->ToElement()->Attribute("cut_expr"))
	  {
	     ;
	  }
	if (child->ToElement()->Attribute("cut_algo"))
	  {
	     vcutalgo = child->ToElement()->Attribute("cut_algo");
	  }
	
	std::map<std::string, boost::any > vmap;
	bool bsave = atoi(vsave.c_str());
	fillCutVector("cut_min", vcutmin, vmap);
	fillCutVector("cut_max", vcutmax, vmap);
	fillCutVector("cut_iso", vcutiso, vmap);
	fillCutVector("cut_algo", vcutalgo, vmap);
	if (vmap.size() > 0)
	  ftree->keep_conf.insert(std::make_pair(vname, vmap));
	ftree->conf.insert(std::make_pair(vname,bsave));
     }
}

TMVA::Reader* FlatTreeProducer::BookLeptonMVAReader(string basePath, string weightFileName, string type)
{
   TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");
   
   reader->AddVariable("relIso-chargedIso/pt",&lepMVA_neuRelIso);
   reader->AddVariable("chargedIso/pt",       &lepMVA_chRelIso);
   reader->AddVariable("min(dr_in,0.5)",      &lepMVA_jetDR);
   reader->AddVariable("min(ptf_in,1.5)",     &lepMVA_jetPtRatio);
   reader->AddVariable("max(CSV_in,0)",       &lepMVA_jetBTagCSV);
   reader->AddVariable("sip3d",               &lepMVA_sip3d);

   if (type == "ele")
   {
       reader->AddVariable("mvaId",           &lepMVA_mvaId);
       reader->AddVariable("innerHits",       &lepMVA_innerHits);
   }
   if (type == "mu")
   {
       reader->AddVariable("log(abs(dxy))",   &lepMVA_dxy);
       reader->AddVariable("log(abs(dz))",    &lepMVA_dz);
   }

   reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

   return reader;
}

FlatTreeProducer::FlatTreeProducer(const edm::ParameterSet& iConfig)
{
   // ###
   // Temporarily redirecting stdout to avoid huge TMVA loading dump
   // ###
   cout << "Temporarily redirecting stdout to avoid huge TMVA loading dump..." << endl;
   stringstream tmpBuffer;
   streambuf* oldStdout = cout.rdbuf(tmpBuffer.rdbuf());

   // ###############
   // #  Load MVAs  #
   // ###############

   string CMSSW_BASE(getenv("CMSSW_BASE")); 
   elecMVA = new EGammaMvaEleEstimator();
   
   string EGammaElectronToolsPath = CMSSW_BASE+"/src/EgammaAnalysis/ElectronTools/";
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/Electrons_BDTG_NonTrigV0_Cat1.weights.xml");
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/Electrons_BDTG_NonTrigV0_Cat2.weights.xml");
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/Electrons_BDTG_NonTrigV0_Cat3.weights.xml");
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/Electrons_BDTG_NonTrigV0_Cat4.weights.xml");
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/Electrons_BDTG_NonTrigV0_Cat5.weights.xml");
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/Electrons_BDTG_NonTrigV0_Cat6.weights.xml");
   
   elecMVA->initialize("BDT", EGammaMvaEleEstimator::kNonTrig, true, elecMVACatWeights);

   string FlatTreeProducerLepMVAPath = CMSSW_BASE+"/src/IPHCFlatTree/FlatTreeProducer/data/lepMVA/";
   mu_reader_high_b   = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/mu_pteta_high_b_BDTG.weights.xml" ,  "mu");
   mu_reader_high_e   = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/mu_pteta_high_e_BDTG.weights.xml" ,  "mu");
   mu_reader_low_b    = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/mu_pteta_low_b_BDTG.weights.xml"  ,  "mu");
   mu_reader_low_e    = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/mu_pteta_low_e_BDTG.weights.xml"  ,  "mu");
   ele_reader_high_cb = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_high_cb_BDTG.weights.xml", "ele");
   ele_reader_high_fb = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_high_fb_BDTG.weights.xml", "ele");
   ele_reader_high_ec = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_high_ec_BDTG.weights.xml", "ele");
   ele_reader_low_cb  = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_low_cb_BDTG.weights.xml" , "ele");
   ele_reader_low_fb  = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_low_fb_BDTG.weights.xml" , "ele");
   ele_reader_low_ec  = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_low_ec_BDTG.weights.xml" , "ele");
   
   // ###
   // Restore stdout
   // ###
   cout.rdbuf(oldStdout);
   cout << "Stdout now restored." << endl;

   // ########################
   // #  Create output tree  #
   // ########################

   TFile& f = fs->file();
   f.SetCompressionAlgorithm(ROOT::kZLIB);
   f.SetCompressionLevel(9);
   ftree = new FlatTree(fs->make<TTree>("tree","tree"));
   
   // #########################
   // #  Read XML config file #
   // #########################

   std::string confFile = iConfig.getParameter<std::string>("confFile");
   ReadConfFile(confFile);
   int buffersize = iConfig.getParameter<int>("bufferSize");
   if (buffersize != 0)
     ftree->CreateBranches(buffersize);
   else
     ftree->CreateBranches();
   
   xmlconf.LoadFile("conf.xml");
   XMLElement* tElement = xmlconf.FirstChildElement("var");

   for( XMLElement* child=tElement;child!=0;child=child->NextSiblingElement() )
     {
	std::string vname = child->ToElement()->Attribute("name");
	std::string vsave = child->ToElement()->Attribute("save");
	bool bsave = atoi(vsave.c_str());
	
	ftree->conf.insert(std::make_pair(vname,bsave));
     }
   
   // ###############################
   // #  Add count & weight histos  #
   // ###############################
   
   hcount = fs->make<TH1D>("hcount","hcount",1,0.,1.);
   hweight = fs->make<TH1D>("hweight","hweight",1,0.,1.);

   // #############################################################
   // #  Read parameters from python file and get consume tokens  #
   // #############################################################

   dataFormat_        = iConfig.getParameter<std::string>("dataFormat");
   isData_            = iConfig.getParameter<bool>("isData");
   triggerBits_       = consumes<edm::TriggerResults>(edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("HLT")));
   triggerPrescales_  = consumes<pat::PackedTriggerPrescales>(edm::InputTag(std::string("patTrigger")));
   vertexToken_       = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexInput"));
   electronToken_     = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronInput"));
   muonToken_         = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonInput"));
   jetToken_          = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetInput"));
   metTokenAOD_       = consumes<std::vector<pat::MET> >(iConfig.getParameter<edm::InputTag>("metInput"));
   metTokenMINIAOD_   = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metInput"));
   rhoToken_          = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInput"));
   genParticlesToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticlesInput"));
}

FlatTreeProducer::~FlatTreeProducer()
{
}

void FlatTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   hcount->SetBinContent(1,hcount->GetBinContent(1)+1);

   ftree->Init();

   // Initial-state info
   edm::Handle<GenEventInfoProduct> genEventInfo;
   iEvent.getByLabel("generator",genEventInfo);

   // Gen particles
   edm::Handle<reco::GenParticleCollection> genParticlesHandle;                                                          
   iEvent.getByToken(genParticlesToken_,genParticlesHandle);

   // Beamspot
   edm::Handle<reco::BeamSpot> bsHandle;
   iEvent.getByLabel("offlineBeamSpot", bsHandle);
   const reco::BeamSpot &beamspot = *bsHandle.product();
 
   // Primary vertex
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vertexToken_,vertices);

   // Triggers
   edm::Handle<edm::TriggerResults> triggerBits;
   iEvent.getByToken(triggerBits_,triggerBits);
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   iEvent.getByToken(triggerPrescales_, triggerPrescales);

   // Pile-up
   edm::Handle<std::vector< PileupSummaryInfo> > pileupInfo;
   iEvent.getByLabel("addPileupInfo",pileupInfo);

   // Rho info
   edm::Handle<double> rhoPtr;
   iEvent.getByToken(rhoToken_,rhoPtr);

   // Packed candidate collection
   edm::Handle<pat::PackedCandidateCollection> pfcands;
   if( dataFormat_ != "AOD" ) iEvent.getByLabel("packedPFCandidates",pfcands);

   // Jets
   edm::Handle<pat::JetCollection> jets;
   iEvent.getByToken(jetToken_,jets);

   // Muons
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_,muons);
   
   // Electrons
   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_,electrons);
   
   // Conversions info
   edm::Handle<reco::ConversionCollection> hConversions;
   iEvent.getByLabel("reducedEgamma","reducedConversions",hConversions);
  
   // ###############################################################
   // #    ____                           _     _        __         #
   // #   / ___| ___ _ __   ___ _ __ __ _| |   (_)_ __  / _| ___    #
   // #  | |  _ / _ \ '_ \ / _ \ '__/ _` | |   | | '_ \| |_ / _ \   #
   // #  | |_| |  __/ | | |  __/ | | (_| | |   | | | | |  _| (_) |  #
   // #   \____|\___|_| |_|\___|_|  \__,_|_|   |_|_| |_|_|  \___/   #
   // #                                                             #
   // ###############################################################                           

   ftree->ev_run = iEvent.id().run();
   ftree->ev_id = iEvent.id().event();
   ftree->ev_lumi = iEvent.id().luminosityBlock();

   // ##########################################################
   // #   ___       _ _   _       _         _        _         #
   // #  |_ _|_ __ (_) |_(_) __ _| |    ___| |_ __ _| |_ ___   #
   // #   | || '_ \| | __| |/ _` | |   / __| __/ _` | __/ _ \  #
   // #   | || | | | | |_| | (_| | |   \__ \ || (_| | ||  __/  #
   // #  |___|_| |_|_|\__|_|\__,_|_|   |___/\__\__,_|\__\___|  #
   // #                                                        #
   // ##########################################################
   
   ftree->mc_weight = genEventInfo->weight();
   ftree->mc_id = genEventInfo->signalProcessID();
   ftree->mc_f1 = genEventInfo->pdf()->id.first;
   ftree->mc_f2 = genEventInfo->pdf()->id.second;
   ftree->mc_x1 = genEventInfo->pdf()->x.first;
   ftree->mc_x2 = genEventInfo->pdf()->x.second;
   ftree->mc_scale = genEventInfo->pdf()->scalePDF;
   if( genEventInfo->binningValues().size() > 0 ) ftree->mc_ptHat = genEventInfo->binningValues()[0];

   hweight->SetBinContent(1,hweight->GetBinContent(1)+ftree->mc_weight);

   // ####################################
   // #   ____  _ _                      #
   // #  |  _ \(_) | ___   _   _ _ __    #
   // #  | |_) | | |/ _ \ | | | | '_ \   #
   // #  |  __/| | |  __/ | |_| | |_) |  #
   // #  |_|   |_|_|\___|  \__,_| .__/   #
   // #                         |_|      #
   // #                                  #                      
   // ####################################

   ftree->mc_pu_Npvi = pileupInfo->size();
   unsigned int vecSize = std::distance(pileupInfo->begin(), pileupInfo->end());
   ftree->mc_pu_BunchCrossing.resize(vecSize);
   ftree->mc_pu_Nzpositions.resize(vecSize);
   ftree->mc_pu_zpositions.resize(vecSize);
   ftree->mc_pu_sumpT_lowpT.resize(vecSize);
   ftree->mc_pu_sumpT_highpT.resize(vecSize);
   ftree->mc_pu_ntrks_lowpT.resize(vecSize);
   ftree->mc_pu_ntrks_highpT.resize(vecSize);
   for(std::vector<PileupSummaryInfo>::const_iterator pvi=pileupInfo->begin();
       pvi!=pileupInfo->end();pvi++)
     {
       unsigned int distance = std::distance(pileupInfo->begin(), pvi);
       signed int n_bc = pvi->getBunchCrossing();
       //ftree->mc_pu_BunchCrossing.push_back(n_bc);
       ftree->mc_pu_BunchCrossing[distance] = n_bc;
       if( n_bc == 0 )
	 {
	   ftree->mc_pu_intime_NumInt = pvi->getPU_NumInteractions();
	   ftree->mc_pu_trueNumInt = pvi->getTrueNumInteractions();
	  }
       else if( n_bc == -1 ) ftree->mc_pu_before_npu = pvi->getPU_NumInteractions();
       else if( n_bc == +1 ) ftree->mc_pu_after_npu  = pvi->getPU_NumInteractions();

       std::vector<float> mc_pu_zpositions;
       std::vector<float> mc_pu_sumpT_lowpT;
       std::vector<float> mc_pu_sumpT_highpT;
       std::vector<int> mc_pu_ntrks_lowpT;
       std::vector<int> mc_pu_ntrks_highpT;

       //       ftree->mc_pu_Nzpositions.push_back(pvi->getPU_zpositions().size());
       ftree->mc_pu_Nzpositions[distance] = pvi->getPU_zpositions().size();
       unsigned int vecSize1 =  pvi->getPU_zpositions().size();
       mc_pu_zpositions.resize(vecSize1);
       mc_pu_sumpT_lowpT.resize(vecSize1);
       mc_pu_sumpT_highpT.resize(vecSize1);
       mc_pu_ntrks_lowpT.resize(vecSize1);
       mc_pu_ntrks_highpT.resize(vecSize1);
       for( unsigned int ipu=0;ipu<pvi->getPU_zpositions().size();ipu++ )
	 {
	   /*mc_pu_zpositions.push_back((pvi->getPU_zpositions())[ipu]);
	   mc_pu_sumpT_lowpT.push_back((pvi->getPU_sumpT_lowpT())[ipu]);
	   mc_pu_sumpT_highpT.push_back((pvi->getPU_sumpT_highpT())[ipu]);
	   mc_pu_ntrks_lowpT.push_back((pvi->getPU_ntrks_lowpT())[ipu]);
	   mc_pu_ntrks_highpT.push_back((pvi->getPU_ntrks_highpT())[ipu]);*/
	   mc_pu_zpositions[ipu] = (pvi->getPU_zpositions())[ipu];
	   mc_pu_sumpT_lowpT[ipu] = (pvi->getPU_sumpT_lowpT())[ipu];
	   mc_pu_sumpT_highpT[ipu] = (pvi->getPU_sumpT_highpT())[ipu];
	   mc_pu_ntrks_lowpT[ipu] = (pvi->getPU_ntrks_lowpT())[ipu];
	   mc_pu_ntrks_highpT[ipu] = (pvi->getPU_ntrks_highpT())[ipu];
	 }

       /*ftree->mc_pu_zpositions.push_back(mc_pu_zpositions);
       ftree->mc_pu_sumpT_lowpT.push_back(mc_pu_sumpT_lowpT);
       ftree->mc_pu_sumpT_highpT.push_back(mc_pu_sumpT_highpT);
       ftree->mc_pu_ntrks_lowpT.push_back(mc_pu_ntrks_lowpT);
       ftree->mc_pu_ntrks_highpT.push_back(mc_pu_ntrks_highpT);*/
       ftree->mc_pu_zpositions[distance] = mc_pu_zpositions;
       ftree->mc_pu_sumpT_lowpT[distance] = mc_pu_sumpT_lowpT;
       ftree->mc_pu_sumpT_highpT[distance] = mc_pu_sumpT_highpT;
       ftree->mc_pu_ntrks_lowpT[distance] = mc_pu_ntrks_lowpT;
       ftree->mc_pu_ntrks_highpT[distance] = mc_pu_ntrks_highpT;
     }

   // ##################################################
   // #   __  __  ____     _____           _   _       #
   // #  |  \/  |/ ___|   |_   _| __ _   _| |_| |__    #
   // #  | |\/| | |         | || '__| | | | __| '_ \   #
   // #  | |  | | |___      | || |  | |_| | |_| | | |  #
   // #  |_|  |_|\____|     |_||_|   \__,_|\__|_| |_|  #
   // #                                                #
   // ##################################################

   bool do_mc_truth_tth = ftree->doWrite("mc_truth_tth");
   bool do_mc_truth_tzq = ftree->doWrite("mc_truth_tzq");

   MCTruth *mc_truth = new MCTruth();

   bool reqMCTruth = 0;
   if( (do_mc_truth_tth || do_mc_truth_tzq ) &&
       !isData_ )
     {
	mc_truth->Init(*ftree);
	if( do_mc_truth_tth ) mc_truth->fillTTHSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
	if( do_mc_truth_tzq ) mc_truth->fillTZQSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
	reqMCTruth = 1;
     }

   bool do_gen_all = ftree->doWrite("gen_all");

   if( do_gen_all &&
       !isData_ )
     {
	if( !reqMCTruth ) mc_truth->Init(*ftree);
	mc_truth->fillGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
     }
  
   // #########################################
   // #   _____     _                         #
   // #  |_   _| __(_) __ _  __ _  ___ _ __   #
   // #    | || '__| |/ _` |/ _` |/ _ \ '__|  #
   // #    | || |  | | (_| | (_| |  __/ |     #
   // #    |_||_|  |_|\__, |\__, |\___|_|     #
   // #               |___/ |___/             #
   // #                                       #
   // #########################################

   std::vector<std::string> triggerIdentifiers_;
   triggerIdentifiers_.push_back("HLT_Ele27_WP80_v*");

   for (unsigned int j = 0; j < triggerIdentifiers_.size(); ++j)
     {
	std::string idName = triggerIdentifiers_[j];
	std::string idNameUnstarred = idName;
	bool isStarred = (idName.find("*")!=std::string::npos);
	if( isStarred ) idNameUnstarred.erase( idName.find("*"), 1 );

	for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
	  {
	     if( (isStarred && names.triggerName(i).find(idNameUnstarred)!=std::string::npos ) ||
		 (!isStarred && names.triggerName(i)==idName) )
		 {
//	     std::cout << "Trigger " << names.triggerName(i) <<
//	       ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
//	       ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") << std::endl;
		 }
	  }
     }
   
   reco::Vertex *primVtx = NULL;   

   if( ! vertices->empty() )
     {
	const reco::Vertex &PV = vertices->front();
	primVtx = (reco::Vertex*)&PV;
     }
   if( primVtx )
     {
	ftree->pv_x = primVtx->position().x();
	ftree->pv_y = primVtx->position().y();
	ftree->pv_z = primVtx->position().z();
     }

   // Rho
   ftree->ev_rho = *rhoPtr;

   // ####################################################
   // #   __  __ _         _               _____ _____   #
   // #  |  \/  (_)___ ___(_)_ __   __ _  | ____|_   _|  #
   // #  | |\/| | / __/ __| | '_ \ / _` | |  _|   | |    #
   // #  | |  | | \__ \__ \ | | | | (_| | | |___  | |    #
   // #  |_|  |_|_|___/___/_|_| |_|\__, | |_____| |_|    #
   // #                            |___/                 # 
   // #                                                  #
   // ####################################################

   // MET
   edm::Handle<std::vector<pat::MET> > metAOD;
   edm::Handle<pat::METCollection> metMINIAOD;
   if( dataFormat_ != "AOD" )
     {
	iEvent.getByToken(metTokenMINIAOD_,metMINIAOD);
	const pat::MET &metv = metMINIAOD->front();
	ftree->met_pt = metv.pt();
	ftree->met_phi = metv.phi();
	ftree->met_sumet = metv.sumEt();
     }
   else
     {
	iEvent.getByToken(metTokenAOD_,metAOD);
	const pat::MET &metv = metAOD->front();
	ftree->met_pt = metv.pt();
	ftree->met_phi = metv.phi();
	ftree->met_sumet = metv.sumEt();
     }
 
   // #################################################
   // #   _____ _           _                         #
   // #  | ____| | ___  ___| |_ _ __ ___  _ __  ___   #
   // #  |  _| | |/ _ \/ __| __| '__/ _ \| '_ \/ __|  #
   // #  | |___| |  __/ (__| |_| | | (_) | | | \__ \  #
   // #  |_____|_|\___|\___|\__|_|  \___/|_| |_|___/  #
   // #                                               #
   // #################################################

   int nElec = electrons->size();
   ftree->el_n = nElec;

   ftree->el_pt.resize(nElec);
   ftree->el_eta.resize(nElec);
   ftree->el_phi.resize(nElec);
   ftree->el_m.resize(nElec);
   ftree->el_E.resize(nElec);
   ftree->el_id.resize(nElec);
   ftree->el_charge.resize(nElec);
   ftree->el_scleta.resize(nElec);
   ftree->el_isGsfCtfScPixChargeConsistent.resize(nElec);
   ftree->el_sigmaIetaIeta.resize(nElec);
   ftree->el_sigmaIphiIphi.resize(nElec);
   ftree->el_hadronicOverEm.resize(nElec);
   ftree->el_dr03TkSumPt.resize(nElec);
   ftree->el_dr03EcalRecHitSumEt.resize(nElec);
   ftree->el_dr03HcalTowerSumEt.resize(nElec);
   ftree->el_numberOfLostHits.resize(nElec);
   ftree->el_fbrem.resize(nElec);
   ftree->el_deltaEtaSuperClusterTrackAtVtx.resize(nElec);
   ftree->el_deltaPhiSuperClusterTrackAtVtx.resize(nElec);
   ftree->el_deltaEtaSeedClusterTrackAtCalo.resize(nElec);
   ftree->el_see.resize(nElec);
   ftree->el_spp.resize(nElec);
   ftree->el_superClusterEtaWidth.resize(nElec);
   ftree->el_superClusterPhiWidth.resize(nElec);
   ftree->el_full5x5_OneMinusE1x5E5x5.resize(nElec);
   ftree->el_OneMinusE1x5E5x5.resize(nElec);
   ftree->el_full5x5_r9.resize(nElec);
   ftree->el_r9.resize(nElec);
   ftree->el_eSuperClusterOverP.resize(nElec);
   ftree->el_IoEmIoP.resize(nElec);
   ftree->el_eleEoPout.resize(nElec);
   ftree->el_PreShowerOverRaw.resize(nElec);
   ftree->el_dB3D.resize(nElec);
   ftree->el_edB3D.resize(nElec);
   ftree->el_mvaNonTrigV0.resize(nElec);
   ftree->el_neutralHadronIso.resize(nElec);
   ftree->el_chargedHadronIso.resize(nElec);
   ftree->el_puChargedHadronIso.resize(nElec);
   ftree->el_ecalIso.resize(nElec);
   ftree->el_hcalIso.resize(nElec);
   ftree->el_particleIso.resize(nElec);
   ftree->el_photonIso.resize(nElec);
   ftree->el_trackIso.resize(nElec);
   ftree->el_pfIso_sumChargedHadronPt.resize(nElec);
   ftree->el_pfIso_sumNeutralHadronEt.resize(nElec);
   ftree->el_pfIso_sumPhotonEt.resize(nElec);
   ftree->el_pfIso_sumPUPt.resize(nElec);
   ftree->el_isLoose.resize(nElec);
   ftree->el_isTight.resize(nElec);
   ftree->el_isRobustLoose.resize(nElec);
   ftree->el_isRobustTight.resize(nElec);
   ftree->el_isRobustHighEnergy.resize(nElec);
   ftree->el_vx.resize(nElec);
   ftree->el_vy.resize(nElec);
   ftree->el_vz.resize(nElec);
   ftree->el_isGsf.resize(nElec);
   ftree->el_passConversionVeto.resize(nElec);
   ftree->el_numberOfHits.resize(nElec);
   ftree->el_dxy.resize(nElec);
   ftree->el_dz.resize(nElec);
   ftree->el_dxyError.resize(nElec);
   ftree->el_dzError.resize(nElec);
   ftree->el_lepMVA.resize(nElec);
   ftree->el_lepMVA_neuRelIso.resize(nElec);
   ftree->el_lepMVA_chRelIso.resize(nElec);
   ftree->el_lepMVA_jetDR.resize(nElec);
   ftree->el_lepMVA_jetPtRatio.resize(nElec);
   ftree->el_lepMVA_jetBTagCSV.resize(nElec);
   ftree->el_lepMVA_sip3d.resize(nElec);
   ftree->el_lepMVA_mvaId.resize(nElec);
   ftree->el_lepMVA_innerHits.resize(nElec);

   for(int ie=0;ie<nElec;ie++)
     {
	const pat::Electron& elec = electrons->at(ie);

	// Skimming electrons with pT < 5 GeV.
	if (elec.pt() < 5) continue;
	
	ftree->el_pt[ie] = elec.pt();
	ftree->el_eta[ie] = elec.eta();
	ftree->el_phi[ie] = elec.phi();
	ftree->el_m[ie] = elec.mass();
	ftree->el_E[ie] = elec.energy();
	ftree->el_id[ie] = elec.pdgId();
	ftree->el_charge[ie] = elec.charge();
	
	ftree->el_scleta[ie] = elec.superCluster()->eta();
	ftree->el_isGsfCtfScPixChargeConsistent[ie] = elec.isGsfCtfScPixChargeConsistent();
	ftree->el_sigmaIetaIeta[ie] = elec.sigmaIetaIeta();
	ftree->el_sigmaIphiIphi[ie] = elec.sigmaIphiIphi();
	ftree->el_hadronicOverEm[ie] = elec.hadronicOverEm();
	ftree->el_dr03TkSumPt[ie] = elec.dr03TkSumPt();
	ftree->el_dr03EcalRecHitSumEt[ie] = elec.dr03EcalRecHitSumEt();
	ftree->el_dr03HcalTowerSumEt[ie] = elec.dr03HcalTowerSumEt();

	int numberOfLostHits = -666;

	if( elec.gsfTrack().isNonnull() )
	  {
	     numberOfLostHits = elec.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::HitCategory::MISSING_INNER_HITS);
	  }

	ftree->el_numberOfLostHits[ie] = numberOfLostHits;

	bool validKF = false;
	reco::TrackRef myTrackRef = elec.closestCtfTrackRef();
	validKF = (myTrackRef.isAvailable());
	validKF = (myTrackRef.isNonnull());

	ftree->el_fbrem[ie] = elec.fbrem();
	float el_kf_normalizedChi2 = ( validKF ) ? myTrackRef->normalizedChi2() : 666.;
	if( validKF ) ftree->el_kf_normalizedChi2.push_back(myTrackRef->normalizedChi2());
	float el_trackerLayersWithMeasurement = ( validKF ) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -666.;
	if( validKF ) ftree->el_trackerLayersWithMeasurement.push_back(myTrackRef->hitPattern().trackerLayersWithMeasurement());
	if( elec.gsfTrack().isNonnull() ) ftree->el_gsf_normalizedChi2.push_back(elec.gsfTrack()->normalizedChi2());
	float el_gsf_normalizedChi2 = ( elec.gsfTrack().isNonnull() ) ? elec.gsfTrack()->normalizedChi2() : 666.;

	ftree->el_deltaEtaSuperClusterTrackAtVtx[ie] = elec.deltaEtaSuperClusterTrackAtVtx();
	ftree->el_deltaPhiSuperClusterTrackAtVtx[ie] = elec.deltaPhiSuperClusterTrackAtVtx();
	ftree->el_deltaEtaSeedClusterTrackAtCalo[ie] = elec.deltaEtaSeedClusterTrackAtCalo();

	ftree->el_see[ie] = elec.full5x5_sigmaIetaIeta();
	ftree->el_spp[ie] = elec.full5x5_sigmaIphiIphi();

	ftree->el_superClusterEtaWidth[ie] = elec.superCluster()->etaWidth();
	ftree->el_superClusterPhiWidth[ie] = elec.superCluster()->phiWidth();

	double OneMinusE1x5E5x5 = (elec.e5x5() != 0.) ? 1.-(elec.e1x5()/elec.e5x5()) : -1.;
	double full5x5_OneMinusE1x5E5x5 = (elec.full5x5_e5x5() != 0.) ? 1.-(elec.full5x5_e1x5()/elec.full5x5_e5x5()) : -1.;

	ftree->el_full5x5_OneMinusE1x5E5x5[ie] = full5x5_OneMinusE1x5E5x5;
	ftree->el_OneMinusE1x5E5x5[ie] = OneMinusE1x5E5x5;

	ftree->el_full5x5_r9[ie] = elec.full5x5_r9();

	ftree->el_r9[ie] = elec.r9();
	ftree->el_eSuperClusterOverP[ie] = elec.eSuperClusterOverP();

	double IoEmIoP = (1.0/elec.ecalEnergy())-(1.0/elec.p());
	ftree->el_IoEmIoP[ie] = IoEmIoP;
	ftree->el_eleEoPout[ie] = elec.eEleClusterOverPout();
	double PreShowerOverRaw = elec.superCluster()->preshowerEnergy()/elec.superCluster()->rawEnergy();
	ftree->el_PreShowerOverRaw[ie] = PreShowerOverRaw;

	// https://github.com/gpetruc/cmg-cmssw/blob/CMG_MiniAOD_Lite_V6_0_from-CMSSW_7_0_6/EgammaAnalysis/ElectronTools/src/EGammaMvaEleEstimator.cc#L1244-1336
	ftree->el_dB3D[ie] = elec.dB(pat::Electron::PV3D);
	ftree->el_edB3D[ie] = elec.edB(pat::Electron::PV3D);

	double mvaNonTrigV0 = elecMVA->mvaValue(ftree->el_fbrem[ie],
						el_kf_normalizedChi2,
						el_trackerLayersWithMeasurement,
						el_gsf_normalizedChi2,
						ftree->el_deltaEtaSuperClusterTrackAtVtx[ie],
						ftree->el_deltaPhiSuperClusterTrackAtVtx[ie],
						ftree->el_deltaEtaSeedClusterTrackAtCalo[ie],
						ftree->el_see[ie],
						ftree->el_spp[ie],
						ftree->el_superClusterEtaWidth[ie],
						ftree->el_superClusterPhiWidth[ie],
						ftree->el_full5x5_OneMinusE1x5E5x5[ie],
						ftree->el_full5x5_r9[ie],
						ftree->el_hadronicOverEm[ie],
						ftree->el_eSuperClusterOverP[ie],
						ftree->el_IoEmIoP[ie],
						ftree->el_eleEoPout[ie],
//						   ftree->ev_rho,
						ftree->el_PreShowerOverRaw[ie],
						ftree->el_scleta[ie],
						ftree->el_pt[ie],
						false);

	ftree->el_mvaNonTrigV0[ie] = mvaNonTrigV0;

	ftree->el_neutralHadronIso[ie] = elec.neutralHadronIso();
	ftree->el_chargedHadronIso[ie] = elec.chargedHadronIso();
	ftree->el_puChargedHadronIso[ie] = elec.puChargedHadronIso();
	ftree->el_ecalIso[ie] = elec.ecalIso();
	ftree->el_hcalIso[ie] = elec.hcalIso();
	ftree->el_particleIso[ie] = elec.particleIso();
	ftree->el_photonIso[ie] = elec.photonIso();
	ftree->el_trackIso[ie] = elec.trackIso();

	ftree->el_pfIso_sumChargedHadronPt[ie] = elec.pfIsolationVariables().sumChargedHadronPt;
	ftree->el_pfIso_sumNeutralHadronEt[ie] = elec.pfIsolationVariables().sumNeutralHadronEt;
	ftree->el_pfIso_sumPhotonEt[ie] = elec.pfIsolationVariables().sumPhotonEt;
	ftree->el_pfIso_sumPUPt[ie] = elec.pfIsolationVariables().sumPUPt;

	// mini-iso
	float miniIso = -666;
	if( dataFormat_ != "AOD" )
	  {
	     miniIso = getPFIsolation(pfcands,dynamic_cast<const reco::Candidate*>(&elec),0.05,0.2,10.,false,false);
	  }
	ftree->el_miniIso[ie] = miniIso;
	       
	ftree->el_isLoose[ie] = elec.electronID("eidLoose");
	ftree->el_isTight[ie] = elec.electronID("eidTight");
	ftree->el_isRobustLoose[ie] = elec.electronID("eidRobustLoose");
	ftree->el_isRobustTight[ie] = elec.electronID("eidRobustTight");
	ftree->el_isRobustHighEnergy[ie] = elec.electronID("eidRobustHighEnergy");

	ftree->el_vx[ie] = elec.vx();
	ftree->el_vy[ie] = elec.vy();
	ftree->el_vz[ie] = elec.vz();

	ftree->el_isGsf[ie] = elec.gsfTrack().isNonnull();

	ftree->el_passConversionVeto[ie] = elec.passConversionVeto();

	int numberOfHits = -666;

	float dxy = -10E+10;
	float dz = -10E+10;
	float dxyError = -666.;
	float dzError = -666.;

	if( elec.gsfTrack().isNonnull() )
	  {
	     numberOfHits = elec.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::HitCategory::MISSING_INNER_HITS);

	     if( primVtx ) dxy = elec.gsfTrack()->dxy(primVtx->position());
	     if( primVtx ) dz = elec.gsfTrack()->dz(primVtx->position());
	     dxyError = elec.gsfTrack()->dxyError();
	     dzError = elec.gsfTrack()->dzError();
	  }

	ftree->el_numberOfHits[ie] = numberOfHits;
	ftree->el_dxy[ie] = dxy;
	ftree->el_dz[ie] = dz;
	ftree->el_dxyError[ie] = dxyError;
	ftree->el_dzError[ie] = dzError;

	double el_pt = elec.pt();
	double el_eta = elec.eta();
	double el_scleta = ftree->el_scleta[ie];
	double el_phi = elec.phi();
	double el_lepMVA = -666.;

	double relIso = (ftree->el_pfIso_sumChargedHadronPt[ie] +
			 std::max(ftree->el_pfIso_sumNeutralHadronEt[ie] +
				  ftree->el_pfIso_sumPhotonEt[ie] -
				  ftree->el_pfIso_sumPUPt[ie]/2.,0.0))/el_pt;

	lepMVA_neuRelIso = relIso - ftree->el_chargedHadronIso[ie]/el_pt;
	lepMVA_chRelIso = ftree->el_chargedHadronIso[ie]/el_pt;
	float drmin = 0.5;
	int jcl = -1;
	for(unsigned int ij=0;ij<jets->size();ij++)
	  {
	     if( jets->at(ij).pt() < 10. ) continue;
	     float dr = GetDeltaR(jets->at(ij).eta(),
				  jets->at(ij).phi(),
				  el_eta,
				  el_phi);
	     if( dr < drmin )
	       {
		  drmin = dr;
		  jcl = ij;
	       }
	  }
	lepMVA_jetDR = drmin;
	lepMVA_jetPtRatio = (jcl >= 0) ? std::min(el_pt/jets->at(jcl).pt(),1.5) : 1.5;
	float csv = (jcl >= 0) ? jets->at(jcl).bDiscriminator("combinedSecondaryVertexBJetTags") : 0.;
	if( csv == -1000. ) csv = jets->at(jcl).bDiscriminator("pfCombinedSecondaryVertexBJetTags");
	lepMVA_jetBTagCSV = std::max(double(csv),0.);
	lepMVA_sip3d = fabs(ftree->el_dB3D[ie]/ftree->el_edB3D[ie]);
	lepMVA_mvaId = mvaNonTrigV0;
	
	if( elec.gsfTrack().isNonnull() ) lepMVA_innerHits = ftree->el_numberOfHits[ie];

	if( el_pt >= 10 && fabs(el_scleta) <= 0.8 ) el_lepMVA = ele_reader_high_cb->EvaluateMVA("BDTG method");
	else if( el_pt >= 10 && fabs(el_scleta) > 0.8 && fabs(el_scleta) <= 1.479 ) el_lepMVA = ele_reader_high_fb->EvaluateMVA("BDTG method");
	else if( el_pt >= 10 && fabs(el_scleta) > 1.479 )                           el_lepMVA = ele_reader_high_ec->EvaluateMVA("BDTG method");
	else if( el_pt < 10  && fabs(el_scleta) <= 0.8 )                            el_lepMVA = ele_reader_low_cb ->EvaluateMVA("BDTG method");
	else if( el_pt < 10  && fabs(el_scleta) > 0.8 && fabs(el_scleta) <= 1.479 ) el_lepMVA = ele_reader_low_fb ->EvaluateMVA("BDTG method");
	else                                                                        el_lepMVA = ele_reader_low_ec ->EvaluateMVA("BDTG method");

        ftree->el_lepMVA[ie] = el_lepMVA;

	if( !isData_ )
	  {
	     reco::GenParticle *genp = new reco::GenParticle();

	     float drmin;
	     bool hasMCMatch = mc_truth->doMatch(iEvent,iSetup,genParticlesHandle,genp,drmin,
						 elec.pt(),elec.eta(),elec.phi(),elec.pdgId());
	     ftree->el_hasMCMatch.push_back(hasMCMatch);
	     if( hasMCMatch )
	       {
		  ftree->el_gen_pt.push_back(genp->pt());
		  ftree->el_gen_eta.push_back(genp->eta());
		  ftree->el_gen_phi.push_back(genp->phi());
		  ftree->el_gen_m.push_back(genp->mass());
		  ftree->el_gen_status.push_back(genp->status());
		  ftree->el_gen_id.push_back(genp->pdgId());
		  ftree->el_gen_charge.push_back(genp->charge());
		  ftree->el_gen_dr.push_back(drmin);
	       }

	     delete genp;	     
	  }

	ftree->el_lepMVA_neuRelIso[ie] = lepMVA_neuRelIso;
	ftree->el_lepMVA_chRelIso[ie] = lepMVA_chRelIso;
	ftree->el_lepMVA_jetDR[ie] = lepMVA_jetDR;
	ftree->el_lepMVA_jetPtRatio[ie] = lepMVA_jetPtRatio;
	ftree->el_lepMVA_jetBTagCSV[ie] = lepMVA_jetBTagCSV;
	ftree->el_lepMVA_sip3d[ie] = lepMVA_sip3d;
	ftree->el_lepMVA_mvaId[ie] = lepMVA_mvaId;
	ftree->el_lepMVA_innerHits[ie] = lepMVA_innerHits;

	bool allowCkfMatch = true;
	float lxyMin = 2.0;
	float probMin = 1e-6;
	uint nHitsBeforeVtxMax = 0;

	bool matchConv = 0;
	if( &beamspot ) matchConv = ConversionTools::hasMatchedConversion(elec,hConversions,beamspot.position(),allowCkfMatch,lxyMin,probMin,nHitsBeforeVtxMax);
	ftree->el_hasMatchedConversion.push_back(matchConv);
     }   
   
   // ####################################
   // #   __  __                         #
   // #  |  \/  |_   _  ___  _ __  ___   #
   // #  | |\/| | | | |/ _ \| '_ \/ __|  #
   // #  | |  | | |_| | (_) | | | \__ \  #
   // #  |_|  |_|\__,_|\___/|_| |_|___/  #
   // #                                  #                                     
   // ####################################

   // Muons

   int nMuon = muons->size();
   ftree->mu_n = nMuon;
   ftree->mu_pt.resize(nMuon);
   ftree->mu_eta.resize(nMuon);
   ftree->mu_phi.resize(nMuon);
   ftree->mu_m.resize(nMuon);
   ftree->mu_E.resize(nMuon);
   ftree->mu_id.resize(nMuon);
   ftree->mu_charge.resize(nMuon);
   ftree->mu_dB3D.resize(nMuon);
   ftree->mu_edB3D.resize(nMuon);
   ftree->mu_neutralHadronIso.resize(nMuon);
   ftree->mu_chargedHadronIso.resize(nMuon);
   ftree->mu_puChargedHadronIso.resize(nMuon);
   ftree->mu_ecalIso.resize(nMuon);
   ftree->mu_hcalIso.resize(nMuon);
   ftree->mu_photonIso.resize(nMuon);
   ftree->mu_trackIso.resize(nMuon);
   ftree->mu_pfIso03_sumChargedHadronPt.resize(nMuon);
   ftree->mu_pfIso03_sumNeutralHadronEt.resize(nMuon);
   ftree->mu_pfIso03_sumPhotonEt.resize(nMuon);
   ftree->mu_pfIso03_sumPUPt.resize(nMuon);
   ftree->mu_isGlobalMuon.resize(nMuon);
   ftree->mu_isTrackerMuon.resize(nMuon);
   ftree->mu_isStandAloneMuon.resize(nMuon);
   ftree->mu_isCaloMuon.resize(nMuon);
   ftree->mu_isPFMuon.resize(nMuon);
   ftree->mu_vx.resize(nMuon);
   ftree->mu_vy.resize(nMuon);
   ftree->mu_vz.resize(nMuon);
   ftree->mu_hasGlobalTrack.resize(nMuon);
   ftree->mu_hasInnerTrack.resize(nMuon);
   ftree->mu_globalTrack_dxy.resize(nMuon);
   ftree->mu_globalTrack_dz.resize(nMuon);
   ftree->mu_globalTrack_dxyError.resize(nMuon);
   ftree->mu_globalTrack_dzError.resize(nMuon);
   ftree->mu_innerTrack_dxy.resize(nMuon);
   ftree->mu_innerTrack_dz.resize(nMuon);
   ftree->mu_innerTrack_dxyError.resize(nMuon);
   ftree->mu_innerTrack_dzError.resize(nMuon);
   ftree->mu_bestTrack_dxy.resize(nMuon);
   ftree->mu_bestTrack_dz.resize(nMuon);
   ftree->mu_bestTrack_dxyError.resize(nMuon);
   ftree->mu_bestTrack_dzError.resize(nMuon);
   ftree->mu_innerTrack_pt.resize(nMuon);
   ftree->mu_innerTrack_ptError.resize(nMuon);
   ftree->mu_numberOfMatches.resize(nMuon);
   ftree->mu_numberOfValidMuonHits.resize(nMuon);
   ftree->mu_lepMVA.resize(nMuon);
   ftree->mu_lepMVA_neuRelIso.resize(nMuon);
   ftree->mu_lepMVA_chRelIso.resize(nMuon);
   ftree->mu_lepMVA_jetDR.resize(nMuon);
   ftree->mu_lepMVA_jetPtRatio.resize(nMuon);
   ftree->mu_lepMVA_jetBTagCSV.resize(nMuon);
   ftree->mu_lepMVA_sip3d.resize(nMuon);
   ftree->mu_lepMVA_dxy.resize(nMuon);
   ftree->mu_lepMVA_dz.resize(nMuon);

   for(int im=0;im<nMuon;im++)
     {
	const pat::Muon& muon = muons->at(im);
	
	// Skimming electrons with pT < 5 GeV.
	if (muon.pt() < 5) continue;

	ftree->mu_pt[im] = muon.pt();
	ftree->mu_eta[im] = muon.eta();
	ftree->mu_phi[im] = muon.phi();
	ftree->mu_m[im] = muon.mass();
	ftree->mu_E[im] = muon.energy();
	ftree->mu_id[im] = muon.pdgId();
	ftree->mu_charge[im] = muon.charge();

	ftree->mu_dB3D[im] = muon.dB(pat::Muon::PV3D);
	ftree->mu_edB3D[im] = muon.edB(pat::Muon::PV3D);

	ftree->mu_neutralHadronIso[im] = muon.neutralHadronIso();
	ftree->mu_chargedHadronIso[im] = muon.chargedHadronIso();
	ftree->mu_puChargedHadronIso[im] = muon.puChargedHadronIso();
	ftree->mu_ecalIso[im] = muon.ecalIso();
	ftree->mu_hcalIso[im] = muon.hcalIso();
	ftree->mu_photonIso[im] = muon.photonIso();
	ftree->mu_trackIso[im] = muon.trackIso();

	ftree->mu_pfIso03_sumChargedHadronPt[im] = muon.pfIsolationR03().sumChargedHadronPt;
	ftree->mu_pfIso03_sumNeutralHadronEt[im] = muon.pfIsolationR03().sumNeutralHadronEt;
	ftree->mu_pfIso03_sumPhotonEt[im] = muon.pfIsolationR03().sumPhotonEt;
	ftree->mu_pfIso03_sumPUPt[im] = muon.pfIsolationR03().sumPUPt;

	// mini-iso
	float miniIso = -666;
	if( dataFormat_ != "AOD" )
	{
	  miniIso = getPFIsolation(pfcands,dynamic_cast<const reco::Candidate*>(&muon),0.05,0.2,10.,false,false);
	}
	ftree->mu_miniIso[im] = miniIso;

	ftree->mu_isGlobalMuon[im] = muon.isGlobalMuon();
	ftree->mu_isTrackerMuon[im] = muon.isTrackerMuon();
	ftree->mu_isStandAloneMuon[im] = muon.isStandAloneMuon();
	ftree->mu_isCaloMuon[im] = muon.isCaloMuon();
	ftree->mu_isPFMuon[im] = muon.isPFMuon();

	ftree->mu_vx[im] = muon.vx();
	ftree->mu_vy[im] = muon.vy();
	ftree->mu_vz[im] = muon.vz();

	ftree->mu_hasGlobalTrack[im] = muon.globalTrack().isNonnull();
	ftree->mu_hasInnerTrack[im] = muon.innerTrack().isNonnull();

	float globalTrack_dxy = -10E+10;
	float globalTrack_dz = -10E+10;
	float globalTrack_dxyError = -666.;
	float globalTrack_dzError = -666.;

	float innerTrack_dxy = -10E+10;
	float innerTrack_dz = -10E+10;
	float innerTrack_dxyError = -666.;
	float innerTrack_dzError = -666.;

	float bestTrack_dxy = -10E+10;
	float bestTrack_dz = -10E+10;
	float bestTrack_dxyError = -666.;
	float bestTrack_dzError = -666.;

	bool isTightMuon = 0;
	if( primVtx ) isTightMuon = muon.isTightMuon(*primVtx);
	ftree->mu_isTightMuon[im] = isTightMuon;

	if( muon.globalTrack().isNonnull() )
	  {
	     globalTrack_dxy = muon.globalTrack()->dxy(muon.globalTrack()->referencePoint());
	     globalTrack_dz = muon.globalTrack()->dz(muon.globalTrack()->referencePoint());
	     globalTrack_dxyError = muon.globalTrack()->dxyError();
	     globalTrack_dzError = muon.globalTrack()->dzError();
	  }
	if( muon.innerTrack().isNonnull() )
	  {
	     if( primVtx ) innerTrack_dxy = muon.innerTrack()->dxy(primVtx->position());
	     if( primVtx ) innerTrack_dz = muon.innerTrack()->dz(primVtx->position());
	     innerTrack_dxyError = muon.innerTrack()->dxyError();
	     innerTrack_dzError = muon.innerTrack()->dzError();
	  }

	bestTrack_dxy = muon.bestTrack()->dxy(muon.bestTrack()->referencePoint());
	bestTrack_dz = muon.bestTrack()->dz(muon.bestTrack()->referencePoint());
	bestTrack_dxyError = muon.bestTrack()->dxyError();
	bestTrack_dzError = muon.bestTrack()->dzError();

	ftree->mu_globalTrack_dxy[im] = globalTrack_dxy;
	ftree->mu_globalTrack_dz[im] = globalTrack_dz;
	ftree->mu_globalTrack_dxyError[im] = globalTrack_dxyError;
	ftree->mu_globalTrack_dzError[im] = globalTrack_dzError;

	ftree->mu_innerTrack_dxy[im] = innerTrack_dxy;
	ftree->mu_innerTrack_dz[im] = innerTrack_dz;
	ftree->mu_innerTrack_dxyError[im] = innerTrack_dxyError;
	ftree->mu_innerTrack_dzError[im] = innerTrack_dzError;

	ftree->mu_bestTrack_dxy[im] = bestTrack_dxy;
	ftree->mu_bestTrack_dz[im] = bestTrack_dz;
	ftree->mu_bestTrack_dxyError[im] = bestTrack_dxyError;
	ftree->mu_bestTrack_dzError[im] = bestTrack_dzError;

	float innerTrack_pt = -666.;
	float innerTrack_ptError = -666.;

	if( muon.innerTrack().isNonnull() )
	  {
	     innerTrack_pt = muon.innerTrack()->pt();
	     innerTrack_ptError = muon.innerTrack()->ptError();
	  }

	ftree->mu_innerTrack_pt[im] = innerTrack_pt;
	ftree->mu_innerTrack_ptError[im] = innerTrack_ptError;

	ftree->mu_numberOfMatches[im] = muon.numberOfMatches();

	int numberOfValidMuonHits = -666;

	if( muon.globalTrack().isNonnull() )
	  {
	     numberOfValidMuonHits = muon.globalTrack()->hitPattern().numberOfValidMuonHits();
	  }

	ftree->mu_numberOfValidMuonHits[im] = numberOfValidMuonHits;

	double mu_pt = muon.pt();
	double mu_eta = muon.eta();
	double mu_phi = muon.phi();
	double mu_lepMVA = -666.;

	double relIso = (muon.pfIsolationR04().sumChargedHadronPt +
			 std::max(muon.pfIsolationR04().sumNeutralHadronEt +
				  muon.pfIsolationR04().sumPhotonEt -
				  muon.pfIsolationR04().sumPUPt/2.,0.0))/mu_pt;

	lepMVA_chRelIso = ftree->mu_chargedHadronIso[im]/mu_pt;
	lepMVA_neuRelIso = relIso - lepMVA_chRelIso;
	float drmin = 0.5;
	int jcl = -1;
	for(unsigned int ij=0;ij<jets->size();ij++)
	  {
	     if( jets->at(ij).pt() < 10. ) continue;
	     float dr = GetDeltaR(jets->at(ij).eta(),
				  jets->at(ij).phi(),
				  mu_eta,
				  mu_phi);
	     if( dr < drmin )
	       {
		  drmin = dr;
		  jcl = ij;
	       }
	  }
	lepMVA_jetDR = drmin;
	lepMVA_jetPtRatio = (jcl >= 0) ? std::min(mu_pt/jets->at(jcl).pt(),1.5) : 1.5;
	float csv = (jcl >= 0) ? jets->at(jcl).bDiscriminator("combinedSecondaryVertexBJetTags") : 0.;
	if( csv == -1000. ) csv = jets->at(jcl).bDiscriminator("pfCombinedSecondaryVertexBJetTags");
	lepMVA_jetBTagCSV = std::max(double(csv),0.);
	lepMVA_sip3d = fabs(ftree->mu_dB3D[im]/ftree->mu_edB3D[im]);
	lepMVA_dxy = log(fabs(ftree->mu_innerTrack_dxy[im]));
	lepMVA_dz = log(fabs(ftree->mu_innerTrack_dz[im]));

	if( mu_pt > 15 && fabs(mu_eta) < 1.5 ) mu_lepMVA = mu_reader_high_b->EvaluateMVA("BDTG method");
	else if( mu_pt > 15 && fabs(mu_eta) >= 1.5 ) mu_lepMVA = mu_reader_high_e->EvaluateMVA("BDTG method");
	else if( mu_pt <= 15 && fabs(mu_eta) < 1.5 ) mu_lepMVA = mu_reader_low_b->EvaluateMVA("BDTG method");
	else mu_lepMVA = mu_reader_low_e->EvaluateMVA("BDTG method");

	ftree->mu_lepMVA[im] = mu_lepMVA;
	
	ftree->mu_lepMVA_neuRelIso[im] = lepMVA_neuRelIso;
	ftree->mu_lepMVA_chRelIso[im] = lepMVA_chRelIso;
	ftree->mu_lepMVA_jetDR[im] = lepMVA_jetDR;
	ftree->mu_lepMVA_jetPtRatio[im] = lepMVA_jetPtRatio;
	ftree->mu_lepMVA_jetBTagCSV[im] = lepMVA_jetBTagCSV;
	ftree->mu_lepMVA_sip3d[im] = lepMVA_sip3d;
	ftree->mu_lepMVA_dxy[im] = lepMVA_dxy;
	ftree->mu_lepMVA_dz[im] = lepMVA_dz;

	if( !isData_ )
	  {
	     reco::GenParticle *genp = new reco::GenParticle();

	     float drmin;
	     bool hasMCMatch = mc_truth->doMatch(iEvent,iSetup,genParticlesHandle,genp,drmin,
						 muon.pt(),muon.eta(),muon.phi(),muon.pdgId());
	     ftree->mu_hasMCMatch.push_back(hasMCMatch);
	     if( hasMCMatch )
	       {
		  ftree->mu_gen_pt.push_back(genp->pt());
		  ftree->mu_gen_eta.push_back(genp->eta());
		  ftree->mu_gen_phi.push_back(genp->phi());
		  ftree->mu_gen_m.push_back(genp->mass());
		  ftree->mu_gen_status.push_back(genp->status());
		  ftree->mu_gen_id.push_back(genp->pdgId());
		  ftree->mu_gen_charge.push_back(genp->charge());
		  ftree->mu_gen_dr.push_back(drmin);
	       }

	     delete genp;	     
	  }
     }   
 
   // ##########################
   // #       _      _         #
   // #      | | ___| |_ ___   #
   // #   _  | |/ _ \ __/ __|  #
   // #  | |_| |  __/ |_\__ \  #
   // #   \___/ \___|\__|___/  #
   // #                        #                        
   // ##########################

   // Jets

   int nJet = jets->size();
   ftree->jet_n = nJet;

   ftree->jet_pt.resize(nJet);
   ftree->jet_eta.resize(nJet);
   ftree->jet_phi.resize(nJet);
   ftree->jet_m.resize(nJet);
   ftree->jet_E.resize(nJet);
   ftree->jet_JBP.resize(nJet);
   ftree->jet_JP.resize(nJet);
   ftree->jet_TCHP.resize(nJet);
   ftree->jet_TCHE.resize(nJet);
   ftree->jet_SSVHE.resize(nJet);
   ftree->jet_SSVHP.resize(nJet);
   ftree->jet_CMVA.resize(nJet);
   ftree->jet_CSVv2.resize(nJet);
   ftree->jet_CSV.resize(nJet);
   ftree->jet_flavour.resize(nJet);
   ftree->jet_neutralHadronEnergy.resize(nJet);
   ftree->jet_neutralEmEnergy.resize(nJet);
   ftree->jet_chargedHadronEnergy.resize(nJet);
   ftree->jet_chargedEmEnergy.resize(nJet);
   ftree->jet_electronEnergy.resize(nJet);
   ftree->jet_muonEnergy.resize(nJet);
   ftree->jet_photonEnergy.resize(nJet);
   ftree->jet_pileupJetId.resize(nJet);

   for(int ij=0;ij<nJet;ij++)
     {
	const pat::Jet& jet = jets->at(ij);

	ftree->jet_pt[ij] = jet.pt();
	ftree->jet_eta[ij] = jet.eta();
	ftree->jet_phi[ij] = jet.phi();
	ftree->jet_m[ij] = jet.mass();
	ftree->jet_E[ij] = jet.energy();
	ftree->jet_JBP[ij] = jet.bDiscriminator("jetBProbabilityBJetTags");
	ftree->jet_JP[ij] = jet.bDiscriminator("jetProbabilityBJetTags");
	ftree->jet_TCHP[ij] = jet.bDiscriminator("trackCountingHighPurBJetTags");
	ftree->jet_TCHE[ij] = jet.bDiscriminator("trackCountingHighEffBJetTags");
	ftree->jet_SSVHE[ij] = jet.bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
	ftree->jet_SSVHP[ij] = jet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
	ftree->jet_CMVA[ij] = jet.bDiscriminator("combinedMVABJetTags");

	float CSVIVF = jet.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags");
	if( CSVIVF == -1000. ) CSVIVF = jet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
	ftree->jet_CSVv2[ij] = CSVIVF;

	float CSV = jet.bDiscriminator("combinedSecondaryVertexBJetTags");
	if( CSV == -1000. ) CSV = jet.bDiscriminator("pfCombinedSecondaryVertexBJetTags");

	ftree->jet_CSV[ij] = CSV;

	ftree->jet_flavour[ij] = jet.partonFlavour();

	ftree->jet_neutralHadronEnergy[ij] = jet.neutralHadronEnergy();
	ftree->jet_neutralEmEnergy[ij] = jet.neutralEmEnergy();
	ftree->jet_chargedHadronEnergy[ij] = jet.chargedHadronEnergy();
	ftree->jet_chargedEmEnergy[ij] = jet.chargedEmEnergy();
	ftree->jet_electronEnergy[ij] = jet.electronEnergy();
	ftree->jet_muonEnergy[ij] = jet.muonEnergy();
	ftree->jet_photonEnergy[ij] = jet.photonEnergy();

	ftree->jet_pileupJetId[ij] = jet.userFloat("pileupJetId:fullDiscriminant");

	const reco::GenJet* genJet = jet.genJet();
	if( genJet )
	  {
	     ftree->jet_gen_pt.push_back(genJet->pt());
	     ftree->jet_gen_eta.push_back(genJet->eta());
	     ftree->jet_gen_phi.push_back(genJet->phi());
	     ftree->jet_gen_m.push_back(genJet->mass());
	     ftree->jet_gen_E.push_back(genJet->energy());

	     ftree->jet_gen_status.push_back(genJet->status());
	     ftree->jet_gen_id.push_back(genJet->pdgId());
	  }
     }

   this->KeepEvent();
   ftree->tree->Fill();
   
   delete mc_truth;
}

void FlatTreeProducer::beginJob()
{
}

void FlatTreeProducer::endJob()
{
}

void FlatTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FlatTreeProducer);

