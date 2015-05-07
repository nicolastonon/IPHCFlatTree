#include <memory>
#include <iostream>
#include <sstream>

#include "FWCore/Common/interface/TriggerNames.h"
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
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "IPHCFlatTree/FlatTreeProducer/interface/tinyxml2.h"

#include "IPHCFlatTree/FlatTreeProducer/interface/FlatTree.hh"
#include "IPHCFlatTree/FlatTreeProducer/interface/MCTruth.hh"

#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"

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

   EGammaMvaEleEstimatorCSA14* elecMVA;
   std::vector<std::string> elecMVACatWeights;

   TMVA::Reader* mu_reader_high_b;
   TMVA::Reader* mu_reader_high_e;
   TMVA::Reader* mu_reader_low;
   TMVA::Reader* mu_reader_medium_b;
   TMVA::Reader* mu_reader_medium_e;
   TMVA::Reader* ele_reader_high_cb;
   TMVA::Reader* ele_reader_high_fb;
   TMVA::Reader* ele_reader_high_ec;
   TMVA::Reader* ele_reader_medium_cb;
   TMVA::Reader* ele_reader_medium_fb;
   TMVA::Reader* ele_reader_medium_ec;
   TMVA::Reader* ele_reader_low;

   float lepMVA_neuRelIso;
   float lepMVA_chRelIso;
   float lepMVA_jetDR;
   float lepMVA_jetPtRatio;
   float lepMVA_jetBTagCSV;
   float lepMVA_sip3d;   
   float lepMVA_dxy;
   float lepMVA_dz;
   float lepMVA_mvaId;
   
   XMLDocument xmlconf;
   
   std::string dataFormat_;
   bool isData_;
   
   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
   edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

   edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
   edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
   edm::EDGetTokenT<pat::MuonCollection> muonToken_;
   edm::EDGetTokenT<pat::TauCollection> tauToken_;
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
   if( !name.compare("n_presel_jets") )
     {
	ftree->n_presel_jets += 1;
     }
   else if( !name.compare("n_presel_electron") )
     {
	ftree->n_presel_electron += 1;
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
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_JP") )
     {
	for( unsigned int i=0;i<ftree->jet_JP.size();++i )
	  if( ftree->jet_JP[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_TCHP") )
     {
	for( unsigned int i=0;i<ftree->jet_TCHP.size();++i )
	  if( ftree->jet_TCHP[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_TCHE") )
     {
	for( unsigned int i=0;i<ftree->jet_TCHE.size();++i )
	  if( ftree->jet_TCHE[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_SSVHP") )
     {
	for( unsigned int i=0;i<ftree->jet_SSVHP.size();++i )
	  if( ftree->jet_SSVHP[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_SSVHE") )
     {
	for( unsigned int i=0;i<ftree->jet_SSVHE.size();++i )
	  if( ftree->jet_SSVHE[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_CMVA") )
     {
	for( unsigned int i=0;i<ftree->jet_CMVA.size();++i )
	  if( ftree->jet_CMVA[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_CSV") )
     {
	for( unsigned int i=0;i<ftree->jet_CSV.size();++i )
	  if( ftree->jet_CSV[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_CSVv2") )
     {
	for( unsigned int i=0;i<ftree->jet_CSVv2.size();++i )
	  if( ftree->jet_CSVv2[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_flavour") )
     {
	for( unsigned int i=0;i<ftree->jet_flavour.size();++i )
	  if( ftree->jet_flavour[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
}

// FIXME : refactorize CheckJet, CheckElectron and CheckMuon
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
						AddValue("n_presel_jets");
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
						AddValue("n_presel_jets");
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
					   AddValue("n_presel_electron");
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
					   AddValue("n_presel_electron");
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
                                    AddValue("n_presel_electron");
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
                                    AddValue("n_presel_electron");
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
   
   reader->AddVariable("LepGood_relIso03-LepGood_chargedHadRelIso03", &lepMVA_neuRelIso);
   reader->AddVariable("LepGood_chargedHadRelIso03",                  &lepMVA_chRelIso);
   reader->AddVariable("min(LepGood_jetDR,0.5)",                      &lepMVA_jetDR);
   reader->AddVariable("min(LepGood_jetPtRatio,1.5)",                 &lepMVA_jetPtRatio);
   reader->AddVariable("max(LepGood_jetBTagCSV,0)",                   &lepMVA_jetBTagCSV);
   reader->AddVariable("LepGood_sip3d",                               &lepMVA_sip3d);
   reader->AddVariable("log(abs(LepGood_dxy))",                       &lepMVA_dxy);
   reader->AddVariable("log(abs(LepGood_dz))",                        &lepMVA_dz);
   if( type == "ele" ) reader->AddVariable("LepGood_mvaIdPhys14",     &lepMVA_mvaId);
   else reader->AddVariable("LepGood_segmentCompatibility",           &lepMVA_mvaId);

   reader->BookMVA("BDTG method", basePath+"/"+weightFileName);
   
   return reader;
}

FlatTreeProducer::FlatTreeProducer(const edm::ParameterSet& iConfig)
{
   // ###
   // Temporarily redirecting stdout to avoid huge TMVA loading dump
   // ###
   cout << "Temporarily redirecting stdout to avoid huge TMVA dump when loading MVA readers..." << endl;
   stringstream tmpBuffer;
   streambuf* oldStdout = cout.rdbuf(tmpBuffer.rdbuf());

   // ###############
   // #  Load MVAs  #
   // ###############

   string CMSSW_BASE(getenv("CMSSW_BASE")); 
   elecMVA = new EGammaMvaEleEstimatorCSA14();
   
   string EGammaElectronToolsPath = CMSSW_BASE+"/src/EgammaAnalysis/ElectronTools/";
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml");
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml");
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml");
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml");
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml");
   elecMVACatWeights.push_back(EGammaElectronToolsPath+"/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml");

   elecMVA->initialize("BDT", EGammaMvaEleEstimatorCSA14::kNonTrigPhys14, true, elecMVACatWeights);

   string FlatTreeProducerLepMVAPath = CMSSW_BASE+"/src/IPHCFlatTree/FlatTreeProducer/data/lepMVA/";
   mu_reader_high_b      = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/mu_pteta_high_b_BDTG.weights.xml" ,  "mu");
   mu_reader_high_e      = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/mu_pteta_high_e_BDTG.weights.xml" ,  "mu");
   mu_reader_low         = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/mu_pteta_low_BDTG.weights.xml"  ,  "mu");
   mu_reader_medium_b    = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/mu_pteta_medium_b_BDTG.weights.xml" ,  "mu");
   mu_reader_medium_e    = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/mu_pteta_medium_e_BDTG.weights.xml" ,  "mu");
   ele_reader_high_cb    = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_high_cb_BDTG.weights.xml", "ele");
   ele_reader_high_fb    = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_high_fb_BDTG.weights.xml", "ele");
   ele_reader_high_ec    = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_high_ec_BDTG.weights.xml", "ele");
   ele_reader_medium_cb  = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_medium_cb_BDTG.weights.xml" , "ele");
   ele_reader_medium_fb  = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_medium_fb_BDTG.weights.xml" , "ele");
   ele_reader_medium_ec  = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_medium_ec_BDTG.weights.xml" , "ele");
   ele_reader_low        = BookLeptonMVAReader(FlatTreeProducerLepMVAPath, "/el_pteta_low_BDTG.weights.xml", "ele");
   
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
   if (buffersize <= 0) buffersize = 32000;
   ftree->CreateBranches(buffersize);
   
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
   triggerObjects_    = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"));
   triggerPrescales_  = consumes<pat::PackedTriggerPrescales>(edm::InputTag(std::string("patTrigger")));
   vertexToken_       = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexInput"));
   electronToken_     = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronInput"));
   muonToken_         = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonInput"));
   tauToken_          = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tauInput"));
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

   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerObjects_, triggerObjects);

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

   // Taus
   edm::Handle<pat::TauCollection> taus;
   iEvent.getByToken(tauToken_,taus);
   
   // Conversions info
   edm::Handle<reco::ConversionCollection> hConversions;
   if( dataFormat_ != "AOD" ) iEvent.getByLabel("reducedEgamma","reducedConversions",hConversions);
  
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
   for(std::vector<PileupSummaryInfo>::const_iterator pvi=pileupInfo->begin();
       pvi!=pileupInfo->end();pvi++)
     {
       signed int n_bc = pvi->getBunchCrossing();
       ftree->mc_pu_BunchCrossing.push_back(n_bc);
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

       ftree->mc_pu_Nzpositions.push_back(pvi->getPU_zpositions().size());
       for( unsigned int ipu=0;ipu<pvi->getPU_zpositions().size();ipu++ )
	 {
	   mc_pu_zpositions.push_back((pvi->getPU_zpositions())[ipu]);
	   mc_pu_sumpT_lowpT.push_back((pvi->getPU_sumpT_lowpT())[ipu]);
	   mc_pu_sumpT_highpT.push_back((pvi->getPU_sumpT_highpT())[ipu]);
	   mc_pu_ntrks_lowpT.push_back((pvi->getPU_ntrks_lowpT())[ipu]);
	   mc_pu_ntrks_highpT.push_back((pvi->getPU_ntrks_highpT())[ipu]);
	 }

       ftree->mc_pu_zpositions.push_back(mc_pu_zpositions);
       ftree->mc_pu_sumpT_lowpT.push_back(mc_pu_sumpT_lowpT);
       ftree->mc_pu_sumpT_highpT.push_back(mc_pu_sumpT_highpT);
       ftree->mc_pu_ntrks_lowpT.push_back(mc_pu_ntrks_lowpT);
       ftree->mc_pu_ntrks_highpT.push_back(mc_pu_ntrks_highpT);
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
   bool do_mc_truth_ttz = ftree->doWrite("mc_truth_ttz");
   bool do_mc_truth_ttw = ftree->doWrite("mc_truth_ttw");
   bool do_mc_truth_tzq = ftree->doWrite("mc_truth_tzq");
   bool do_mc_truth_thq = ftree->doWrite("mc_truth_thq");

   MCTruth *mc_truth = new MCTruth();

   bool reqMCTruth = 0;
   if( (
	do_mc_truth_tth || 
	do_mc_truth_tzq ||
        do_mc_truth_ttz ||
	do_mc_truth_ttw ||
	do_mc_truth_thq
       ) &&
       !isData_ )
     {
	mc_truth->Init(*ftree);
	if( do_mc_truth_tth ) mc_truth->fillTTHSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
	if( do_mc_truth_ttz ) mc_truth->fillTTZSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
	if( do_mc_truth_ttw ) mc_truth->fillTTWSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
	if( do_mc_truth_tzq ) mc_truth->fillTZQSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
	if( do_mc_truth_thq ) mc_truth->fillTHQSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
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
   //triggerIdentifiers_.push_back("HLT_Ele27_WP80_v*");
   triggerIdentifiers_.push_back("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1");
   triggerIdentifiers_.push_back("HLT_Mu17_Mu8_v1");
   triggerIdentifiers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1");
   triggerIdentifiers_.push_back("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1");

   for (unsigned int j = 0; j < triggerIdentifiers_.size(); ++j)
   {
	std::string idName = triggerIdentifiers_[j];
	std::string idNameUnstarred = idName;
	bool isStarred = (idName.find("*")!=std::string::npos);
	if( isStarred ) idNameUnstarred.erase( idName.find("*"), 1 );

	  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
	  {
	     //if( (isStarred && names.triggerName(i).find(idNameUnstarred)!=std::string::npos ) ||
		 //(!isStarred && names.triggerName(i)==idName) )
		 //{
	     //   std::cout << "[" << i << "] " << (triggerBits->accept(i) ? "1" : "0") << "  " << names.triggerName(i)  << std::endl;
		 //}

         ftree->trigger.push_back(i);
         ftree->trigger_pass.push_back(triggerBits->accept(i) ? true : false);
         ftree->trigger_prescale.push_back(triggerPrescales->getPrescaleForIndex(i));
	  }
   }
 
   //std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
   for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
        obj.unpackPathNames(names);
        //std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
        // Print trigger object collection and type
        //std::cout << "\t   Collection: " << obj.collection() << std::endl;
        //std::cout << "\t   Type IDs:   ";

        for (unsigned h = 0; h < obj.filterIds().size(); ++h)
        {
            ftree->triggerobject_isTriggerL1Mu.push_back(obj.filterIds()[h] == -81 ? true : false);
            ftree->triggerobject_isTriggerL1NoIsoEG.push_back(obj.filterIds()[h] == -82 ? true : false);
            ftree->triggerobject_isTriggerL1IsoEG.push_back(obj.filterIds()[h] == -83 ? true : false);
            ftree->triggerobject_isTriggerL1CenJet.push_back(obj.filterIds()[h] == -84 ? true : false);
            ftree->triggerobject_isTriggerL1ForJet.push_back(obj.filterIds()[h] == -85 ? true : false);
            ftree->triggerobject_isTriggerL1TauJet.push_back(obj.filterIds()[h] == -86 ? true : false);
            ftree->triggerobject_isTriggerL1ETM.push_back(obj.filterIds()[h] == -87 ? true : false);
            ftree->triggerobject_isTriggerL1ETT.push_back(obj.filterIds()[h] == -88 ? true : false);
            ftree->triggerobject_isTriggerL1HTT.push_back(obj.filterIds()[h] == -89 ? true : false);
            ftree->triggerobject_isTriggerL1HTM.push_back(obj.filterIds()[h] == -90 ? true : false);
            ftree->triggerobject_isTriggerL1JetCounts.push_back(obj.filterIds()[h] == -91 ? true : false);
            ftree->triggerobject_isTriggerL1HfBitCounts.push_back(obj.filterIds()[h] == -92 ? true : false);
            ftree->triggerobject_isTriggerL1HfRingEtSums.push_back(obj.filterIds()[h] == -93 ? true : false);
            ftree->triggerobject_isTriggerL1TechTrig.push_back(obj.filterIds()[h] == -94 ? true : false);
            ftree->triggerobject_isTriggerL1Castor.push_back(obj.filterIds()[h] == -95 ? true : false);
            ftree->triggerobject_isTriggerL1BPTX.push_back(obj.filterIds()[h] == -96 ? true : false);
            ftree->triggerobject_isTriggerL1GtExternal.push_back(obj.filterIds()[h] == -97 ? true : false);

            ftree->triggerobject_isHLT_TriggerPhoton.push_back(obj.filterIds()[h] == 81 ? true : false);
            ftree->triggerobject_isHLT_TriggerElectron.push_back(obj.filterIds()[h] == 82 ? true : false);
            ftree->triggerobject_isHLT_TriggerMuon.push_back(obj.filterIds()[h] == 83 ? true : false);
            ftree->triggerobject_isHLT_TriggerTau.push_back(obj.filterIds()[h] == 84 ? true : false);
            ftree->triggerobject_isHLT_TriggerJet.push_back(obj.filterIds()[h] == 85 ? true : false);
            ftree->triggerobject_isHLT_TriggerBJet.push_back(obj.filterIds()[h] == 86 ? true : false);
            ftree->triggerobject_isHLT_TriggerMET.push_back(obj.filterIds()[h] == 87 ? true : false);
            ftree->triggerobject_isHLT_TriggerTET.push_back(obj.filterIds()[h] == 88 ? true : false);
            ftree->triggerobject_isHLT_TriggerTHT.push_back(obj.filterIds()[h] == 89 ? true : false);
            ftree->triggerobject_isHLT_TriggerMHT.push_back(obj.filterIds()[h] == 90 ? true : false);
            ftree->triggerobject_isHLT_TriggerTrack.push_back(obj.filterIds()[h] == 91 ? true : false);
            ftree->triggerobject_isHLT_TriggerCluster.push_back(obj.filterIds()[h] == 92 ? true : false);
            ftree->triggerobject_isHLT_TriggerMETSig.push_back(obj.filterIds()[h] == 93 ? true : false);
            ftree->triggerobject_isHLT_TriggerELongit.push_back(obj.filterIds()[h] == 94 ? true : false);
            ftree->triggerobject_isHLT_TriggerMHTSig.push_back(obj.filterIds()[h] == 95 ? true : false);
            ftree->triggerobject_isHLT_TriggerHLongit.push_back(obj.filterIds()[h] == 96 ? true : false);
        }

        //for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
        //std::cout << std::endl;
        // Print associated trigger filters
        //std::cout << "\t   Filters:    ";
        //for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
        //std::cout << std::endl;
        std::vector<std::string> pathNamesAll  = obj.pathNames(false);
        std::vector<std::string> pathNamesLast = obj.pathNames(true);
        // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
        // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
        // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
        //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
//            bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
//            bool isL3   = obj.hasPathName( pathNamesAll[h], false, true ); 
//            bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
//            bool isNone = obj.hasPathName( pathNamesAll[h], false, false ); 
            //std::cout << "   " << pathNamesAll[h];
            //if (isBoth) std::cout << "(L,3)";
            //if (isL3 && !isBoth) std::cout << "(*,3)";
            //if (isLF && !isBoth) std::cout << "(L,*)";
            //if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
        }
        //std::cout << std::endl;
        
        ftree->triggerobject_pt.push_back(obj.pt());
        ftree->triggerobject_eta.push_back(obj.eta());
        ftree->triggerobject_phi.push_back(obj.phi());
    }
//    std::cout << std::endl;
   
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
	
	ftree->pv_ndof = primVtx->ndof();
	ftree->pv_rho = primVtx->position().Rho();
	ftree->pv_isFake = primVtx->isFake();
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
   
   for(int ie=0;ie<nElec;ie++)
     {
	const pat::Electron& elec = electrons->at(ie);

	// Skimming electrons with pT < 5 GeV.
	if (elec.pt() < 5) continue;
        ftree->el_pt.push_back(elec.pt());
        ftree->el_eta.push_back(elec.eta());
        ftree->el_phi.push_back(elec.phi());
        ftree->el_m.push_back(elec.mass());
        ftree->el_E.push_back(elec.energy());
        ftree->el_id.push_back(elec.pdgId());
        ftree->el_charge.push_back(elec.charge());

        ftree->el_scleta.push_back(elec.superCluster()->eta());
        ftree->el_isGsfCtfScPixChargeConsistent.push_back(elec.isGsfCtfScPixChargeConsistent());
        ftree->el_sigmaIetaIeta.push_back(elec.sigmaIetaIeta());
        ftree->el_sigmaIphiIphi.push_back(elec.sigmaIphiIphi());
        ftree->el_hadronicOverEm.push_back(elec.hadronicOverEm());
        ftree->el_dr03TkSumPt.push_back(elec.dr03TkSumPt());
        ftree->el_dr03EcalRecHitSumEt.push_back(elec.dr03EcalRecHitSumEt());
        ftree->el_dr03HcalTowerSumEt.push_back(elec.dr03HcalTowerSumEt());

	int numberOfLostHits = -666;

	if( elec.gsfTrack().isNonnull() )
	  {
	     numberOfLostHits = elec.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::HitCategory::MISSING_INNER_HITS);
	  }

	ftree->el_numberOfLostHits.push_back(numberOfLostHits);

	ftree->el_fbrem.push_back(elec.fbrem());

        ftree->el_deltaEtaSuperClusterTrackAtVtx.push_back(elec.deltaEtaSuperClusterTrackAtVtx());
        ftree->el_deltaPhiSuperClusterTrackAtVtx.push_back(elec.deltaPhiSuperClusterTrackAtVtx());
        ftree->el_deltaEtaSeedClusterTrackAtCalo.push_back(elec.deltaEtaSeedClusterTrackAtCalo());

        ftree->el_see.push_back(elec.full5x5_sigmaIetaIeta());
        ftree->el_spp.push_back(elec.full5x5_sigmaIphiIphi());

        ftree->el_superClusterEtaWidth.push_back(elec.superCluster()->etaWidth());
        ftree->el_superClusterPhiWidth.push_back(elec.superCluster()->phiWidth());

	double OneMinusE1x5E5x5 = (elec.e5x5() != 0.) ? 1.-(elec.e1x5()/elec.e5x5()) : -1.;
	double full5x5_OneMinusE1x5E5x5 = (elec.full5x5_e5x5() != 0.) ? 1.-(elec.full5x5_e1x5()/elec.full5x5_e5x5()) : -1.;

        ftree->el_full5x5_OneMinusE1x5E5x5.push_back(full5x5_OneMinusE1x5E5x5);
        ftree->el_OneMinusE1x5E5x5.push_back(OneMinusE1x5E5x5);

        ftree->el_full5x5_r9.push_back(elec.full5x5_r9());
        ftree->el_r9.push_back(elec.r9());

        ftree->el_eSuperClusterOverP.push_back(elec.eSuperClusterOverP());

        double IoEmIoP = (1.0/elec.ecalEnergy())-(1.0/elec.p());
        ftree->el_IoEmIoP.push_back(IoEmIoP);
        ftree->el_eleEoPout.push_back(elec.eEleClusterOverPout());
        double PreShowerOverRaw = elec.superCluster()->preshowerEnergy()/elec.superCluster()->rawEnergy();
        ftree->el_PreShowerOverRaw.push_back(PreShowerOverRaw);
	ftree->el_ecalEnergy.push_back(elec.ecalEnergy());
	
	// https://github.com/gpetruc/cmg-cmssw/blob/CMG_MiniAOD_Lite_V6_0_from-CMSSW_7_0_6/EgammaAnalysis/ElectronTools/src/EGammaMvaEleEstimator.cc#L1244-1336
        ftree->el_dB3D.push_back(elec.dB(pat::Electron::PV3D));
        ftree->el_edB3D.push_back(elec.edB(pat::Electron::PV3D));

	double mvaNonTrigV0 = elecMVA->mvaValue(elec,
						false);
	
        ftree->el_mvaNonTrigV0.push_back(mvaNonTrigV0);

        ftree->el_neutralHadronIso.push_back(elec.neutralHadronIso());
        ftree->el_chargedHadronIso.push_back(elec.chargedHadronIso());
        ftree->el_puChargedHadronIso.push_back(elec.puChargedHadronIso());
        ftree->el_ecalIso.push_back(elec.ecalIso());
        ftree->el_hcalIso.push_back(elec.hcalIso());
        ftree->el_particleIso.push_back(elec.particleIso());
        ftree->el_photonIso.push_back(elec.photonIso());
        ftree->el_trackIso.push_back(elec.trackIso());

        ftree->el_pfIso_sumChargedHadronPt.push_back(elec.pfIsolationVariables().sumChargedHadronPt);
        ftree->el_pfIso_sumNeutralHadronEt.push_back(elec.pfIsolationVariables().sumNeutralHadronEt);
        ftree->el_pfIso_sumPhotonEt.push_back(elec.pfIsolationVariables().sumPhotonEt);
        ftree->el_pfIso_sumPUPt.push_back(elec.pfIsolationVariables().sumPUPt);

	// mini-iso
	float miniIso = -666;
	if( dataFormat_ != "AOD" )
	  {
	     miniIso = getPFIsolation(pfcands,dynamic_cast<const reco::Candidate*>(&elec),0.05,0.2,10.,false,false);
	  }
	ftree->el_miniIso.push_back(miniIso);
        ftree->el_isLoose.push_back(elec.electronID("eidLoose"));
        ftree->el_isTight.push_back(elec.electronID("eidTight"));
        ftree->el_isRobustLoose.push_back(elec.electronID("eidRobustLoose"));
        ftree->el_isRobustTight.push_back(elec.electronID("eidRobustTight"));
        ftree->el_isRobustHighEnergy.push_back(elec.electronID("eidRobustHighEnergy"));

        ftree->el_vx.push_back(elec.vx());
        ftree->el_vy.push_back(elec.vy());
        ftree->el_vz.push_back(elec.vz());

        ftree->el_isGsf.push_back(elec.gsfTrack().isNonnull());	
	
        ftree->el_passConversionVeto.push_back(elec.passConversionVeto());
	       
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

        ftree->el_numberOfHits.push_back(numberOfHits);
        ftree->el_dxy.push_back(dxy);
        ftree->el_dz.push_back(dz);
        ftree->el_dxyError.push_back(dxyError);
        ftree->el_dzError.push_back(dzError);

	double el_pt = elec.pt();
	double el_eta = elec.eta();
	double el_scleta = ftree->el_scleta.back();
	double el_phi = elec.phi();
	double el_lepMVA = -666.;

	double relIso = (ftree->el_pfIso_sumChargedHadronPt.back() +
			 std::max(ftree->el_pfIso_sumNeutralHadronEt.back() +
				  ftree->el_pfIso_sumPhotonEt.back() -
				  ftree->el_pfIso_sumPUPt.back()/2.,0.0))/el_pt;

	lepMVA_neuRelIso = relIso - ftree->el_chargedHadronIso.back()/el_pt;
	lepMVA_chRelIso = ftree->el_chargedHadronIso.back()/el_pt;
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
	lepMVA_sip3d = fabs(ftree->el_dB3D.back()/ftree->el_edB3D.back());
	lepMVA_dxy = log(fabs(ftree->el_dxy.back()));
	lepMVA_dz = log(fabs(ftree->el_dz.back()));
	lepMVA_mvaId = mvaNonTrigV0;

	if( el_pt < 10 ) el_lepMVA = ele_reader_low->EvaluateMVA("BDTG method");
	else if( el_pt >= 10 && el_pt < 25 )
	  {
	     if( fabs(el_scleta) < 0.8 ) el_lepMVA = ele_reader_medium_cb->EvaluateMVA("BDTG method");
	     else if( fabs(el_scleta) >= 0.8 && fabs(el_scleta) < 1.479 ) el_lepMVA = ele_reader_medium_fb->EvaluateMVA("BDTG method");
	     else if( fabs(el_scleta) >= 1.479 ) el_lepMVA = ele_reader_medium_ec->EvaluateMVA("BDTG method");
	  }	
	else if( el_pt >= 25 )
	  {
	     if( fabs(el_scleta) < 0.8 ) el_lepMVA = ele_reader_high_cb->EvaluateMVA("BDTG method");
	     else if( fabs(el_scleta) >= 0.8 && fabs(el_scleta) < 1.479 ) el_lepMVA = ele_reader_high_fb->EvaluateMVA("BDTG method");
	     else if( fabs(el_scleta) >= 1.479 ) el_lepMVA = ele_reader_high_ec->EvaluateMVA("BDTG method");
	  }	
	
	ftree->el_lepMVA.push_back(el_lepMVA);

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

        ftree->el_lepMVA_neuRelIso.push_back(lepMVA_neuRelIso);
        ftree->el_lepMVA_chRelIso.push_back(lepMVA_chRelIso);
        ftree->el_lepMVA_jetDR.push_back(lepMVA_jetDR);
        ftree->el_lepMVA_jetPtRatio.push_back(lepMVA_jetPtRatio);
        ftree->el_lepMVA_jetBTagCSV.push_back(lepMVA_jetBTagCSV);
        ftree->el_lepMVA_sip3d.push_back(lepMVA_sip3d);
	ftree->el_lepMVA_dxy.push_back(lepMVA_dxy);
	ftree->el_lepMVA_dz.push_back(lepMVA_dz);
        ftree->el_lepMVA_mvaId.push_back(lepMVA_mvaId);

	bool allowCkfMatch = true;
	float lxyMin = 2.0;
	float probMin = 1e-6;
	uint nHitsBeforeVtxMax = 0;

	if( dataFormat_ != "AOD" )
	  {	     
	     bool matchConv = 0;
	     if( &beamspot ) matchConv = ConversionTools::hasMatchedConversion(elec,hConversions,beamspot.position(),allowCkfMatch,lxyMin,probMin,nHitsBeforeVtxMax);
	     ftree->el_hasMatchedConversion.push_back(matchConv);
	  }	
     }   
   ftree->el_n = ftree->el_pt.size();
   
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
   for(int im=0;im<nMuon;im++)
     {
	const pat::Muon& muon = muons->at(im);
	
	// Skimming muons with pT < 5 GeV.
	if (muon.pt() < 5) continue;

        ftree->mu_pt.push_back(muon.pt());
        ftree->mu_eta.push_back(muon.eta());
        ftree->mu_phi.push_back(muon.phi());
        ftree->mu_m.push_back(muon.mass());
        ftree->mu_E.push_back(muon.energy());
        ftree->mu_id.push_back(muon.pdgId());
        ftree->mu_charge.push_back(muon.charge());

        ftree->mu_dB3D.push_back(muon.dB(pat::Muon::PV3D));
        ftree->mu_edB3D.push_back(muon.edB(pat::Muon::PV3D));

        ftree->mu_neutralHadronIso.push_back(muon.neutralHadronIso());
        ftree->mu_chargedHadronIso.push_back(muon.chargedHadronIso());
        ftree->mu_puChargedHadronIso.push_back(muon.puChargedHadronIso());
        ftree->mu_ecalIso.push_back(muon.ecalIso());
        ftree->mu_hcalIso.push_back(muon.hcalIso());
        ftree->mu_photonIso.push_back(muon.photonIso());
        ftree->mu_trackIso.push_back(muon.trackIso());

        ftree->mu_pfIso03_sumChargedHadronPt.push_back(muon.pfIsolationR03().sumChargedHadronPt);
        ftree->mu_pfIso03_sumNeutralHadronEt.push_back(muon.pfIsolationR03().sumNeutralHadronEt);
        ftree->mu_pfIso03_sumPhotonEt.push_back(muon.pfIsolationR03().sumPhotonEt);
        ftree->mu_pfIso03_sumPUPt.push_back(muon.pfIsolationR03().sumPUPt);

        ftree->mu_pfIso04_sumChargedHadronPt.push_back(muon.pfIsolationR04().sumChargedHadronPt);
        ftree->mu_pfIso04_sumNeutralHadronEt.push_back(muon.pfIsolationR04().sumNeutralHadronEt);
        ftree->mu_pfIso04_sumPhotonEt.push_back(muon.pfIsolationR04().sumPhotonEt);
        ftree->mu_pfIso04_sumPUPt.push_back(muon.pfIsolationR04().sumPUPt);
	
        ftree->mu_isGlobalMuon.push_back(muon.isGlobalMuon());
        ftree->mu_isTrackerMuon.push_back(muon.isTrackerMuon());
        ftree->mu_isStandAloneMuon.push_back(muon.isStandAloneMuon());
        ftree->mu_isCaloMuon.push_back(muon.isCaloMuon());
        ftree->mu_isPFMuon.push_back(muon.isPFMuon());

        ftree->mu_vx.push_back(muon.vx());
        ftree->mu_vy.push_back(muon.vy());
        ftree->mu_vz.push_back(muon.vz());

        ftree->mu_hasGlobalTrack.push_back(muon.globalTrack().isNonnull());
        ftree->mu_hasInnerTrack.push_back(muon.innerTrack().isNonnull());
	ftree->mu_hasTrack.push_back(muon.track().isNonnull());

	// mini-iso
	float miniIso = -666;
	if( dataFormat_ != "AOD" )
	{
	  miniIso = getPFIsolation(pfcands,dynamic_cast<const reco::Candidate*>(&muon),0.05,0.2,10.,false,false);
	}
	ftree->mu_miniIso.push_back(miniIso);

	ftree->mu_segmentCompatibility.push_back(muon.segmentCompatibility());

	float globalTrack_dxy = -10E+10;
	float globalTrack_dz = -10E+10;
	float globalTrack_dxyError = -666.;
	float globalTrack_dzError = -666.;	
	float globalTrack_normalizedChi2 = -666.;
	
	float combinedQuality_chi2LocalPosition = -666.;
	float combinedQuality_trkKink = -666.;

	float innerTrack_dxy = -10E+10;
	float innerTrack_dz = -10E+10;
	float innerTrack_dxyError = -666.;
	float innerTrack_dzError = -666.;
	float innerTrack_normalizedChi2 = -666.;
	
	float innerTrack_validFraction = -666.;

	float bestTrack_dxy = -10E+10;
	float bestTrack_dz = -10E+10;
	float bestTrack_dxyError = -666.;
	float bestTrack_dzError = -666.;
	float bestTrack_normalizedChi2 = -666.;

	bool isTightMuon = 0;
	if( primVtx ) isTightMuon = muon.isTightMuon(*primVtx);
	ftree->mu_isTightMuon.push_back(isTightMuon);

	ftree->mu_bestTrackType.push_back(muon.bestTrackType());
	
	if( muon.globalTrack().isNonnull() )
	  {
	     globalTrack_dxy = muon.globalTrack()->dxy(muon.globalTrack()->referencePoint());
	     globalTrack_dz = muon.globalTrack()->dz(muon.globalTrack()->referencePoint());
	     globalTrack_dxyError = muon.globalTrack()->dxyError();
	     globalTrack_dzError = muon.globalTrack()->dzError();	     
	     globalTrack_normalizedChi2 = muon.globalTrack()->normalizedChi2();
	     
	     combinedQuality_chi2LocalPosition = muon.combinedQuality().chi2LocalPosition;
	     combinedQuality_trkKink = muon.combinedQuality().trkKink;
	  }
	if( muon.innerTrack().isNonnull() )
	  {
	     if( primVtx ) innerTrack_dxy = muon.innerTrack()->dxy(primVtx->position());
	     if( primVtx ) innerTrack_dz = muon.innerTrack()->dz(primVtx->position());
	     innerTrack_dxyError = muon.innerTrack()->dxyError();
	     innerTrack_dzError = muon.innerTrack()->dzError();
	     innerTrack_normalizedChi2 = muon.innerTrack()->normalizedChi2();
	     
	     innerTrack_validFraction = muon.innerTrack()->validFraction();
	  }

	bestTrack_dxy = muon.bestTrack()->dxy(muon.bestTrack()->referencePoint());
	bestTrack_dz = muon.bestTrack()->dz(muon.bestTrack()->referencePoint());
	bestTrack_dxyError = muon.bestTrack()->dxyError();
	bestTrack_dzError = muon.bestTrack()->dzError();
	bestTrack_normalizedChi2 = muon.bestTrack()->normalizedChi2();

        ftree->mu_globalTrack_dxy.push_back(globalTrack_dxy);
        ftree->mu_globalTrack_dz.push_back(globalTrack_dz);
        ftree->mu_globalTrack_dxyError.push_back(globalTrack_dxyError);
        ftree->mu_globalTrack_dzError.push_back(globalTrack_dzError);	
	ftree->mu_globalTrack_normalizedChi2.push_back(globalTrack_normalizedChi2);
	
	ftree->mu_combinedQuality_chi2LocalPosition.push_back(combinedQuality_chi2LocalPosition);
	ftree->mu_combinedQuality_trkKink.push_back(combinedQuality_trkKink);

        ftree->mu_innerTrack_dxy.push_back(innerTrack_dxy);
        ftree->mu_innerTrack_dz.push_back(innerTrack_dz);
        ftree->mu_innerTrack_dxyError.push_back(innerTrack_dxyError);
        ftree->mu_innerTrack_dzError.push_back(innerTrack_dzError);
	ftree->mu_innerTrack_normalizedChi2.push_back(innerTrack_normalizedChi2);
	
	ftree->mu_innerTrack_validFraction.push_back(innerTrack_validFraction);

        ftree->mu_bestTrack_dxy.push_back(bestTrack_dxy);
        ftree->mu_bestTrack_dz.push_back(bestTrack_dz);
        ftree->mu_bestTrack_dxyError.push_back(bestTrack_dxyError);
        ftree->mu_bestTrack_dzError.push_back(bestTrack_dzError);
	ftree->mu_bestTrack_normalizedChi2.push_back(bestTrack_normalizedChi2);

	int trackerLayersWithMeasurement = -666;
	if( muon.track().isNonnull() )
	  {
	     trackerLayersWithMeasurement = muon.track()->hitPattern().trackerLayersWithMeasurement();
	  }	
	
	ftree->mu_track_trackerLayersWithMeasurement.push_back(trackerLayersWithMeasurement);
	
	float innerTrack_pt = -666.;
	float innerTrack_ptError = -666.;
	int innerTrack_numberOfValidPixelHits = -666;
		
	if( muon.innerTrack().isNonnull() )
	  {
	     innerTrack_pt = muon.innerTrack()->pt();
	     innerTrack_ptError = muon.innerTrack()->ptError();
	     innerTrack_numberOfValidPixelHits = muon.innerTrack()->hitPattern().numberOfValidPixelHits();
	  }

        ftree->mu_innerTrack_pt.push_back(innerTrack_pt);
        ftree->mu_innerTrack_ptError.push_back(innerTrack_ptError);
	ftree->mu_innerTrack_numberOfValidPixelHits.push_back(innerTrack_numberOfValidPixelHits);

        ftree->mu_numberOfMatches.push_back(muon.numberOfMatches());
	
	ftree->mu_numberOfMatchedStations.push_back(muon.numberOfMatchedStations());

	int numberOfValidMuonHits = -666;

	if( muon.globalTrack().isNonnull() )
	  {
	     numberOfValidMuonHits = muon.globalTrack()->hitPattern().numberOfValidMuonHits();
	  }

        ftree->mu_numberOfValidMuonHits.push_back(numberOfValidMuonHits);

	double mu_pt = muon.pt();
	double mu_eta = muon.eta();
	double mu_phi = muon.phi();
	double mu_lepMVA = -666.;

	double relIso = (muon.pfIsolationR03().sumChargedHadronPt +
			 std::max(muon.pfIsolationR03().sumNeutralHadronEt +
				  muon.pfIsolationR03().sumPhotonEt -
				  muon.pfIsolationR03().sumPUPt/2.,0.0))/mu_pt;

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
	lepMVA_sip3d = fabs(ftree->mu_dB3D.back()/ftree->mu_edB3D.back());
	lepMVA_dxy = log(fabs(ftree->mu_innerTrack_dxy.back()));
	lepMVA_dz = log(fabs(ftree->mu_innerTrack_dz.back()));
	lepMVA_mvaId = ftree->mu_segmentCompatibility.back();

	if( mu_pt < 10 ) mu_lepMVA = mu_reader_low->EvaluateMVA("BDTG method");
	else if( mu_pt >= 10 && mu_pt < 25 )
	  {
	     if( fabs(mu_eta) < 1.5 ) mu_lepMVA = mu_reader_medium_b->EvaluateMVA("BDTG method");
	     else if( fabs(mu_eta) >= 1.5 ) mu_lepMVA = mu_reader_medium_e->EvaluateMVA("BDTG method");
	  }	
	else if( mu_pt >= 25 )
	  {
	     if( fabs(mu_eta) < 1.5 ) mu_lepMVA = mu_reader_high_b->EvaluateMVA("BDTG method");
	     else if( fabs(mu_eta) >= 1.5 ) mu_lepMVA = mu_reader_high_e->EvaluateMVA("BDTG method");
	  }	

	ftree->mu_lepMVA.push_back(mu_lepMVA);
	
        ftree->mu_lepMVA_neuRelIso.push_back(lepMVA_neuRelIso);
        ftree->mu_lepMVA_chRelIso.push_back(lepMVA_chRelIso);
        ftree->mu_lepMVA_jetDR.push_back(lepMVA_jetDR);
        ftree->mu_lepMVA_jetPtRatio.push_back(lepMVA_jetPtRatio);
        ftree->mu_lepMVA_jetBTagCSV.push_back(lepMVA_jetBTagCSV);
        ftree->mu_lepMVA_sip3d.push_back(lepMVA_sip3d);
        ftree->mu_lepMVA_dxy.push_back(lepMVA_dxy);
        ftree->mu_lepMVA_dz.push_back(lepMVA_dz);
	ftree->mu_lepMVA_mvaId.push_back(lepMVA_mvaId);

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
   ftree->mu_n = ftree->mu_pt.size();

   // Taus

   int nTau = taus->size();
   for(int it=0;it<nTau;it++)
     {
	const pat::Tau& tau = taus->at(it);
	
	// Skimming taus with pT < 5 GeV. (should do nothing for miniAOD where pT > 18 GeV is applied)
	if (tau.pt() < 5) continue;

	ftree->tau_pt.push_back(tau.pt());
	ftree->tau_eta.push_back(tau.eta());
	ftree->tau_phi.push_back(tau.phi());
	ftree->tau_m.push_back(tau.mass());
	ftree->tau_E.push_back(tau.energy());
	ftree->tau_id.push_back(tau.pdgId());
	ftree->tau_charge.push_back(tau.charge());
	
	// https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Taus
	
	float tau_leadingTrackPt = -666;
	int tau_leadingTrackCharge = -666;
	
	ftree->tau_hasLeadChargedHadrCand.push_back(tau.leadChargedHadrCand().isNonnull());
	
	if( tau.leadChargedHadrCand().isNonnull() )
	  {
	    tau_leadingTrackPt = tau.leadChargedHadrCand()->pt();
	    tau_leadingTrackCharge = tau.leadChargedHadrCand()->charge();
	  }	
	ftree->tau_leadingTrackPt.push_back(tau_leadingTrackPt);
	ftree->tau_leadingTrackCharge.push_back(tau_leadingTrackCharge);
	
	ftree->tau_decayMode.push_back(tau.decayMode());
//	ftree->tau_decayModeFindingOldDMs.push_back(tau.tauID("decayModeFindingOldDMs"));
	ftree->tau_decayModeFindingNewDMs.push_back(tau.tauID("decayModeFindingNewDMs"));
	
	ftree->tau_puCorrPtSum.push_back(tau.tauID("puCorrPtSum"));
	ftree->tau_neutralIsoPtSum.push_back(tau.tauID("neutralIsoPtSum"));
	ftree->tau_chargedIsoPtSum.push_back(tau.tauID("chargedIsoPtSum"));
	ftree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
	
	ftree->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
	ftree->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
	ftree->tau_byTightCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
	
	ftree->tau_againstMuonLoose3.push_back(tau.tauID("againstMuonLoose3"));
	ftree->tau_againstMuonTight3.push_back(tau.tauID("againstMuonTight3"));
	
	ftree->tau_againstElectronVLooseMVA5.push_back(tau.tauID("againstElectronVLooseMVA5"));
	ftree->tau_againstElectronLooseMVA5.push_back(tau.tauID("againstElectronLooseMVA5"));
	ftree->tau_againstElectronMediumMVA5.push_back(tau.tauID("againstElectronMediumMVA5"));

	ftree->tau_pfEssential_jet_pt.push_back(tau.pfEssential().p4Jet_.pt());
	ftree->tau_pfEssential_jet_eta.push_back(tau.pfEssential().p4Jet_.eta());
	ftree->tau_pfEssential_jet_phi.push_back(tau.pfEssential().p4Jet_.phi());
	ftree->tau_pfEssential_jet_m.push_back(tau.pfEssential().p4Jet_.mass());

	ftree->tau_pfEssential_jetCorr_pt.push_back(tau.pfEssential().p4CorrJet_.pt());
	ftree->tau_pfEssential_jetCorr_eta.push_back(tau.pfEssential().p4CorrJet_.eta());
	ftree->tau_pfEssential_jetCorr_phi.push_back(tau.pfEssential().p4CorrJet_.phi());
	ftree->tau_pfEssential_jetCorr_m.push_back(tau.pfEssential().p4CorrJet_.mass());
	
	float tau_pfEssential_sv_x = -666;
	float tau_pfEssential_sv_y = -666;
	float tau_pfEssential_sv_z = -666;

	ftree->tau_pfEssential_hasSV.push_back(tau.pfEssential().sv_.isNonnull());
	
	if( tau.pfEssential().sv_.isNonnull() )
	  {	     
	     tau_pfEssential_sv_x = tau.pfEssential().sv_->x();
	     tau_pfEssential_sv_y = tau.pfEssential().sv_->y();
	     tau_pfEssential_sv_z = tau.pfEssential().sv_->z();
	  }
	
	ftree->tau_pfEssential_sv_x.push_back(tau_pfEssential_sv_x);
	ftree->tau_pfEssential_sv_y.push_back(tau_pfEssential_sv_y);
	ftree->tau_pfEssential_sv_z.push_back(tau_pfEssential_sv_z);

	ftree->tau_pfEssential_flightLengthSig.push_back(tau.pfEssential().flightLengthSig_);
	ftree->tau_pfEssential_dxy.push_back(tau.pfEssential().dxy_);
	ftree->tau_pfEssential_dxy_error.push_back(tau.pfEssential().dxy_error_);
	ftree->tau_pfEssential_dxy_Sig.push_back(tau.pfEssential().dxy_Sig_);

	/*	ftree->tau_pfEssential_flightLengthSig.push_back(tau.pfEssential().flightLengthSig);
	ftree->tau_pfEssential_dxy.push_back(tau.pfEssential().dxy);
	ftree->tau_pfEssential_dxy_error.push_back(tau.pfEssential().dxy_error);
	ftree->tau_pfEssential_dxy_Sig.push_back(tau.pfEssential().dxy_Sig);*/
     }   
   ftree->tau_n = ftree->tau_pt.size();   
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
   for(int ij=0;ij<nJet;ij++)
     {
	const pat::Jet& jet = jets->at(ij);

	ftree->jet_pt.push_back(jet.pt());
	ftree->jet_eta.push_back(jet.eta());
	ftree->jet_phi.push_back(jet.phi());
	ftree->jet_m.push_back(jet.mass());
	ftree->jet_E.push_back(jet.energy());
	ftree->jet_JBP.push_back(jet.bDiscriminator("jetBProbabilityBJetTags"));
	ftree->jet_JP.push_back(jet.bDiscriminator("jetProbabilityBJetTags"));
	ftree->jet_TCHP.push_back(jet.bDiscriminator("trackCountingHighPurBJetTags"));
	ftree->jet_TCHE.push_back(jet.bDiscriminator("trackCountingHighEffBJetTags"));
	ftree->jet_SSVHE.push_back(jet.bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
	ftree->jet_SSVHP.push_back(jet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
	ftree->jet_CMVA.push_back(jet.bDiscriminator("combinedMVABJetTags"));
	
	ftree->jet_chargedMultiplicity.push_back(jet.chargedMultiplicity());
	ftree->jet_neutralMultiplicity.push_back(jet.neutralMultiplicity());
	ftree->jet_chargedHadronMultiplicity.push_back(jet.chargedHadronMultiplicity());
	
	ftree->jet_jecFactorUncorrected.push_back(jet.jecFactor("Uncorrected"));
	ftree->jet_jecFactorL1FastJet.push_back(jet.jecFactor("L1FastJet"));
	ftree->jet_jecFactorL2Relative.push_back(jet.jecFactor("L2Relative"));
	ftree->jet_jecFactorL3Absolute.push_back(jet.jecFactor("L3Absolute"));
	
	ftree->jet_ntrk.push_back(jet.associatedTracks().size());

	float CSVIVF = jet.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags");
	if( CSVIVF == -1000. ) CSVIVF = jet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
	ftree->jet_CSVv2.push_back(CSVIVF);

	float CSV = jet.bDiscriminator("combinedSecondaryVertexBJetTags");
	if( CSV == -1000. ) CSV = jet.bDiscriminator("pfCombinedSecondaryVertexBJetTags");

	ftree->jet_CSV.push_back(CSV);

	ftree->jet_flavour.push_back(jet.partonFlavour());

	ftree->jet_neutralHadronEnergy.push_back(jet.neutralHadronEnergy());
	ftree->jet_neutralEmEnergy.push_back(jet.neutralEmEnergy());
	ftree->jet_chargedHadronEnergy.push_back(jet.chargedHadronEnergy());
	ftree->jet_chargedEmEnergy.push_back(jet.chargedEmEnergy());
	ftree->jet_electronEnergy.push_back(jet.electronEnergy());
	ftree->jet_muonEnergy.push_back(jet.muonEnergy());
	ftree->jet_photonEnergy.push_back(jet.photonEnergy());

	ftree->jet_pileupJetId.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));

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

