import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
import os

options = VarParsing('analysis')
options.register('isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.int,'Run on real data')

process = cms.Process("FlatTree")
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#    'root://sbgse1.in2p3.fr//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/20000/14587980-CB7E-E411-A0F4-001E67397701.root'
#'root://xrootd.unl.edu//store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/007B37D4-8B70-E411-BC2D-0025905A6066.root'
    'root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/testFiles/MiniAOD/ttH_ev_2.root'
#'root://xrootd.unl.edu//store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2010CF0F-CB71-E411-9331-002481E0D678.root'
#        'root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/testFiles/TZqSynchMoreStat/001C0B68-536A-E311-B25F-002590D0B066.root'
    )
)

#process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
#from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
#process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')
#process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
#from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
#process.ak4PFJetsCHS = ak4PFJets.clone(src = 'pfCHS', doAreaFastjet = True)

#bTagDiscriminators = ['pfTrackCountingHighEffBJetTags',
#'pfTrackCountingHighPurBJetTags',
#'pfJetProbabilityBJetTags',
#'pfJetBProbabilityBJetTags',
#'pfSimpleSecondaryVertexHighEffBJetTags',
#'pfSimpleSecondaryVertexHighPurBJetTags',
#'pfCombinedSecondaryVertexBJetTags',
#'pfCombinedInclusiveSecondaryVertexV2BJetTags'
#]

#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJetCollection(
#  process,
#  jetSource = cms.InputTag('ak4PFJetsCHS'),
#  pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
#  pfCandidates = cms.InputTag('packedPFCandidates'),
#  svSource = cms.InputTag('slimmedSecondaryVertices'),
#  btagDiscriminators = bTagDiscriminators,
#  jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
#  genJetCollection = cms.InputTag('ak4GenJetsNoNu')
#)

#getattr(process,'patJetPartons').particles = cms.InputTag('prunedGenParticles')
#getattr(process,'patJetPartonMatch').matched = cms.InputTag('prunedGenParticles')
#if hasattr(process,'pfInclusiveSecondaryVertexFinderTagInfos'):
#    getattr(process,'pfInclusiveSecondaryVertexFinderTagInfos').extSVCollection = cms.InputTag('slimmedSecondaryVertices')
#getattr(process,'patJets').addAssociatedTracks = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
#getattr(process,'patJets').addJetCharge = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
#getattr(process,'selectedPatJets').cut = cms.string('pt > 10') # to match the selection for slimmedJets in MiniAOD
#from PhysicsTools.PatAlgos.tools.pfTools import *
#adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

process.FlatTree = cms.EDAnalyzer('FlatTreeProducer',
                                  dataFormat = cms.string("MINIAOD"),
                                  isData = cms.bool(options.isData),
                                  vertexInput = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                  electronInput = cms.InputTag("slimmedElectrons"),
                                  muonInput = cms.InputTag("slimmedMuons"),
#                                  jetInput = cms.InputTag("selectedPatJets"),
                                  jetInput = cms.InputTag("slimmedJets"),
                                  metInput = cms.InputTag("slimmedMETs"),
                                  rhoInput = cms.InputTag("fixedGridRhoFastjetAll"),
#                                  rhoInput = cms.InputTag("fixedGridRhoAll"),
                                  genParticlesInput = cms.InputTag("prunedGenParticles")
)


process.TFileService = cms.Service("TFileService",
fileName = cms.string("output.root")
#fileName = cms.string(options.outputFile)
)

process.p = cms.Path(
                     process.FlatTree)
