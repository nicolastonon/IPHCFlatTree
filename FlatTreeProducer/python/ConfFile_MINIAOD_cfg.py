import FWCore.ParameterSet.Config as cms

#####################
#  Options parsing  #
#####################

from FWCore.ParameterSet.VarParsing import VarParsing
import os, sys

options = VarParsing('analysis')
options.register('isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.int,'Run on real data')
options.register('confFile', 'conf.xml', VarParsing.multiplicity.singleton, VarParsing.varType.string, "Flattree variables configuration")
options.register('bufferSize', 32000, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Buffer size for branches of the flat tree")
options.parseArguments()

##########################
#  Global configuration  #
##########################

process = cms.Process("FlatTree")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

########################
#  Additional modules  #
########################

process.load('IPHCFlatTree.FlatTreeProducer.genJetFlavorMatching')

###########
#  Input  #
###########

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"), # WARNING / FIXME for test only !
    fileNames = cms.untracked.vstring(
        '/store/mc/Phys14DR/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/484D51C6-2673-E411-8AB0-001E67398412.root'
        #'root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/testFiles/MiniAOD/ttH_ev_2.root'
    )
)

############
#  Output  #
############

process.TFileService = cms.Service("TFileService", fileName = cms.string("output.root"))

#############################
#  Flat Tree configuration  #
#############################

process.FlatTree = cms.EDAnalyzer('FlatTreeProducer',

                  dataFormat        = cms.string("MINIAOD"),

                  bufferSize        = cms.int32(options.bufferSize),
                  confFile          = cms.string(options.confFile),

                  isData            = cms.bool(options.isData),

                  vertexInput       = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  electronInput     = cms.InputTag("slimmedElectrons"),
                  muonInput         = cms.InputTag("slimmedMuons"),
                  tauInput          = cms.InputTag("slimmedTaus"),
                  jetInput          = cms.InputTag("slimmedJets"),
                  genJetInput       = cms.InputTag("slimmedGenJets"),
                  jetFlavorMatchTokenInput = cms.InputTag("jetFlavourMatch"),
                  metInput          = cms.InputTag("slimmedMETs"),
                  rhoInput          = cms.InputTag("fixedGridRhoFastjetAll"),
                  genParticlesInput = cms.InputTag("prunedGenParticles"),
                  objects           = cms.InputTag("selectedPatTrigger")
)

##########
#  Path  #
##########

process.p = cms.Path(process.genJetFlavourAlg*
                     process.FlatTree)
