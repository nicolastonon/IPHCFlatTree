import FWCore.ParameterSet.Config as cms

#####################
#  Options parsing  #
#####################

from FWCore.ParameterSet.VarParsing import VarParsing
import os

options = VarParsing('analysis')
options.register('isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.int,'Run on real data')

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

###########
#  Input  #
###########

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/testFiles/MiniAOD/ttH_ev_2.root'
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

    isData            = cms.bool(options.isData),
    vertexInput       = cms.InputTag("offlineSlimmedPrimaryVertices"),
    electronInput     = cms.InputTag("slimmedElectrons"),
    muonInput         = cms.InputTag("slimmedMuons"),
    jetInput          = cms.InputTag("slimmedJets"),
    metInput          = cms.InputTag("slimmedMETs"),
    rhoInput          = cms.InputTag("fixedGridRhoFastjetAll"),
    genParticlesInput = cms.InputTag("prunedGenParticles")
)

##########
#  Path  #
##########

process.p = cms.Path(process.FlatTree)
