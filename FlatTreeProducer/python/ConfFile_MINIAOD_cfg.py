import FWCore.ParameterSet.Config as cms

#####################
#  Options parsing  #
#####################

from FWCore.ParameterSet.VarParsing import VarParsing
import os, sys

options = VarParsing('analysis')
options.register('isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Run on real data')
options.register('applyJEC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply JEC corrections')
options.register('confFile', 'conf.xml', VarParsing.multiplicity.singleton, VarParsing.varType.string, "Flattree variables configuration")
options.register('bufferSize', 32000, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Buffer size for branches of the flat tree")
options.parseArguments()

##########################
#  Global configuration  #
##########################

process = cms.Process("FlatTree")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# remove verbose from patTrigger due to missing L1 prescales for some trigger paths
#process.MessageLogger.suppressWarning.append('patTrigger')
#process.MessageLogger.cerr.FwkJob.limit=1
#process.MessageLogger.cerr.ERROR = cms.untracked.PSet( limit = cms.untracked.int32(0) )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

if options.isData:
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
    from Configuration.AlCa.GlobalTag import GlobalTag
    process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1' # DATA 50ns
else:
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
    from Configuration.AlCa.GlobalTag import GlobalTag
    process.GlobalTag.globaltag = 'MCRUN2_74_V9A::All' # MC 50ns
    #process.GlobalTag.globaltag = 'MCRUN2_74_V9::All' # MC 25ns
    
########################
#  Additional modules  #
########################

if not options.isData:
    process.load('IPHCFlatTree.FlatTreeProducer.genJetFlavorMatching')

#####################
#  JES corrections  #
#####################

jetsName="slimmedJets"

if options.applyJEC:
    
    process.load("CondCore.DBCommon.CondDBCommon_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.jec = cms.ESSource("PoolDBESSource",
         DBParameters = cms.PSet(
           messageLevel = cms.untracked.int32(0)
           ),
         timetype = cms.string('runnumber'),
         toGet = cms.VPSet(
         cms.PSet(
               record = cms.string('JetCorrectionsRecord'),
               tag    = cms.string('JetCorrectorParametersCollection_Summer15_50nsV2_MC_AK4PFchs'),
               label  = cms.untracked.string('AK4PFchs')
               ),
         ## here you add as many jet types as you need
         ## note that the tag name is specific for the particular sqlite file
         ),
         connect = cms.string('sqlite:Summer15_50nsV2_MC.db')
         # uncomment above tag lines and this comment to use MC JEC
         # connect = cms.string('sqlite:Summer12_V7_MC.db')
    )
    ## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
    process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
    process.patJetCorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
      src = cms.InputTag("slimmedJets"),
      levels = ['L1FastJet',
            'L2Relative',
            'L3Absolute'],
            payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!
#    process.patJetCorrFactorsReapplyJEC.levels.append('L2L3Residual')

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
    process.patJetsReapplyJEC = process.patJetsUpdated.clone(
      jetSource = cms.InputTag("slimmedJets"),
      jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
      )
    process.JEC = cms.Sequence( process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC )

    jetsName="patJetsReapplyJEC"

# egamma
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)

# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_50ns_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

# MET
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
# for now only HBHE should be rerun on MINIAOD
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

###########
#  Input  #
###########

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"), # WARNING / FIXME for test only !
    fileNames = cms.untracked.vstring(
#    'root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/testFiles/MiniAOD/ttH_ev_2.root'
    'file:MC.root'
#    'file:DATA.root'
        #'file:/opt/sbg/data/safe1/cms/xcoubez/PhD/Analysis/WZAnalysisX/TriggerAgain/KirillFlatTreeStandalone/CMSSW_7_2_3_MantaRayXavier/CMSSW_7_2_3/src/IPHCFlatTree/FlatTreeProducer/test/InputRootFile/step2.root'
        #'/store/mc/Phys14DR/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/484D51C6-2673-E411-8AB0-001E67398412.root'
#        'root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/testFiles/MiniAOD/ttH_ev_2.root'
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
                  
                  vertexInput              = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  electronInput            = cms.InputTag("slimmedElectrons"),
                  electronPATInput         = cms.InputTag("slimmedElectrons"),
                  
                  eleVetoCBIdMap           = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V1-standalone-veto"),
                  eleLooseCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V1-standalone-loose"),
                  eleMediumCBIdMap         = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V1-standalone-medium"),
                  eleTightCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V1-standalone-tight"),
                  eleHEEPCBIdMap           = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),

                  eleMediumMVAIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                  eleTightMVAIdMap         = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                  mvaValuesMap             = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                  mvaCategoriesMap         = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),

                  filterTriggerNames       = cms.untracked.vstring("*"),
                  
                  muonInput                = cms.InputTag("slimmedMuons"),
                  tauInput                 = cms.InputTag("slimmedTaus"),
                  jetInput                 = cms.InputTag(jetsName),
                  jetPuppiInput            = cms.InputTag("slimmedJetsPuppi"),
                  genJetInput              = cms.InputTag("slimmedGenJets"),
                  jetFlavorMatchTokenInput = cms.InputTag("jetFlavourMatch"),
                  metInput                 = cms.InputTag("slimmedMETs"),
                  metPuppiInput            = cms.InputTag("slimmedMETsPuppi"),
                  rhoInput                 = cms.InputTag("fixedGridRhoFastjetAll"),
                  genParticlesInput        = cms.InputTag("prunedGenParticles"),
                  objects                  = cms.InputTag("selectedPatTrigger")
)

##########
#  Path  #
##########

if not options.isData:
    if options.applyJEC:
        process.p = cms.Path(
        process.HBHENoiseFilterResultProducer+
        process.ApplyBaselineHBHENoiseFilter+
        process.JEC+
        process.electronMVAValueMapProducer+
        process.egmGsfElectronIDSequence+
        process.genJetFlavourAlg+
        process.FlatTree)
    else:
        process.p = cms.Path(
        process.HBHENoiseFilterResultProducer+
        process.ApplyBaselineHBHENoiseFilter+
        process.electronMVAValueMapProducer+
        process.egmGsfElectronIDSequence+
        process.genJetFlavourAlg+
        process.FlatTree)        
else:
    if options.applyJEC:
        process.p = cms.Path(
        process.HBHENoiseFilterResultProducer+
        process.ApplyBaselineHBHENoiseFilter+
        process.JEC+
        process.electronMVAValueMapProducer+
        process.egmGsfElectronIDSequence+
        process.FlatTree)
    else:
        process.p = cms.Path(
        process.HBHENoiseFilterResultProducer+
        process.ApplyBaselineHBHENoiseFilter+
        process.electronMVAValueMapProducer+
        process.egmGsfElectronIDSequence+
        process.FlatTree)
    
    
