import FWCore.ParameterSet.Config as cms

#####################
#  Options parsing  #
#####################

from FWCore.ParameterSet.VarParsing import VarParsing
import os, sys

options = VarParsing('analysis')
options.register('isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Run on real data')
options.register('applyMETFilters',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply MET filters')
options.register('applyJEC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply JEC corrections')
# runBTag option is not fully functional - please don't use it
options.register('runBTag',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Run b-tagging')
options.register('runAK10',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Add AK10 jets')
options.register('runQG',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Run QGTagger')
options.register('fillMCScaleWeight',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Fill PDF weights')
options.register('fillPUInfo',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Fill PU info')
options.register('nPDF', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "nPDF")
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

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

if options.isData:
    process.GlobalTag.globaltag = '76X_dataRun2_v15'
else:
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'

### to activate the new JP calibration: using the data base
if options.isData:
    trkProbaCalibTag = "JPcalib_Data76X_2015D_v1"
else:
    trkProbaCalibTag = "JPcalib_MC76X_25ns_v1"
    
process.GlobalTag.toGet = cms.VPSet(
cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
         tag = cms.string(trkProbaCalibTag),
         connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
         )
)

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

# Re-apply JEC to AK4
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
process.patJetCorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet','L2Relative','L3Absolute'],
    payload = 'AK4PFchs' 
)

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = process.patJetsUpdated.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
)

# Re-apply JEC to AK8
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
process.patJetCorrFactorsReapplyJECAK8 = process.patJetCorrFactorsUpdated.clone(
    src = cms.InputTag("slimmedJetsAK8"),
    levels = ['L1FastJet','L2Relative','L3Absolute'],
    payload = 'AK8PFchs' 
)

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJECAK8 = process.patJetsUpdated.clone(
    jetSource = cms.InputTag("slimmedJetsAK8"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECAK8"))
)

jetsNameAK4="patJetsReapplyJEC"
jetsNameAK8="patJetsReapplyJECAK8"
#jetsNameAK10="patJetsReapplyJECAK10"
jetsNameAK10="selectedPatJetsAK10PFCHS"

if options.runBTag:
    bTagDiscriminators = [
        'pfTrackCountingHighEffBJetTags',
        'pfTrackCountingHighPurBJetTags',
        'pfJetProbabilityBJetTags',
        'pfJetBProbabilityBJetTags',
        'pfSimpleSecondaryVertexHighEffBJetTags',
        'pfSimpleSecondaryVertexHighPurBJetTags',
        'pfCombinedSecondaryVertexV2BJetTags',
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfCombinedMVABJetTags'
    ]

    bTagInfos = [
        'pfImpactParameterTagInfos',
        'pfSecondaryVertexTagInfos',
        'pfInclusiveSecondaryVertexFinderTagInfos',
        'softPFMuonsTagInfos',
        'softPFElectronsTagInfos'
    ]

    ##########################
    #  Rebuild ak4PFJetsCHS  #
    ##########################

    # https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchool2015BTaggingExercise#Part_II_Remaking_b_tag_discrimin
    process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
    from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
    process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')

    process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    process.ak4PFJetsCHS = ak4PFJets.clone(src = 'pfCHS', doAreaFastjet = True)
    
    from PhysicsTools.PatAlgos.tools.jetTools import *
    switchJetCollection(
        process,
        jetSource = cms.InputTag('ak4PFJetsCHS'),
        pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
        genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
        genParticles = cms.InputTag('prunedGenParticles'),
        btagDiscriminators = bTagDiscriminators,
        btagInfos = bTagInfos
    )
    getattr(process,'selectedPatJets').cut = cms.string('pt > 10') # to match the selection for slimmedJets in MiniAOD
    process.selectedPatJetsAK4 = process.selectedPatJets.clone(src = 'selectedPatJets')

    from PhysicsTools.PatAlgos.tools.pfTools import *
    adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

    ##########################
    #  Rebuild ak8PFJetsCHS  #
    ##########################

    # https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchool2015BTaggingExercise#Subjet_b_tagging_and_subjet_flav
    process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))

    from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
    process.ak8GenJetsNoNu = ak4GenJets.clone(
        rParam = cms.double(0.8),
        src = cms.InputTag("packedGenParticlesForJetsNoNu")
    )

    from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
    from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
    process.ak8GenJetsNoNuPruned = ak4GenJets.clone(
        SubJetParameters,
        rParam = cms.double(0.8),
        src = cms.InputTag("packedGenParticlesForJetsNoNu"),
        usePruning = cms.bool(True),
        writeCompound = cms.bool(True),
        jetCollInstanceName=cms.string("SubJets")
    )

    process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    process.ak8PFJetsCHS = ak4PFJets.clone(
        rParam = cms.double(0.8),
        src = cms.InputTag("pfCHS"),
        doAreaFastjet = cms.bool(True),
        jetPtMin = cms.double(50.)
    )

    from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
    process.ak8PFJetsCHSPruned = ak5PFJetsPruned.clone(
        rParam = cms.double(0.8),
        src = cms.InputTag("pfCHS"),
        doAreaFastjet = cms.bool(True),
        writeCompound = cms.bool(True),
        jetCollInstanceName=cms.string("SubJets"),
        jetPtMin = cms.double(50.)
    )

    from PhysicsTools.PatAlgos.tools.jetTools import *

    addJetCollection(
        process,
        labelName = 'AK8PFCHS',
        jetSource = cms.InputTag('ak8PFJetsCHS'),
        algo = 'ak',  # needed for jet flavor clustering
        rParam = 0.8, # needed for jet flavor clustering
        pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagDiscriminators = bTagDiscriminators,
        btagInfos = bTagInfos,
        jetCorrections = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
        genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
        genParticles = cms.InputTag('prunedGenParticles')
    )

    addJetCollection(
        process,
        labelName = 'AK8PFCHSPruned',
        jetSource = cms.InputTag('ak8PFJetsCHSPruned'),
        btagDiscriminators = ['None'],
        jetCorrections = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
        genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
        genParticles = cms.InputTag('prunedGenParticles'),
        getJetMCFlavour = False # jet flavor disabled
    )

    addJetCollection(
        process,
        labelName = 'AK8PFCHSPrunedSubjets',
        jetSource = cms.InputTag('ak8PFJetsCHSPruned','SubJets'),
        algo = 'ak',  # needed for subjet flavor clustering
        rParam = 0.8, # needed for subjet flavor clustering
        pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagDiscriminators = bTagDiscriminators,
        btagInfos = bTagInfos,
        jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
        genJetCollection = cms.InputTag('ak8GenJetsNoNuPruned','SubJets'),
        genParticles = cms.InputTag('prunedGenParticles'),
        explicitJTA = True,  # needed for subjet b tagging
        svClustering = True, # needed for subjet b tagging
        fatJets=cms.InputTag('ak8PFJetsCHS'),             # needed for subjet flavor clustering
        groomedFatJets=cms.InputTag('ak8PFJetsCHSPruned') # needed for subjet flavor clustering
    )

    process.selectedPatJetsAK8PFCHSPrunedPacked = cms.EDProducer("BoostedJetMerger",
        jetSrc=cms.InputTag("selectedPatJetsAK8PFCHSPruned"),    
        subjetSrc=cms.InputTag("selectedPatJetsAK8PFCHSPrunedSubjets")
    )

    process.packedPatJetsAK8PFCHS = cms.EDProducer("JetSubstructurePacker",
        jetSrc = cms.InputTag('selectedPatJetsAK8PFCHS'),
        distMax = cms.double(0.8),
        algoTags = cms.VInputTag( cms.InputTag('selectedPatJetsAK8PFCHSPrunedPacked') ),
        algoLabels = cms.vstring( 'Pruned' ),
        fixDaughters = cms.bool(False)
    )

    process.selectedPatJetsAK8 = process.selectedPatJetsAK8PFCHS.clone(src = 'selectedPatJetsAK8PFCHS')

    from PhysicsTools.PatAlgos.tools.pfTools import *
    adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))
    
    jetsNameAK4="selectedPatJetsAK4"
    jetsNameAK8="selectedPatJetsAK8"
    jetsNameAK10="selectedPatJetsAK10"

########################
#  Additional modules  #
########################

#if not options.isData:
#    if not options.runBTag:
#        process.load('IPHCFlatTree.FlatTreeProducer.genJetFlavorMatching')

# egamma
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)

# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
my_id_modules = [
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_50ns_V2_cff',
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff'
]

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

# MET
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
    reverseDecision = cms.bool(False)
)

process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
    reverseDecision = cms.bool(False)
)

#####################
# MET Significance  #
#####################
process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")
from RecoMET.METProducers.testInputFiles_cff import recoMETtestInputFiles


#######################
# AK10 collection     #
#######################
if options.runAK10:
    from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
    jetToolbox( process, 'ak10', 'ak10JetSubs', 'out', runOnMC=(not options.isData),
                addPruning=True, addSoftDrop=True , addPrunedSubjets=True, addSoftDropSubjets=True,
                JETCorrPayload='AK3Pachs', subJETCorrPayload='AK10PFchs', JETCorrLevels=['L1FastJet', 'L2Relative', 'L3Absolute'],
                addNsub=True, maxTau=6, addTrimming=True, addFiltering=True,
                addEnergyCorrFunc=True, maxECF=5 )    
                
#######################
# Quark gluon tagging #
#######################
if options.runQG:
    qgDatabaseVersion = 'v1' # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

    from CondCore.DBCommon.CondDBSetup_cfi import *
    QGPoolDBESSource = cms.ESSource("PoolDBESSource",
          CondDBSetup,
          toGet = cms.VPSet(),
          connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
    )

    for type in ['AK4PFchs','AK4PFchs_antib']:
        QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
              record = cms.string('QGLikelihoodRcd'),
              tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
              label  = cms.untracked.string('QGL_'+type)
        )))

    process.load('RecoJets.JetProducers.QGTagger_cfi')
    process.QGTagger.srcJets          = cms.InputTag(jetsNameAK4) # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
    process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs') # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

###########
#  Input  #
###########

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"), # WARNING / FIXME for test only !
    fileNames = cms.untracked.vstring(
    'root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/testFiles/MiniAOD/ttH_76X.root'
#    '/store/data/Run2015D/SingleMuon/MINIAOD/05Oct2015-v1/40000/9A470821-676F-E511-8AF3-0025905A606A.root'
#    '/store/mc/RunIISpring15MiniAODv2/ttHToNonbb_M125_13TeV_powheg_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/02FE2DB6-D06D-E511-8BC7-0025905C431C.root'
            )
)

############
#  Output  #
############

process.TFileService = cms.Service("TFileService", fileName = cms.string("output.root"))

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True)	 # needed for ak10 computation (JMEAnalysis/JetToolbox)
)
        
#############################
#  Flat Tree configuration  #
#############################

process.FlatTree = cms.EDAnalyzer('FlatTreeProducer',

                  dataFormat        = cms.string("MINIAOD"),

                  bufferSize        = cms.int32(options.bufferSize),
                  confFile          = cms.string(options.confFile),

                  isData            = cms.bool(options.isData),
                  applyMETFilters   = cms.bool(options.applyMETFilters),
                  fillMCScaleWeight = cms.bool(options.fillMCScaleWeight),
                  fillPUInfo	    = cms.bool(options.fillPUInfo),
                  nPDF              = cms.int32(options.nPDF),
                  
                  vertexInput              = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  electronInput            = cms.InputTag("slimmedElectrons"),
                  electronPATInput         = cms.InputTag("slimmedElectrons"),

                  eleVetoCBIdMap           = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                  eleLooseCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                  eleMediumCBIdMap         = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                  eleTightCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                  eleHEEPCBIdMap           = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),

                  eleMediumMVAIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                  eleTightMVAIdMap         = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                  mvaValuesMap             = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                  mvaCategoriesMap         = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),

                  filterTriggerNames       = cms.untracked.vstring(
                  "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*",
                  "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*",
                  "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
                  "HLT_DoubleMu33NoFiltersNoVtx_v*",
                  "HLT_DoubleMu38NoFiltersNoVtx_v*",
                  "HLT_Ele22_eta2p1_WPLoose_Gsf_v*",
                  "HLT_Ele22_eta2p1_WPTight_Gsf_v*",
                  "HLT_Ele22_eta2p1_WP75_Gsf_v*",
                  "HLT_Ele23_WPLoose_Gsf_v*",
                  "HLT_Ele23_WP75_Gsf_v*",
                  "HLT_Ele27_eta2p1_WP75_Gsf_v*",
                  "HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v*",
                  "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
                  "HLT_Ele27_eta2p1_WPTight_Gsf_v*",
                  "HLT_Ele32_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v*",
                  "HLT_Ele32_eta2p1_WPLoose_Gsf_v*",
                  "HLT_Ele32_eta2p1_WPTight_Gsf_v*",
                  "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
                  "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*",
                  "HLT_IsoMu17_eta2p1_v*",
                  "HLT_DoubleIsoMu17_eta2p1_v*",
                  "HLT_IsoMu18_v*",
                  "HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v*",
                  "HLT_IsoMu20_v*",
                  "HLT_IsoMu20_eta2p1_v*",
                  "HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v*",
                  "HLT_IsoMu24_eta2p1_v*",
                  "HLT_IsoMu27_v*",
                  "HLT_IsoTkMu20_v*",
                  "HLT_IsoTkMu20_eta2p1_v*",
                  "HLT_IsoTkMu24_eta2p1_v*",
                  "HLT_IsoTkMu27_v*",
                  "HLT_Mu17_Mu8_v*",
                  "HLT_Mu17_Mu8_DZ_v*",
                  "HLT_Mu17_Mu8_SameSign_DZ_v*",
                  "HLT_Mu20_Mu10_v*",
                  "HLT_Mu20_Mu10_DZ_v*",
                  "HLT_Mu20_Mu10_SameSign_DZ_v*",
                  "HLT_Mu17_TkMu8_DZ_v*",
                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
                  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
                  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
                  "HLT_Mu27_TkMu8_v*",
                  "HLT_Mu30_TkMu11_v*",
                  "HLT_Mu40_TkMu11_v*",
                  "HLT_Mu20_v*",
                  "HLT_TkMu20_v*",
                  "HLT_Mu24_eta2p1_v*",
                  "HLT_TkMu24_eta2p1_v*",
                  "HLT_Mu27_v*",
                  "HLT_TkMu27_v*",
                  "HLT_Mu50_v*",
                  "HLT_Mu55_v*",
                  "HLT_Mu45_eta2p1_v*",
                  "HLT_Mu50_eta2p1_v*",
                  "HLT_PFJet40_v*",
                  "HLT_PFJet60_v*",
                  "HLT_PFJet80_v*",
                  "HLT_PFJet140_v*",
                  "HLT_PFJet200_v*",
                  "HLT_PFJet260_v*",
                  "HLT_PFJet320_v*",
                  "HLT_PFJet400_v*",
                  "HLT_PFJet450_v*",
                  "HLT_PFJet500_v*",
                  "HLT_Mu8_TrkIsoVVL_v*",
                  "HLT_Mu17_TrkIsoVVL_v*",
                  "HLT_Mu24_TrkIsoVVL_v*",
                  "HLT_Mu34_TrkIsoVVL_v*",
                  "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
                  "HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
                  "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
                  "HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
                  "HLT_BTagMu_DiJet20_Mu5_v*",
                  "HLT_BTagMu_DiJet40_Mu5_v*",
                  "HLT_BTagMu_DiJet70_Mu5_v*",
                  "HLT_BTagMu_DiJet110_Mu5_v*",
                  "HLT_BTagMu_Jet300_Mu5_v*",
                  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                  "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                  "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*",
                  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*",
                  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Mu8_v*",
                  "HLT_Mu17_v*",
                  "HLT_Mu24_v*",
                  "HLT_Mu34_v*",
                  "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*",
                  "HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v*",
                  "HLT_Ele18_CaloIdM_TrackIdM_PFJet30_v*",
                  "HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v*",
                  "HLT_Ele33_CaloIdM_TrackIdM_PFJet30_v*",
                  "HLT_Mu300_v*",
                  "HLT_Mu350_v*",
                  "HLT_PFHT450_SixJet40_PFBTagCSV0p72_v*",
                  "HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v*",
                  "HLT_PFHT450_SixJet40_PFBTagCSV_v*",
                  "HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v*"
                  ),
                  
                  muonInput                = cms.InputTag("slimmedMuons"),
                  tauInput                 = cms.InputTag("slimmedTaus"),
                  jetInput                 = cms.InputTag(jetsNameAK4),
                  jetPuppiInput            = cms.InputTag("slimmedJetsPuppi"),
                  ak8jetInput              = cms.InputTag(jetsNameAK8),
                  ak10jetInput             = cms.InputTag(jetsNameAK10),
                  genJetInput              = cms.InputTag("slimmedGenJets"),
                  jetFlavorMatchTokenInput = cms.InputTag("jetFlavourMatch"),
                  metInput                 = cms.InputTag("slimmedMETs"),
                  metPuppiInput            = cms.InputTag("slimmedMETsPuppi"),
                  metNoHFInput             = cms.InputTag("slimmedMETsNoHF"),
                  metSigInput              = cms.InputTag("METSignificance"),
                  metCovInput              = cms.InputTag("METSignificance","METCovariance"),
                  rhoInput                 = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                  genParticlesInput        = cms.InputTag("prunedGenParticles"),
                  genEventInfoInput        = cms.InputTag("generator"),
                  LHEEventProductInput     = cms.InputTag("externalLHEProducer"),
                  bsInput                  = cms.InputTag("offlineBeamSpot"),
                  pfcandsInput             = cms.InputTag("packedPFCandidates"),
                  hConversionsInput        = cms.InputTag("reducedEgamma","reducedConversions"),
                  puInfoInput		   = cms.InputTag("slimmedAddPileupInfo"),
#                  puInfoInput		   = cms.InputTag("addPileupInfo"),
                  objects                  = cms.InputTag("selectedPatTrigger")
)

##########
#  Path  #
##########

#process.metFilters = cms.Sequence()
#if options.applyMETFilters:
#    process.metFilters = cms.Sequence(process.ApplyBaselineHBHENoiseFilter+
#                                      process.ApplyBaselineHBHEIsoNoiseFilter)

process.runQG = cms.Sequence()
if options.runQG:
    process.runQG = cms.Sequence(process.QGTagger)

process.runBTag = cms.Sequence()
if options.runBTag:
    process.runBTag = cms.Sequence(process.selectedPatJetsAK8PFCHS+
                                   process.selectedPatJetsAK8PFCHSPrunedPacked)
                                   
process.p = cms.Path(
#                     process.metFilters+
                     process.electronMVAValueMapProducer+
                     process.egmGsfElectronIDSequence+
                     process.patJetCorrFactorsReapplyJEC+process.patJetsReapplyJEC+
                     process.METSignificance+
                     process.runBTag+
                     process.runQG+
                     process.FlatTree
                    )
