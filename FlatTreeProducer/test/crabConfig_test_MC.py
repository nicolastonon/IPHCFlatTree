from IPHCFlatTree.FlatTreeProducer.ConfFile_MINIAOD_cfg import *
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'THQ_Hincl_13TeV_madgraph_pythia8_TuneCUETP8M1'
config.section_('JobType')
config.JobType.psetName = 'runFlatTreeMINIAOD_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['conf.xml','Summer16_23Sep2016V3_MC.db','Summer16_23Sep2016AllV3_DATA.db']
#config.JobType.outputFiles = ['output.root']
config.JobType.pyCfgParams = ['isData=0','runAK10=0']
config.section_('Data')
config.Data.totalUnits = -1
#config.Data.totalUnits = 10
#config.Data.unitsPerJob = 2
config.Data.unitsPerJob = 30
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.inputDataset = '/THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDBS = 'phys03'
config.Data.outputDatasetTag = 'RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
config.Data.outLFNDirBase = '/store/user/ntonon/FlatTree/Walrus-patch2/'
config.Data.lumiMask = 'PROD/GRL/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_FR_IPHC'
