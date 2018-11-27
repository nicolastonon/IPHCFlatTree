from IPHCFlatTree.FlatTreeProducer.ConfFile_MINIAOD_cfg import *
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'test'
config.section_('JobType')

config.JobType.psetName = 'runFlatTreeMINIAOD_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['conf.xml','Fall17_17Nov2017_V6_MC.db','Fall17_17Nov2017BCDEF_V6_DATA.db']
#config.JobType.outputFiles = ['output.root']
config.JobType.pyCfgParams = ['isData=0','runAK10=0']
config.section_('Data')

config.Data.totalUnits = -1 #nof files (or lumisection) to analyze in total
#config.Data.totalUnits = 10

config.Data.unitsPerJob = 1

config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'

config.Data.publication = False

#NB : when running interactively, files will be taken from ../python/ConfFile_MINIAOD_cfg.py, in the 'Input' section ! Add filenames there !
config.Data.inputDataset = '/THW_5f_Hincl_13TeV_madgraph_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'

#config.Data.inputDBS = 'phys03'

config.Data.outputDatasetTag = 'test'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
config.Data.outLFNDirBase = '/store/user/ntonon/FlatTree/tHq2017'

#config.Data.lumiMask = '' #FIXME -- FIND 2017 JSON!

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_FR_IPHC'
