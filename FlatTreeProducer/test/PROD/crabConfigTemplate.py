from WMCore.Configuration import Configuration

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'REQUESTNAME'
config.section_('JobType')

config.JobType.psetName = '../runFlatTreeMINIAOD_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['../conf.xml','../Summer16_23Sep2016V3_MC.db','../Summer16_23Sep2016AllV3_DATA.db']
#config.JobType.outputFiles = ['output.root']
config.JobType.pyCfgParams = ['isData=1','runAK10=0']
config.section_('Data')

config.Data.totalUnits = -1 #nof files (or lumisection) to analyze in total
#config.Data.totalUnits = 10

#config.Data.unitsPerJob = 1 #nof files (or lumisections) in each job
config.Data.unitsPerJob = 10

#config.Data.splitting = 'FileBased'
config.Data.splitting = 'LumiBased'

config.Data.publication = False

config.Data.inputDataset = 'INPUTDATASET'
#config.Data.inputDBS = 'phys03'

config.Data.outputDatasetTag = 'PUBLISHDATANAME'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
config.Data.outLFNDirBase = 'OUTLFN'

config.Data.lumiMask = 'GRL/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_FR_IPHC'

#config.Data.ignoreLocality = True
