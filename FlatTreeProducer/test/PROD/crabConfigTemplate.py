from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'REQUESTNAME'
config.section_('JobType')
config.JobType.psetName = '../runFlatTreeMINIAOD_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['../conf.xml', '../Summer15_25nsV7_MC_L2Relative_AK8PFchs.txt', '../Summer15_25nsV7_MC_L3Absolute_AK8PFchs.txt', '../Summer15_25nsV7_DATA_L2L3Residual_AK8PFchs.txt']
#config.JobType.outputFiles = ['output.root']
config.JobType.pyCfgParams = ['isData=0','runAK10=1','runQG=1']
config.section_('Data')
config.Data.totalUnits = -1 #@MJ@ TODO
#config.Data.totalUnits = 10
config.Data.unitsPerJob = 9000
#config.Data.unitsPerJob = 20
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.publication = False
config.Data.inputDataset = 'INPUTDATASET'
#config.Data.inputDBS = 'phys03'
config.Data.outputDatasetTag = 'PUBLISHDATANAME'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
config.Data.outLFNDirBase = 'OUTLFN'
#config.Data.lumiMask = 'GRL/Cert_246908-258159_13TeV_PromptReco_Collisions15_25ns_JSON_v3.txt'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_FR_IPHC'
