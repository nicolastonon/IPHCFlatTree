from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'REQUESTNAME'
config.section_('JobType')
config.JobType.psetName = '../runFlatTreeMINIAOD_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['../conf.xml','../Summer15_50nsV2_MC.db']
#config.JobType.outputFiles = ['output.root']
config.JobType.pyCfgParams = ['isData=0','applyJEC=1']
config.section_('Data')
config.Data.totalUnits = 1
#config.Data.totalUnits = 10
#config.Data.unitsPerJob = 2
config.Data.unitsPerJob = 20
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.inputDataset = 'INPUTDATASET'
config.Data.publishDataName = 'PUBLISHDATANAME'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
config.Data.outLFNDirBase = 'OUTLFN'
config.Data.lumiMask = 'GRL/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt' 
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_FR_IPHC'
