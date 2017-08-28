#########################
#Author: Sam Higginbotham
########################
from WMCore.Configuration import Configuration
config = Configuration()


#name='Pt15to30'
config.section_("General")
config.General.requestName = 'PCC_AlCaLumiPixels_Run2017B_PIXONLY_RawPCC_Random_LS_v3'
config.General.workArea = 'RawPCCRandom2017B'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'raw_corr_Random_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['correctionsC.db','CorrectionHisto.root']


#config.JobType.inputFiles = ['dttf_config.db']

config.section_("Data")
config.Data.inputDataset = '/AlCaLumiPixels/Run2017B-AlCaPCCRandom-PromptReco-v1/ALCARECO'
#config.Data.lumiMask = ''
#config.Data.runRange='297283,297278,297280,297281,297271,297227,297230,297276,297265,297266'
config.Data.ignoreLocality = True
#useParent = True


config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.unitsPerJob = 500
#config.Data.totalUnits = -1
#config.Data.publishDbsUrl = 'test'
config.Data.outputDatasetTag = 'PCC_AlCaLumiPixels_Run2017B_RawPCC_Random_PIXONLY_LS'
config.Data.outLFNDirBase = '/store/group/comm_luminosity/PCC/ForLumiSystematics/2017/Aug24'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist=['T2_FR_CCIN2P3','T2_IT_Pisa','T2_UK_London_IC','T2_HU_Budapest']
#config.Site.whitelist=['T2_FR_CCIN2P3']
