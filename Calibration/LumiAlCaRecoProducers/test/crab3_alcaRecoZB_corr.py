#########################
#Author: Sam Higginbotham
########################
from WMCore.Configuration import Configuration
config = Configuration()


config.section_("General")
config.General.requestName = 'PCC_Run2018B_Corrections'
config.General.workArea = 'AlCaRecoZeroBias2018'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'lumi_alcaZB_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['rawPCC.csv']

config.section_("Data")
config.Data.inputDataset = '/AlCaLumiPixels/Run2018B-AlCaPCCZeroBias-PromptReco-v1/ALCARECO'
#config.Data.runRange=''
config.Data.ignoreLocality = True
#useParent = True


config.Data.inputDBS = 'global'
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'FileBased'
config.Data.publication = False
#config.Data.unitsPerJob = 1
config.Data.outputDatasetTag = 'PCC_AlCaLumiPixels_Run2018C_1kLS_NoZeroes'
config.Data.outLFNDirBase = '/store/group/comm_luminosity/PCC/ForLumiComputations/2018/'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
