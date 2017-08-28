#########################
#Author: Sam Higginbotham
'''

* File Name : crab3_raw_lumi.py

* Purpose :to run the raw_lumi_ZeroBias_cfg.py on the sets of data that are AlCaPCCZeroBias. 

* Creation Date : 28-08-2017

* Last Modified :

'''
#########################

from WMCore.Configuration import Configuration
config = Configuration()


#name='Pt15to30'
config.section_("General")
config.General.requestName = 'PCC_AlCaLumiPixels_Run2017B_PIXONLY_RawPCC_ZeroBias_LS_v3'
config.General.workArea = 'RawPCCZerobias2017BTEST'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'raw_lumi_ZeroBias_cfg.py'
config.JobType.allowUndistributedCMSSW = True
#For the correctionsC.db files needed in the lumiproducer 
config.JobType.inputFiles = ['/eos/cms/store/group/comm_luminosity/PCC/ForLumiSystematics/2017/Aug24/AlCaLumiPixels/PCC_AlCaLumiPixels_Run2017B_RawPCC_Random_PIXONLY_LS/170826_064251/0000/']
config.JobType.outputFiles = ['lumi.csv']


#config.JobType.inputFiles = ['dttf_config.db']

config.section_("Data")
config.Data.inputDataset = '/AlCaLumiPixels/Run2017B-AlCaPCCZeroBias-PromptReco-v1/ALCARECO'
#config.Data.lumiMask = ''
config.Data.runRange='297411'#'297283,297278,297280,297281,297271,297227,297230,297276,297265,297266'
config.Data.ignoreLocality = True
#useParent = True


config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.unitsPerJob = 500
#config.Data.totalUnits = -1
#config.Data.publishDbsUrl = 'test'
config.Data.outputDatasetTag = 'PCC_AlCaLumiPixels_Run2017B_RawPCC_ZeroBias_PIXONLY_LS_TEST'
config.Data.outLFNDirBase = '/store/group/comm_luminosity/PCC/ForLumiComputation/2017/Aug24'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist=['T2_FR_CCIN2P3','T2_IT_Pisa','T2_UK_London_IC','T2_HU_Budapest']
#config.Site.whitelist=['T2_FR_CCIN2P3']
