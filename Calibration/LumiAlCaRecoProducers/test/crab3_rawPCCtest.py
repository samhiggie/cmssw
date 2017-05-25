#########################
#Author: Sam Higginbotham
'''

* File Name :

* Purpose :

* Creation Date : 10-05-2017

* Last Modified :

'''
#########################
from WMCore.Configuration import Configuration
config = Configuration()


#name='Pt15to30'
config.section_("General")
config.General.requestName = 'PCC_ZeroBias_DataCert_150820'
config.General.workArea = 'taskManagement'

config.section_("JobType")
config.JobType.pluginName = 'ANALYSIS'
config.JobType.psetName = 'alcaPCCZeroBias_ALCA.py'
config.JobType.allowUndistributedCMSSW = True


#config.JobType.inputFiles = ['dttf_config.db']

config.section_("Data")
config.Data.inputDataset = '/AlCaLumiPixels1/Run2016H-v1/RAW'
config.Data.lumiMask = ''

config.Data.ignoreLocality = True
#useParent = True


config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.publication = False
config.Data.unitsPerJob = 10
#config.Data.totalUnits = -1
#config.Data.publishDbsUrl = 'test'
config.Data.outputDatasetTag = 'PCC_ZeroBias_DataCert_150820'
config.Data.outLFNDirBase = '/store/group/comm_luminosity/PCC/ALCA'
config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist=['T2_FR_CCIN2P3','T2_IT_Pisa','T2_UK_London_IC','T2_HU_Budapest']
#config.Site.whitelist=['T2_FR_CCIN2P3']



