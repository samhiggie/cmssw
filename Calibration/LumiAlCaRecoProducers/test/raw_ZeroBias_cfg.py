#########################
#Author: Sam Higginbotham
#Purpose: To investigate the AlCaPCCProducer input and output. 
#########################
import FWCore.ParameterSet.Config as cms

process = cms.Process("rawRECO")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_1.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_10.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_100.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_101.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_102.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_103.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_104.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_105.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_106.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_107.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_108.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_109.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_11.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_110.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_111.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_112.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_113.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_114.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_115.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_116.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_117.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_118.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_119.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_12.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_120.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_121.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_122.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_123.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_124.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_125.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_126.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_127.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_128.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_129.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_13.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_130.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_131.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_132.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_133.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_134.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_135.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_136.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_137.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_138.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_139.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_14.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_140.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_141.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_142.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_143.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_144.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_145.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_146.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_147.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_148.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_149.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_15.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_150.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_151.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_152.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_153.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_154.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_155.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_156.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_157.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_158.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_159.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_16.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_160.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_161.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_162.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_163.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_164.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_165.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_166.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_167.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_168.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_169.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_17.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_170.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_171.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_172.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_173.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_174.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_175.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_176.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_177.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_178.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_179.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_18.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_180.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_181.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_182.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_183.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_184.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_185.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_186.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_187.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_188.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_189.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_19.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_190.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_191.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_192.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_193.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_194.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_195.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_196.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_197.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_198.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_199.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_2.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_20.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_200.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_201.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_202.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_203.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_204.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_205.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_206.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_207.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_208.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_209.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_21.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_210.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_211.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_212.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_213.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_214.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_215.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_216.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_217.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_218.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_219.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_22.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_220.root','file:/eos/cms/store/group/comm_luminosity/PCC/ALCA/AlCaLumiPixels1/PCC_ZeroBias_DataCert_150820/170511_182353/0000/AlCaPCCZeroBias_221.root')
)
#Added process to select the appropriate events 
process.OutALCARECOPromptCalibProdPCC = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOPromptCalibProdPCC')
    ),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_rawPCCProd_*_*')
)

#Make sure that variables match in producer.cc and .h
process.rawPCCProd = cms.EDProducer("RawPCCProducer",
    RawPCCProducerParameters = cms.PSet(
        #Mod factor to count lumi and the string to specify output 
        PCCobLabel = cms.string("alcaPCCProducerZeroBias"),
        ProdInst = cms.string("alcaPCCZeroBias"),
        resetEveryNLumi = cms.untracked.int32(1),
        trigstring = cms.untracked.string("rawPCCtest"), 
        #Below is a list of module IDs that will be ignored in calculation of luminosity
        modVeto=cms.vint32(302057476,302057496,302121488,302121492,302121732,302122248,302122500,302122768,302122780,302122784,302123024,302123036,302123292,302123296,302123780,302123800,302123804,302123808,302124036,302124040,302124292,302124296,302124300,302124304,302124308,302124548,302124552,302124556,302124560,302124812,302125848,302126084,302126344,302126364,302127136,302127876,302127888,302128156,302128388,302128392,302128416,302128644,302186768,302186776,302188552,302188820,302189332,302190868,302191112,302191388,302191648,302192412,302193164,302193424,302194184,302194192,302194200,302194440,302194716,302194968,302194972,302195220,302195224,302195228,302195232,302195476,302195480,302195484,302195488,302195732,302195736,302195740,302195744,302196000,302196996,302197536,344012048,344013072,344014096,344021252,344021256,344021260,344021264,344021508,344021512,344021516,344022276,344022280,344022284,344022288,344022532,344022536,344022540,344023300,344023304,344023308,344023312,344023556,344023560,344023564,344071692,344074764,344080644,344080648,344080652,344080656,344082692,344082696,344082700,344082704,352395524,352395528,352395784,352395788,352396556,352396560,352398604,352398608,352407044,352461320,352462084,352462088,352462092,352462096,352462348,352465164,352465168)    
    )
)


#From the end path, this is where we specify format for our output.
process.ALCARECOStreamPromptCalibProdPCC = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOPromptCalibProdPCC')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('ALCAPROMPT'),
        filterName = cms.untracked.string('PromptCalibProdPCC')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('rawPCC_ZB.root'),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_rawPCCProd_*_*')
)


#
process.alcaLumi = cms.Sequence(process.rawPCCProd)

#This is the key sequence that we are adding first...
process.seqALCARECOPromptCalibProdPCC = cms.Sequence(process.rawPCCProd)

process.pathALCARECOPromptCalibProdPCC = cms.Path(process.seqALCARECOPromptCalibProdPCC)

#process.seqALCARECOLumiPixels = cms.Sequence(process.siPixelDigisForLumi+process.siPixelClustersForLumi)

#process.pathALCARECOLumiPixels = cms.Path(process.seqALCARECOLumiPixels)

process.ALCARECOStreamPromptCalibProdOutPath = cms.EndPath(process.ALCARECOStreamPromptCalibProdPCC)

process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring('FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary'),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(100000)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring('warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)
#added line for additional output summary `
#process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),SkipEvent = cms.untracked.vstring('ProductNotFound') )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))


process.schedule = cms.Schedule(*[ process.pathALCARECOPromptCalibProdPCC, process.ALCARECOStreamPromptCalibProdOutPath ])
