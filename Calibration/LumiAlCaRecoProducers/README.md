#Project Title 
These are the instructions to setup the environment for and to produce Luminosity using the Pixel Cluster Counting method. Please follow the installation instructions before proceeding to the run steps. 


### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be


```

cmsrel CMSSW_9_2_0 
cd CMSSW_9_2_0/src
cmsenv
git cms-init 
git cms-addpkg Calibration/LumiAlCaRecoProducers 
git cms-addpkg DataFormats/Luminosity 
git cms-addpkg CondFormats/LumiCorrections
git cms-addpkg CondFormats/DataRecord
git cms-addpkg CondCore/LumiCorrectionsPlugin/ 

scram b -j 16

```



## Running


First complete the corrections using the Random triggered data set. 
To change the configuration you may edit the file:
```
Calibration/LumiAlCaRecoProducers/test/raw_corr_Random_cfg.py 
```
You may also test that it is working locally by doing:
```
cmsRun Calibration/LumiAlCaRecoProducers/test/raw_corr_Random_cfg.py 
```
To run a crab job using this file you may use:
```
crab submit -c  Calibration/LumiAlCaRecoProducers/test/crab3_raw_corr.py
```
If you decide to edit the naming conventions for the input and output please make sure that they match the corresponding string in the 
 JobType.outputFiles line within the crab job. 
orocess.outpath = cms.EndPath(process.out)

Second complete the application of the corrections using the ZeroBias triggered data set. 
To change the configuration you may edit the file:
```
Calibration/LumiAlCaRecoProducers/test/raw_lumi_ZeroBias_cfg.py
```
You may also test that it is working locally by doing:
```
cmsRun Calibration/LumiAlCaRecoProducers/test/raw_lumi_ZeroBias_cfg.py 
```
To run a crab job using this file you may use:
```
crab submit -c  Calibration/LumiAlCaRecoProducers/test/crab3_raw_lumi.py
```

If you decide to edit the naming conventions for the input and output please make sure that they match the corresponding string in the 
 JobType.outputFiles line and the JobType.inputFiles line within the crab job.


## Built With

Reminder: CMSSW_9_2_0

