Repository for collecting recipes for standard physics objects for analysis of simulated events with the CMS phase 2 detector.
=========================


Contains (in the future)

a) scripts/recipes to set up the correct CMSSW environment

b) a collection of latest recipes and a CMSSW config to apply them to get collections of standard objects suitable for physics analysis as CMSSW producers (e.g. recommendedTightMuons = tightMuonProducer(slimmedMuons) etc.)

c) an ntuple production configuration file that produces an ntuple that can be analysed in a similar manner als Delphes samples (with the DAnalysis framework)


Installation
--------------

```bash
cmsrel CMSSW_9_1_1_patch1
cd CMSSW_9_1_1_patch1/src
cmsenv
git cms-addpkg RecoEgamma/EgammaIsolationAlgos
cd RecoEgamma
git clone git@github.com:nsmith-/Phase2InterimID.git
cd ..
git clone git@github.com:jkiesele/PhaseTwoAnalysis.git
cp PhaseTwoAnalysis/RecoEgammaFix/* RecoEgamma/EgammaIsolationAlgos/plugins/
scram b -j8
```

How to run PAT on RECO datasets
----------------

The `PatProducer` folder contains a configuration file to produce miniAOD files from RECO files. Interactively, after updating the list of input files, one can run
```bash
cmsRun miniAOD-prod_PAT.py
```
A skeleton of crab configuration file is also provided in this folder. The following fields need to be updated:
   * `config.General.requestName` 
   * `config.Data.inputDataset`
   * `config.Data.outLFNDirBase`

before running
```bash
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crabConfig.py
```

Plotting basic distributions from RECO collections
-----------------

A basic EDAnalyzer is available in the `BasicRecoDistrib` folder. Several private functions handle electron and forward muon ID. Lepton isolation is computed with a simple loop over neighbouring particles and there is no b-tagging information. Normalization to luminosity is not handled. More details are given in the `implementation` section of the `.cc` file.
After updating the list of input files, the analyzer can be run interactively from the `test` subfolder :
```bash
cmsRun ConfFile_cfg.py
```
Befor the EDAnalyzer, PUPPI is run on the fly and jets are re-clustered. The MET is also recomputed but not exactly with the official recipe (that needs PAT collections).

Plots in a pdf format can be obtained by running:
```bash
root -l 
.L plotIt.C++
plotIt()
```

A skeleton of crab configuration file is also provided. The following fields need to be updated:
   * `config.General.requestName` 
   * `config.Data.inputDataset`
   * `config.Data.outLFNDirBase`

before running
```bash
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crabConfig.py
```


Plotting basic distributions from PAT collections
-----------------

A basic EDAnalyzer is available in the `BasicPatDistrib` folder. Several private functions handle central electron and forward muon ID. A flag (`useDeepCSV`) can be set to true in the configuration file to use deepCSV rather than CSVv2 as b-tagging discriminant. Normalization to luminosity is not handled. More details are given in the `implementation` section of the `.cc` file.
After updating the list of input files, the analyzer can be run interactively from the `test` subfolder :
```bash
cmsRun ConfFile_cfg.py
```

Plots in a pdf format can be obtained by running:
```bash
root -l 
.L plotIt.C++
plotIt()
```

A skeleton of crab configuration file is also provided. The following fields need to be updated:
   * `config.General.requestName` 
   * `config.Data.inputDataset`
   * `config.Data.outLFNDirBase`

before running
```bash
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crabConfig.py
```

Producing flat ntuples
-----------------

Flat ntuples can be produced in the `NTupler` folder, either from PAT or RECO events, by running interactively:
```bash
cmsRun scripts/produceNtuples_cfg.py skim=False/True outFilename=MiniEvents.root inputFormat=RECO/PAT
```

The `skim` flag can be used to reduce the size of the output files. A histogram containing the number of events before the skim is then stored in the output files. By default, events are required to contain at least 1 lepton and 2 jets, but this can be easily modified ll.71-97 of `src/produceNtuples_cfg.py`.

The structure of the output tree can be seen/modified in `interface/MiniEvent.h` and `src/MiniEvent.cc`.

The main analyzers are:
   * `plugins/MiniFromPat.cc` -- to run over PAT events 
   * `plugins/MiniFromReco.cc` -- to run over RECO events 

Details on the object definitions are given in the `implementation` section.

A skeleton of crab configuration file is also provided. The following fields need to be updated:
   * `config.General.requestName` 
   * `config.Data.inputDataset`
   * `config.Data.outLFNDirBase`

before running
```bash
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crabConfig.py
```

To adjust the input parameters of `scripts/produceNtuples_cfg.py`, the three following fields need to be modified consistently:
   * `config.JobType.pyCfgParams`
   * `config.JobType.inputFiles`
   * `config.JobType.outputFiles`
