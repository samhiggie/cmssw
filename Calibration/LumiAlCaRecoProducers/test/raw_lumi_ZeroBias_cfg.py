#########################
#Author: Sam Higginbotham
#########################
import FWCore.ParameterSet.Config as cms

process = cms.Process("rawRECO")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'file:/eos/cms/store/data/Run2017B/AlCaLumiPixels/ALCARECO/AlCaPCCRandom-PromptReco-v1/000/297/411/00000/0AC6C13D-F259-E711-AB3A-02163E019E8E.root',
    #'file:/eos/cms/store/data/Run2017B/AlCaLumiPixels/ALCARECO/AlCaPCCRandom-PromptReco-v1/000/297/411/00000/1E2F6E5E-0C5A-E711-88B4-02163E01A420.root',
    'file:/eos/cms/store/data/Run2017B/AlCaLumiPixels/ALCARECO/AlCaPCCZeroBias-PromptReco-v1/000/297/411/00000/2062F49E-0B5A-E711-8D71-02163E011F08.root',
    'file:/eos/cms/store/data/Run2017B/AlCaLumiPixels/ALCARECO/AlCaPCCZeroBias-PromptReco-v1/000/297/411/00000/2885844A-F559-E711-885F-02163E0146CA.root',
    'file:/eos/cms/store/data/Run2017B/AlCaLumiPixels/ALCARECO/AlCaPCCZeroBias-PromptReco-v1/000/297/411/00000/6A2A8967-F159-E711-A818-02163E01A3B8.root')
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
        trigstring = cms.untracked.string("rawPCCZeroBias"), 
        #Below is a list of module IDs that will be ignored in calculation of luminosity
        modVeto=cms.vint32(303042564,303042568,303042572,303042576,303042580,303042584,303042588,303042592,303046660,303046664,303046668,303046672,303046676,303046680,303046684,303046688,303050756,303050760,303050764,303050768,303050772,303050776,303050780,303050784,303054852,303054856,303054860,303054864,303054868,303054872,303054876,303054880,303058948,303058952,303058956,303058960,303058964,303058968,303058972,303058976,303063044,303063048,303063052,303063056,303063060,303063064,303063068,303063072,303067140,303067144,303067148,303067152,303067156,303067160,303067164,303067168,303071236,303071240,303071244,303071248,303071252,303071256,303071260,303071264,303075332,303075336,303075340,303075344,303075348,303075352,303075356,303075360,303079428,303079432,303079436,303079440,303079444,303079448,303079452,303079456,303083524,303083528,303083532,303083536,303083540,303083544,303083548,303083552,303087620,303087624,303087628,303087632,303087636,303087640,303087644,303087648)
    )
)

process.load("CondCore.CondDB.CondDB_cfi")
#process.GlobalTag.globaltag = '92X_dataRun2_Prompt_v4'
#process.load("CondCore.DBCommon.CondDBSetup_cfi")

#process.CondDB.connect = 'sqlite_file:/afs/cern.ch/work/s/shigginb/cmssw/CMSSW_9_2_0/src/testcorrLumi.db'

process.PoolDBESSource = cms.ESSource("PoolDBESSource",
    #process.CondDB,
    #process.CondDBSetup,
    DumpStat=cms.untracked.bool(True),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('LumiCorrectionsRcd'),
        tag = cms.string("TestCorrections")
    )),
    connect = cms.string('sqlite_file:GlobalCorr.db')
)


#Make sure that variables match in producer.cc and .h
process.lumiPCCProd = cms.EDProducer("LumiPCCProducer",
    LumiPCCProducerParameters = cms.PSet(
        #Mod factor to count lumi and the string to specify output 
        PCCobLabel = cms.string("rawPCCProd"),
        ProdInst = cms.string("rawPCCZeroBias"),
        resetEveryNLumi = cms.untracked.int32(1),
        trigstring = cms.untracked.string("rawPCCtest"), 
        label = cms.string("lumi.csv"), 
)
)

#From the end path, this is where we specify format for our output.
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('raw_lumi_PCC_ZB.root'),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_rawPCCProd_*_*',
        'keep *_corrPCCProd_*_*')
)

process.path1 = cms.Path(process.rawPCCProd+process.lumiPCCProd)


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
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.outpath = cms.EndPath(process.out)

