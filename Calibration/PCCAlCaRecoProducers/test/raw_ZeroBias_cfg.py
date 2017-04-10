#########################
#Author: Sam Higginbotham
#Purpose: To investigate the AlCaPCCProducer input and output. 
#########################
import FWCore.ParameterSet.Config as cms

process = cms.Process("rawRECO")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/shigginb/cmssw/CMSSW_9_1_0_pre1/src/ProdPCC_Random_100.root')
)

#Added process to select the appropriate events 
process.OutALCARECOPromptCalibProdPCC = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOPromptCalibProdPCC')
    ),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_rawPCCProd_*_*')
)
#Troubleshooting the Product Not Found issue:
process.content = cms.EDAnalyzer("EventContentAnalyzer")

#Make sure that variables match in producer.cc and .h
process.rawPCCProd = cms.EDProducer("rawPCCProducer",
#process.rawPCCProd = cms.EDAnalyzer("rawPCCProducer",
    rawPCCProducerParameters = cms.PSet(
        WriteToDB = cms.bool(False),
        #Mod factor to count lumi and the string to specify output 
        PCCobLabel = cms.string("alcaPCCProducer"),
        ProdInst = cms.string(""),
        resetEveryNLumi = cms.untracked.int32(1),
        trigstring = cms.untracked.string("rawPCCtest") 
    ),
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
        'keep *_rawPCCProd_*_*', 
        'keep *_MEtoEDMConvertSiStrip_*_*')
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
