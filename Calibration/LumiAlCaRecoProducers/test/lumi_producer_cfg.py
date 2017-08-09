#########################
#Author: Sam Higginbotham
#Purpose: To investigate the AlCaPCCProducer input and output. 
#########################
import FWCore.ParameterSet.Config as cms

process = cms.Process("LUMI")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/shigginb/cmssw/CMSSW_9_2_0/src/rawPCC_297411_ZB.root'),
    processingMode = cms.untracked.string('RunsAndLumis')
)
#Added process to select the appropriate events 
process.OutALCARECOPromptCalibProdPCC = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOPromptCalibProdPCC')
    ),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_lumiPCCProd_*_*')
)

#Trial to see if an SQLite file can be READ locally - will use this for the LumiPCCProducer..
process.load("CondCore.CondDB.CondDB_cfi")
#process.GlobalTag.globaltag = '92X_dataRun2_Prompt_v4'

process.CondDB.connect = 'sqlite_file:/afs/cern.ch/work/s/shigginb/cmssw/CMSSW_9_2_0/src/testcorrLumi.db'

process.PoolDBESSource = cms.ESSource("PoolDBESSource",
    process.CondDB,
    DumpStat=cms.untracked.bool(True),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('MyLumiCorrectionsRcd'),
        tag = cms.string("TestCorrections")
    )),
)


#Make sure that variables match in producer.cc and .h
process.lumiPCCProd = cms.EDProducer("LumiPCCProducer",
    LumiPCCProducerParameters = cms.PSet(
        #Mod factor to count lumi and the string to specify output 
        PCCobLabel = cms.string("rawPCCProd"),
        ProdInst = cms.string("rawPCZeroBias"),
        resetEveryNLumi = cms.untracked.int32(1),
        trigstring = cms.untracked.string("rawPCCtest"), 
)# ,
 #   toGet = cms.VPSet(cms.PSet( 
 #       record = cms.string('MyLumiCorrectionsRcd'),
 #       tag = cms.string("TestCorrections")
 #   ))
)

#process.get = cms.EDAnalyzer("EventSetupRecordDataGetter",
#    toGet = cms.VPSet(cms.PSet(
#        record = cms.string('Corrections'),
#        data = cms.vstring('TestLSBasedCorrLumi')
#    )),
#    verbose = cms.untracked.bool(True)
#)
#
# A data source must always be defined. We don't need it, so here's a dummy one.
#process.source = cms.Source("EmptyIOVSource",
#    timetype = cms.string('lumiid'),
#    firstValue = cms.uint64(1),
#    lastValue = cms.uint64(1),
#    interval = cms.uint64(1)
#)
#
#process.path = cms.Path(process.get)



#process.GlobalTag.toGet.append(
# cms.PSet(
#   connect = cms.string('sqlite_file:/afs/cern.ch/work/s/shigginb/cmssw/CMSSW_9_1_1_patch1/src/testcorrLumi.db'),
#   record = cms.string("Corrections"),
#   tag = cms.string("TestLSBasedCorrLumi")),
#)





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
    fileName = cms.untracked.string('correctedLumi.root'),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_lumiPCCProd_*_*')
)


#
process.alcaLumi = cms.Sequence(process.lumiPCCProd)

#This is the key sequence that we are adding first...
process.seqALCARECOPromptCalibProdPCC = cms.Sequence(process.lumiPCCProd)

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
            reportEvery = cms.untracked.int32(1000000)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1000000)
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
