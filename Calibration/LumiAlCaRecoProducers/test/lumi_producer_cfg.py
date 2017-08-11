#########################
#Author: Sam Higginbotham
#Purpose: To investigate the AlCaPCCProducer input and output. 
#########################
import FWCore.ParameterSet.Config as cms

process = cms.Process("LUMI")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/shigginb/cmssw/CMSSW_9_2_0/src/rawPCC_297411_ZB.root'),
    processingMode=cms.untracked.string('RunsLumisAndEvents'),

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
    connect = cms.string('sqlite_file:/afs/cern.ch/work/s/shigginb/cmssw/CMSSW_9_2_0/src/testcorrLumi.db')
)


#Make sure that variables match in producer.cc and .h
process.lumiPCCProd = cms.EDProducer("LumiPCCProducer",
    LumiPCCProducerParameters = cms.PSet(
        #Mod factor to count lumi and the string to specify output 
        PCCobLabel = cms.string("rawPCCProd"),
        ProdInst = cms.string("rawPCZeroBias"),
        resetEveryNLumi = cms.untracked.int32(1),
        trigstring = cms.untracked.string("rawPCCtest"), 
        label = cms.untracked.string("/afs/cern.ch/work/s/shigginb/cmssw/CMSSW_9_2_0/src/297411test.csv"), 
)# ,
 #   toGet = cms.VPSet(cms.PSet( 
 #       record = cms.string('LumiCorrectionsRcd'),
 #       tag = cms.string("TestCorrections")
 #   ))
)



#process.source = cms.Source("EmptySource",
#                            processingMode=cms.untracked.string('RunsLumisAndEvents'),
#                            firstRun = cms.untracked.uint32(297411),
#                            lastRun = cms.untracked.uint32(297411),
#                            numberEventsInLuminosityBlock = cms.untracked.uint32(1),
#                            numberEventsInRun=cms.untracked.uint32(10),
#                            firstLuminosityBlock = cms.untracked.uint32(1)                           
#                            )
#



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
