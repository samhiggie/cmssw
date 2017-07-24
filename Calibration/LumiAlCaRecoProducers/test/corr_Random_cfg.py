#########################
#Author: Sam Higginbotham
#Purpose: To investigate the AlCaPCCProducer input and output. 
#########################
import FWCore.ParameterSet.Config as cms

process = cms.Process("corrRECO")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/shigginb/cmssw/CMSSW_9_2_0/src/rawPCC_297411_RD.root')
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/shigginb/cmssw/CMSSW_9_2_0/src/rawPCC_297227_RD.root')
)
#Added process to select the appropriate events 
process.OutALCARECOPromptCalibProdPCC = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOPromptCalibProdPCC')
    ),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_corrPCCProd_*_*')
)

#Make sure that variables match in producer.cc and .h
process.corrPCCProd = cms.EDProducer("CorrPCCProducer",
    CorrPCCProducerParameters = cms.PSet(
        #Mod factor to count lumi and the string to specify output 
        inLumiObLabel = cms.string("rawPCCProd"),
        ProdInst = cms.string("rawPCRandom"),
        resetEveryNLumi=cms.int32(50),
        trigstring = cms.untracked.string("corrPCCRand"), 
        type2_a= cms.double(0.00086),
        type2_b= cms.double(0.014),
    )
)

#output file service for the histograms 
#process.TFileService = cms.Service("TFileService",
#                                       fileName = cms.string('histodemo2.root')
#                                   )

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
    fileName = cms.untracked.string('corrPCC_297411_RD.root'),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_corrPCCProd_*_*')
)

#Output for the Database
process.load("CondCore.CondDB.CondDB_cfi")
#process.load("CondCore.DBCommon.CondDBCommon_cfi")


process.CondDB.connect = "sqlite_file:testcorrLumi.db"
#process.source = cms.Source("EmptySource")

#process.source = cms.Source("EmptyIOVSource",
#    timetype = cms.string('lumiid'),
#    firstValue = cms.uint64(1),
#    lastValue = cms.uint64(1),
#    interval = cms.uint64(1)
#)

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDB,
    toPut = cms.VPSet(
    cms.PSet(record = cms.string('MyLumiCorrectionsRcd'),tag = cms.string('TestCorrections'))#,
     ),
    loadBlobStreamer = cms.untracked.bool(False),
    timetype   = cms.untracked.string('lumiid')
)

#
process.alcaLumi = cms.Sequence(process.corrPCCProd)

#This is the key sequence that we are adding first...
process.seqALCARECOPromptCalibProdPCC = cms.Sequence(process.corrPCCProd)

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
