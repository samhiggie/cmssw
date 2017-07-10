import FWCore.ParameterSet.Config as cms

myana = cms.EDAnalyzer('BasicPatDistrib',
        vertices      = cms.InputTag("offlineSlimmedPrimaryVertices"),
        electrons     = cms.InputTag("slimmedElectrons"),
        beamspot      = cms.InputTag("offlineBeamSpot"),
        conversions   = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
        muons         = cms.InputTag("slimmedMuons"),
        jets          = cms.InputTag("slimmedJetsPuppi"),
        useDeepCSV    = cms.bool(False),
        mets          = cms.InputTag("slimmedMETsPuppi"),
        genParts      = cms.InputTag("packedGenParticles"),
        genJets       = cms.InputTag("slimmedGenJets"),
        #Added photons!
        photons       = cms.InputTag("slimmedPhotons"),
        # Some extra cuts you might wish to make
        #  before histograms/TTrees are filled.
        # Minimum Et
        minPhotonEt     = cms.double(3.630),
        # Minimum and max abs(eta)
        minPhotonAbsEta = cms.double(0.0),
        maxPhotonAbsEta = cms.double(3.0),
        # Minimum R9 = E(3x3) / E(SuperCluster)
        minPhotonR9     = cms.double(0.3),
        # Maximum HCAL / ECAL Energy
        maxPhotonHoverE = cms.double(0.0597),
        
        #Min IetaIeta
        minIetaIeta = cms.double(0.01031),

        # Optionally produce a TTree of photons (set to False or True).
        createPhotonTTree  = cms.bool(True),
)
