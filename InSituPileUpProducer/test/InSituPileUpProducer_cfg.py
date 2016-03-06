import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")

# Number of events to be generated
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Include DQMStore, needed by the famosSimHits
process.DQMStore = cms.Service( "DQMStore")

# Include the RandomNumberGeneratorService definition
process.load("IOMC.RandomEngine.IOMC_cff")

# Generate ttbar events
process.source = cms.Source("EmptySource")
process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(14000.0),
    PythiaParameters = cms.PSet(
    processParameters = cms.vstring(
    'Main:timesAllowErrors    = 10000',
    'WeakSingleBoson:ffbar2gmZ = on',
    'PhaseSpace:mHatMin = 60.',
    'PhaseSpace:mHatMax = -1',
    '23:onMode = 0',
    '23:onIfMatch = 13 -13'
    'Tune:pp 2',
    'Tune:ee 3'),
    parameterSets = cms.vstring('processParameters')
    )
)

# Famos sequences (NO HLT)
#process.load("FastSimulation.Configuration.CommonInputs_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('FastSimulation.Configuration.Geometries_cff')
process.load("FastSimulation.Configuration.FamosSequences_cff")

# Parametrized magnetic field (new mapping, 4.0 and 3.8T)
#process.load("Configuration.StandardSequences.MagneticField_40T_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = True

# If you want to turn on/off pile-up
process.load("FastSimulation.InSituPileUpProducer.MixingInSitu_cff")
process.RandomNumberGeneratorService.inSituPileUpProducer = cms.PSet(
        initialSeed = cms.untracked.uint32(12345),
        engineName = cms.untracked.string('TRandom3')
    )
process.inSituPileUpProducer.PileUpGenerator.averageNumber = 140

#process.load('FastSimulation.PileUpProducer.mix_2012_Startup_inTimeOnly_cff')
# You may not want to simulate everything for your study
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True

# Get frontier conditions    - not applied in the HCAL, see below
from HLTrigger.Configuration.AutoCondGlobalTag import AutoCondGlobalTag
process.GlobalTag = AutoCondGlobalTag(process.GlobalTag,'auto:startup_GRun')
# Allow reading of the tracker geometry from the DB
process.load('CalibTracker/Configuration/Tracker_DependentRecords_forGlobalTag_nofakes_cff')

# Apply ECAL miscalibration
process.ecalRecHit.doMiscalib = True

# Apply Tracker misalignment
process.famosSimHits.ApplyAlignment = True
process.misalignedTrackerGeometry.applyAlignment = True
process.misalignedDTGeometry.applyAlignment = True
process.misalignedCSCGeometry.applyAlignment = True

#  Attention ! for the HCAL IDEAL==STARTUP
#process.caloRecHits.RecHitsFactory.HCAL.Refactor = 1.0
#process.caloRecHits.RecHitsFactory.HCAL.Refactor_mean = 1.0
#process.caloRecHits.RecHitsFactory.HCAL.fileNameHcal = "hcalmiscalib_0.0.xml"

# Famos with everything !
#process.p1 = cms.Path(process.ProductionFilterSequence*process.famosWithEverything)
#process.source = cms.Source("EmptySource")

process.p1 = cms.Path(process.generator*process.famosWithEverything)

# To write out events
process.load("FastSimulation.Configuration.EventContent_cff")
process.o1 = cms.OutputModule("PoolOutputModule",
    process.AODSIMEventContent,
    fileName = cms.untracked.string('InSituPileUpProducerTest-PU140.root')
)
process.outpath = cms.EndPath(process.o1)

# Keep output to a nice level
# process.Timing =  cms.Service("Timing")
# process.MessageLogger.destinations = cms.untracked.vstring("pyDetailedInfo.txt","cout")
# process.MessageLogger.categories.append("FamosManager")
# process.MessageLogger.cout = cms.untracked.PSet(threshold=cms.untracked.string("INFO"),
#                                                 default=cms.untracked.PSet(limit=cms.untracked.int32(0)),
#                                                 FamosManager=cms.untracked.PSet(limit=cms.untracked.int32(100000)))


# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

# Get process summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Get modules used
dump_file = open("dump_file.py", "w")
dump_file.write(process.dumpPython())
