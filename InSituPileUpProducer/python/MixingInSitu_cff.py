from FastSimulation.Configuration.mixNoPU_cfi import *
from FastSimulation.InSituPileUpProducer.InSituPileUpProducer_cfi import *
# PileupSummaryInfo
from SimGeneral.PileupInformation.AddPileupSummary_cfi import *
addPileupInfo.PileupMixingLabel = 'inSituPileUpProducer'
addPileupInfo.simHitLabel = 'famosSimHits'

famosMixing = cms.Sequence(
    inSituPileUpProducer+
    addPileupInfo
)
