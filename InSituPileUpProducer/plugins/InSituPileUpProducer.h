#ifndef FastSimulation_InSituPileUpProducer_InSituPileUpProducer_H
#define FastSimulation_InSituPileUpProducer_InSituPileUpProducer_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupMixingContent.h"

#include <vector>
#include <string>
#include <fstream>

class ParameterSet;
class Event;
class EventSetup;

class TH1F;
class TFile;
class TTree;
class TBranch;
class PUEvent;

class PrimaryVertexGenerator;
class RandomEngine;

namespace Pythia8 {
  class Pythia;
}

class InSituPileUpProducer : public edm::EDProducer
{

 public:

  explicit InSituPileUpProducer(edm::ParameterSet const & p);
  virtual ~InSituPileUpProducer();
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void produce(edm::Event & e, const edm::EventSetup & c);

 private:

  PrimaryVertexGenerator* theVertexGenerator;

  double averageNumber_;
  const RandomEngine* random;
  bool usePoisson_;
  double minimumPileUpParticlePT;

  Pythia8::Pythia* pileUpPythia;

  TH1F * hprob;
  std::vector<int> dataProbFunctionVar;
  std::vector<double> dataProb;
  int varSize;
  int probSize;

};

#endif
