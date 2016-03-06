#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "FastSimDataFormats/PileUpEvents/interface/PUEvent.h"

#include "FastSimulation/InSituPileUpProducer/plugins/InSituPileUpProducer.h"
#include "FastSimulation/Event/interface/BetaFuncPrimaryVertexGenerator.h"
#include "FastSimulation/Event/interface/GaussianPrimaryVertexGenerator.h"
#include "FastSimulation/Event/interface/FlatPrimaryVertexGenerator.h"
#include "FastSimulation/Event/interface/NoPrimaryVertexGenerator.h"
#include "FastSimulation/Utilities/interface/RandomEngine.h"

#include "HepMC/GenEvent.h"

#include "Pythia.h"

#include "TROOT.h"
#include "TH1F.h"

#include <iostream>
#include <memory>
#include <sys/stat.h>
#include <cmath>

InSituPileUpProducer::InSituPileUpProducer(edm::ParameterSet const & p) : 
  averageNumber_(0),
  random(0),
  usePoisson_(true),
  minimumPileUpParticlePT(0),
  pileUpPythia(0),
  hprob(0)
{    

  // This producer produces an object PileupMixingContent, needed by PileupSummaryInfo
  produces<PileupMixingContent>();
  // This producer produces a HepMCProduct, with all pileup vertices/particles
  produces<edm::HepMCProduct>("PileUpEvents");

  // Initialize the random number generator service
  edm::Service<edm::RandomNumberGenerator> rng;
  if ( ! rng.isAvailable() ) {
    throw cms::Exception("Configuration")
      << "PileUpProducer requires the RandomGeneratorService\n"
         "which is not present in the configuration file.\n"
         "You must add the service in the configuration file\n"
         "or remove the module that requires it";
  }  
  random = new RandomEngine(&(*rng));
  gRandom->SetSeed(rng->mySeed());

  // The pile-up event generation condition
  const edm::ParameterSet& pup = p.getParameter<edm::ParameterSet>("PileUpGenerator");
  usePoisson_ = pup.getParameter<bool>("usePoisson");
  if (usePoisson_) {
    std::cout << " FastSimulation/InSituPileUpProducer -> poissonian distribution" << std::endl;
    averageNumber_ = pup.getParameter<double>("averageNumber");
  } else {//take distribution from configuration
    dataProbFunctionVar = pup.getParameter<std::vector<int> >("probFunctionVariable");
    dataProb = pup.getParameter<std::vector<double> >("probValue");
    varSize = (int) dataProbFunctionVar.size();
    probSize = (int) dataProb.size();
    //    std::cout << " FastSimulation/InSituPileUpProducer -> varSize = " << varSize  << std::endl;
    //    std::cout << " FastSimulation/InSituPileUpProducer -> probSize = " << probSize  << std::endl;
    
    std::cout << " FastSimulation/InSituPileUpProducer -> distribution from configuration file "  << std::endl;
    if (probSize < varSize){
      for (int i=0; i<(varSize - probSize); i++) dataProb.push_back(0);
      edm::LogWarning("") << " The probability function data will be completed with " 
			  << (varSize - probSize) <<" values `0';"
                          << " the number of the P(x) data set after adding those 0's is " << dataProb.size();
      probSize = dataProb.size();
    }
    // Create an histogram with the data from the probability function provided by the user  
    int xmin = (int) dataProbFunctionVar[0];
    int xmax = (int) dataProbFunctionVar[varSize-1]+1;  // need upper edge to be one beyond last value
    int numBins = varSize;
    std::cout << " FastSimulation/InSituPileUpProducer -> An histogram will be created with " << numBins 
	      << " bins in the range ("<< xmin << "," << xmax << ")." << std::endl;
    hprob = new TH1F("h","Histo from the user's probability function",numBins,xmin,xmax); 
    LogDebug("") << " FastSimulation/InSituPileUpProducer -> Filling histogram with the following data:";
    for (int j=0; j < numBins ; j++){
      LogDebug("") << " x = " << dataProbFunctionVar[j ]<< " P(x) = " << dataProb[j];
      hprob->Fill(dataProbFunctionVar[j]+0.5,dataProb[j]); // assuming integer values for the bins, fill bin centers, not edges 
    }

    // Check if the histogram is normalized
    if (((hprob->Integral() - 1) > 1.0e-02) && ((hprob->Integral() - 1) < -1.0e-02)) throw cms::Exception("BadHistoDistribution") << "The histogram should be normalized!" << std::endl;
    
    // Get the averageNumber from the histo 
    averageNumber_ = hprob->GetMean();

  }

  // Initialize Pythia8 using config file
  
  edm::FileInPath fip("FastSimulation/InSituPileUpProducer/data/"+
		      pup.getParameter<std::string>("PileUpPythiaConfigFile"));
  std::string pileUpPythiaConfigFile = fip.fullPath();
  pileUpPythia = new Pythia8::Pythia();
  pileUpPythia->readFile(pileUpPythiaConfigFile.c_str());
  pileUpPythia->rndm.init(rng->mySeed());
  pileUpPythia->init();

  minimumPileUpParticlePT = pup.getParameter<double>("MinimumPileUpParticlePT");

  // Initialize the primary vertex generator
  const edm::ParameterSet& vtx = p.getParameter<edm::ParameterSet>("VertexGenerator");
  std::string vtxType = vtx.getParameter<std::string>("type");
  if ( vtxType == "Gaussian" ) 
    theVertexGenerator = new GaussianPrimaryVertexGenerator(vtx,random);
  else if ( vtxType == "Flat" ) 
    theVertexGenerator = new FlatPrimaryVertexGenerator(vtx,random);
  else if ( vtxType == "BetaFunc" )
    theVertexGenerator = new BetaFuncPrimaryVertexGenerator(vtx,random);
  else
    theVertexGenerator = new NoPrimaryVertexGenerator();

  if (averageNumber_ > 0.)
    {
      std::cout << " FastSimulation/InSituPileUpProducer ->" << std::endl
		<< " MinBias events produced using Pythia8 with configuration " 
		<< pileUpPythiaConfigFile << std::endl;
      std::cout << " with an average number of events of " << averageNumber_ << std::endl;
    }
  else std::cout << " FastSimulation/InSituPileUpProducer -> No pileup " << std::endl;


}

InSituPileUpProducer::~InSituPileUpProducer() { 
  delete theVertexGenerator;
  if (hprob) delete hprob;
}

void InSituPileUpProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
  gROOT->cd();
}

void InSituPileUpProducer::endRun(edm::Run const&, edm::EventSetup const&)
{ 
  gROOT->cd();
}
 
void InSituPileUpProducer::produce(edm::Event & iEvent, const edm::EventSetup & es)
{
  // Create the GenEvent and the HepMCProduct
  std::auto_ptr<edm::HepMCProduct> pu_product(new edm::HepMCProduct());  
  HepMC::GenEvent* evt = new HepMC::GenEvent();
  
  // How many pile-up events?
  int PUevts; float truePUevts;
  if (usePoisson_) {
    PUevts = (int) random->poissonShoot(averageNumber_);
    truePUevts = (float) averageNumber_;
  }
  else {
    float d = (float) hprob->GetRandom();
    PUevts = (int) d;
    truePUevts = d;
  }
  //  std::cout << "PUevts = " << PUevts << std::endl;

  // Save this information in the PileupMixingContent object
  // IMPORTANT: the bunch crossing number is always 0 because FastSim has no out-of-time PU
  std::auto_ptr< PileupMixingContent > PileupMixing_;

  std::vector<int> bunchCrossingList;
  bunchCrossingList.push_back(0);

  std::vector<int> numInteractionList;
  numInteractionList.push_back(PUevts);
  
  std::vector<float> trueInteractionList;
  trueInteractionList.push_back(truePUevts);
  
  PileupMixing_ = std::auto_ptr< PileupMixingContent >(new PileupMixingContent(bunchCrossingList,numInteractionList,trueInteractionList));
  iEvent.put(PileupMixing_);

  // Generate N events using Pythia8
  for (int puEvent = 0; puEvent < PUevts; ) {
    if(!pileUpPythia->next()) continue;
    
    // Smear the primary vertex and express it in mm (stupid GenEvent convention...)
    theVertexGenerator->generate();
    HepMC::FourVector smearedVertex =  
      HepMC::FourVector(theVertexGenerator->X()*10.,
			theVertexGenerator->Y()*10.,
			theVertexGenerator->Z()*10.,
			0.);
    HepMC::GenVertex* aVertex = new HepMC::GenVertex(smearedVertex);
    evt->add_vertex(aVertex);

    // Some rotation around the z axis, for more randomness
    double theAngle = random->flatShoot() * 2. * 3.14159265358979323;
    double cAngle = std::cos(theAngle);
    double sAngle = std::sin(theAngle);
    
    // Loop on particles
    for (int iTrack=0; iTrack<pileUpPythia->event.size(); ++iTrack ) {
      Pythia8::Particle& particle = pileUpPythia->event[iTrack];
      if(particle.status() > 0) {
	if(particle.isVisible()) {
	  if(particle.pT() > minimumPileUpParticlePT) {
	    // Create a FourVector, with rotation 
	    HepMC::FourVector myPart( cAngle * particle.px() + sAngle * particle.py(),
				      -sAngle * particle.px() + cAngle * particle.py(),
				      particle.pz(), particle.e());
	    // Add a GenParticle
	    HepMC::GenParticle* aGenParticle = new HepMC::GenParticle(myPart,particle.id());
	    aVertex->add_particle_out(aGenParticle);
	  }
	}
      }
    }
    // End of particle loop
    
    puEvent++;
    
  }
  // End of pile-up event loop

  // Fill the HepMCProduct from the GenEvent
  if ( evt )  { 
    pu_product->addHepMCData( evt );
    // Boost in case of beam crossing angle
    TMatrixD* boost = theVertexGenerator->boost();
    if ( boost ) pu_product->boostToLab(boost,"momentum");
  }

  // Put the HepMCProduct onto the event
  iEvent.put(pu_product,"PileUpEvents");
  // delete evt;

}

DEFINE_FWK_MODULE(InSituPileUpProducer);
