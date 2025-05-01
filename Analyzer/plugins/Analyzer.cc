// -*- C++ -*-
//
// Package:    pOAnalysis/Analyzer
// Class:      Analyzer
//
/**\class Analyzer Analyzer.cc pOAnalysis/Analyzer/plugins/Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Pitt
//         Created:  Sun, 27 Apr 2025 08:28:13 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "TTree.h"
#include "TH1.h"

#include "pOAnalysis/Analyzer/interface/MiniEvent.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace reco;
using reco::TrackCollection;

class Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Analyzer(const edm::ParameterSet&);
  ~Analyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<DeDxDataValueMap> DeDxDataToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif

  TTree *tree_;
  MiniEvent_t ev_;
  TH1F *h_counter;
  
  edm::Service<TFileService> fs;
  
  // apply filter
  bool applyFilt_;	
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Analyzer::Analyzer(const edm::ParameterSet& iConfig) :
	tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
    pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfCands"))),
    DeDxDataToken_(consumes<DeDxDataValueMap>(iConfig.getUntrackedParameter<edm::InputTag>("DeDxData"))),
	generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
	generatorlheToken_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
    applyFilt_( iConfig.getParameter<bool>("applyFilt") )
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
  h_counter = fs->make<TH1F>("counter", ";Counter;Events",4,0.5,4.5);
  h_counter->GetXaxis()->SetBinLabel(1,"Total");
  h_counter->GetXaxis()->SetBinLabel(2,"nPV=1");
  h_counter->GetXaxis()->SetBinLabel(3,"Veto Elastic");
  h_counter->GetXaxis()->SetBinLabel(4,"CMS-TOTEM matching");
  tree_ = fs->make<TTree>("tree","tree with selected events");
  createMiniEventTree(tree_,ev_);
}

Analyzer::~Analyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  
  h_counter->Fill(1,1);
  
  ev_.isData  = iEvent.isRealData();
  ev_.run = iEvent.id().run();
  ev_.lumi = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event();
  
  // GEN level
  edm::Handle<LHEEventProduct> evet;
  iEvent.getByToken(generatorlheToken_, evet);

  edm::Handle<GenEventInfoProduct> evt;
  iEvent.getByToken( generatorToken_,evt);

  if(evet.isValid()) {
	  ev_.typevt = evet->hepeup().IDPRUP;
  }
  else if (evt.isValid()){
	  ev_.typevt = evt->signalProcessID();
  }
  
  // Get dE/dx collection
  Handle<DeDxDataValueMap> dEdxTrackHandle;
  iEvent.getByToken(DeDxDataToken_, dEdxTrackHandle);
  DeDxDataValueMap dEdxTrack = *dEdxTrackHandle.product();
  
  //PF candidates
  Handle<pat::PackedCandidateCollection> pfcandHandle;
  iEvent.getByToken(pfToken_,pfcandHandle);
  pat::PackedCandidateCollection pfcands = *pfcandHandle.product();
  
  ev_.ntrk = 0;
  for(unsigned int i=0; i<pfcands.size(); i++){
	
	auto pf_ref  = Ref<pat:: PackedCandidateCollection>( pfcandHandle, i );
	
    if(pf_ref->charge()==0) continue;
	
	ev_.trk_p[ev_.ntrk] = pf_ref->p() ;
	ev_.trk_pt[ev_.ntrk] = pf_ref->pt() ;
	ev_.trk_eta[ev_.ntrk] = pf_ref->eta() ;
	ev_.trk_phi[ev_.ntrk] = pf_ref->phi() ;
	ev_.trk_q[ev_.ntrk] = pf_ref->charge() ;
	
	ev_.trk_dedx[ev_.ntrk] = dEdxTrack[pf_ref].dEdx();
	
	// det dE/dX from the track:
	
	ev_.ntrk++;
	
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif

  tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void Analyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void Analyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
