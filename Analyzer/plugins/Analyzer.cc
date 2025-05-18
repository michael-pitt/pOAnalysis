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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Math/interface/deltaR.h"

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
  //edm::EDGetTokenT<TrackCollection> tracksToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> tracksToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<DeDxDataValueMap> DeDxDataToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;
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
	//tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
	tracksToken_(consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
	pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfCands"))),
	DeDxDataToken_(consumes<DeDxDataValueMap>(iConfig.getUntrackedParameter<edm::InputTag>("DeDxData"))),
	generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
	generatorlheToken_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
	genParticlesToken_(consumes<pat::PackedGenParticleCollection>(edm::InputTag("packedGenParticles"))),
	prunedGenParticlesToken_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"))),
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
  
  ev_.weight = 1;
  
  // GEN level
  edm::Handle<LHEEventProduct> evet;
  iEvent.getByToken(generatorlheToken_, evet);

  edm::Handle<GenEventInfoProduct> evt;
  iEvent.getByToken( generatorToken_,evt);

  if(evet.isValid()) {
	  ev_.typevt = evet->hepeup().IDPRUP;
	  ev_.weight = evet->hepeup().XWGTUP;
  }
  else if (evt.isValid()){
	  ev_.typevt = evt->signalProcessID();
	  ev_.weight = evt->weight();
  }
  
  // Get dE/dx collection
  Handle<DeDxDataValueMap> dEdxTrackHandle;
  iEvent.getByToken(DeDxDataToken_, dEdxTrackHandle);
  DeDxDataValueMap dEdxTrack = *dEdxTrackHandle.product();
  
  //PF candidates
  Handle<pat::PackedCandidateCollection> pfcandHandle;
  iEvent.getByToken(pfToken_,pfcandHandle);
  pat::PackedCandidateCollection pfcands = *pfcandHandle.product();

  //lostTracks 
  Handle<pat::PackedCandidateCollection> tracksHandle;
  iEvent.getByToken(tracksToken_,tracksHandle);
  pat::PackedCandidateCollection tracks = *tracksHandle.product();

  //genParticles
  edm::Handle<pat::PackedGenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_,genParticles);
  if (!genParticles.isValid()) {
    std::cerr << "Error: packedGenParticles collection is invalid" << std::endl;
  }

 
  ev_.ntrk = 0;
  Float_t dRmatch = 99.; // added
  Float_t ptmatch = 99.; // added
  Int_t matchedPdgId = 0; // added
  bool hasReco = false; // added 

  // loop over pf candidates
  for(unsigned int i=0; i<pfcands.size(); i++){
	  
	// check that the total number of reco tracks does not exceed the MAXTRACKS
	if(ev_.MAXTRACKS==ev_.ntrk){
          std::cout << "ERROR: number of reconstructed tracks reach the maximum of MAXTRACKS =  "<<ev_.MAXTRACKS<<", the ntrk loop is terminated"<<std::endl;
          std::cout <<"\t\t... consider increasing MAXTRACKS !!!"<<std::endl;
          break;
	}

	// only continue if a track is present
	auto pf_ref  = Ref<pat:: PackedCandidateCollection>( pfcandHandle, i );
        if(pf_ref->charge()==0 || !(pf_ref->hasTrackDetails())) continue;
	
        // fill basic track info
	ev_.trk_p[ev_.ntrk] = pf_ref->p() ;
	ev_.trk_pt[ev_.ntrk] = pf_ref->pt() ;
	ev_.trk_eta[ev_.ntrk] = pf_ref->eta() ;
	ev_.trk_phi[ev_.ntrk] = pf_ref->phi() ;
	ev_.trk_q[ev_.ntrk] = pf_ref->charge() ;

        ev_.trk_dxy[ev_.ntrk] = pf_ref->dxy() ;
        ev_.trk_dz[ev_.ntrk] = pf_ref->dz() ;
        ev_.trk_numberOfPixelHits[ev_.ntrk] = pf_ref->numberOfPixelHits();
        ev_.trk_numberOfHits[ev_.ntrk] = pf_ref->pseudoTrack().hitPattern().numberOfValidHits();
	
	ev_.trk_dedx[ev_.ntrk] = dEdxTrack[pf_ref].dEdx();

        //now look for gen matching information
        dRmatch = 99.;
        ptmatch = 99.;
        // loop over gen particles
        for (size_t q = 0; q < genParticles->size(); ++q) {
            const auto &p = (*genParticles)[q];
            if (p.status() != 1) continue;
              // if particle is close
              if (deltaR(p, *pf_ref) < dRmatch){
                  // if dRmatch is already below a certain value, we should decide based on pt agreement as well
                  if (dRmatch<0.1){
                    if ( deltaR(p, *pf_ref) + 0.5*std::abs(pf_ref->pt() - p.pt()) < dRmatch + 0.5*std::abs(ptmatch-pf_ref->pt()) ){
                      dRmatch=deltaR(p, *pf_ref);
                      matchedPdgId=p.pdgId();
                      ptmatch=p.pt();
                    }
                  }
                  // else just fill
                  else{
                    dRmatch=deltaR(p, *pf_ref);
                    matchedPdgId=p.pdgId();
                    ptmatch=p.pt();
                  }
                }
            }
        // fill gen matching info
        ev_.trk_matchedPdgId[ev_.ntrk] = matchedPdgId;
        ev_.trk_dRmatch[ev_.ntrk] = dRmatch;
        ev_.trk_genPt[ev_.ntrk] = ptmatch;	
	ev_.ntrk++;
	
  }


  // GEN particles: general info and tracking efficiency info
  ev_.gen_ntrk = 0;

  if(genParticles.isValid()){
      for (size_t i = 0; i < genParticles->size(); ++i){
		  
	  // check that the total number of reco tracks does not exceed the MAXGENTRACKS
	  if(ev_.MAXGENTRACKS==ev_.gen_ntrk){
            std::cout << "ERROR: number of reconstructed tracks reach the maximum of MAXGENTRACKS =  "<<ev_.MAXGENTRACKS<<", the gen_ntrk loop is terminated"<<std::endl;
            std::cout <<"\t\t... consider increasing MAXGENTRACKS !!!"<<std::endl;
            break;
	  }
	
	  auto const& p = genParticles->at(i);

	  if(p.charge()==0) continue;
          if (p.status() != 1) continue;

          // write gen info
	  ev_.gen_trk_pt[ev_.gen_ntrk] = p.pt();
          ev_.gen_trk_eta[ev_.gen_ntrk] = p.eta();
	  ev_.gen_trk_id[ev_.gen_ntrk] = p.pdgId(); 

          hasReco = false;
          // check if track matches
          for(unsigned int h=0; h<tracks.size(); h++){
                auto pf_ref2  = Ref<pat::PackedCandidateCollection>( tracksHandle, h );
                if(pf_ref2->charge()==0 || !(pf_ref2->hasTrackDetails()) ) continue;
                if (deltaR(p, *pf_ref2) < 0.1 && std::abs(pf_ref2->pt() - p.pt()) < 0.5){
                        hasReco = true;
                        continue; 
                        }      
                }
          if (hasReco==false){
          for(unsigned int h=0; h<pfcands.size(); h++){
                auto pf_ref3  = Ref<pat::PackedCandidateCollection>( pfcandHandle, h );
                if(pf_ref3->charge()==0 || !(pf_ref3->hasTrackDetails()) ) continue;

                if (deltaR(p, *pf_ref3) < 0.1 && std::abs(pf_ref3->pt() - p.pt()) < 0.5){
                        hasReco = true;
                        continue;
                        }
                }
          }
          if (hasReco){ ev_.gen_trk_hasReco[ev_.gen_ntrk] = 1;}
          else{ ev_.gen_trk_hasReco[ev_.gen_ntrk] = 0;}

	  ev_.gen_ntrk++;
	  }
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
