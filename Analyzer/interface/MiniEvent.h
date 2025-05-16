#ifndef _minievent_h_
#define _minievent_h_

#include "TTree.h"

struct MiniEvent_t
{
  MiniEvent_t()
  {
    ntrk=0;
  }

  static const int MAXTRACKS     =  400;
  static const int MAXGENTRACKS  =  1000;  
  static const int MAXPROTONS    =  4;

  Bool_t isData;
  UInt_t run,lumi;
  ULong64_t event;

  Float_t weight;
  
  // Vertex info
  Int_t nvtx;

  //track info
  Int_t ntrk, trk_nMeasure[MAXTRACKS], trk_nSaturMeasure[MAXTRACKS], trk_nMeasureLayer[MAXTRACKS];
  Int_t trk_q[MAXTRACKS], trk_isPi[MAXTRACKS], trk_isK[MAXTRACKS], trk_isP[MAXTRACKS];
  Float_t trk_p[MAXTRACKS], trk_pt[MAXTRACKS], trk_eta[MAXTRACKS], trk_phi[MAXTRACKS], trk_dedx[MAXTRACKS], trk_dedxerr[MAXTRACKS];
  Float_t trk_dxy[MAXTRACKS], trk_dz[MAXTRACKS];
  
  //gen track info
  Int_t gen_ntrk;
  Float_t gen_trk_pt[MAXGENTRACKS];
  Int_t gen_trk_id[MAXGENTRACKS];
  
  // Gen info
  Int_t typevt;
  
};

void createMiniEventTree(TTree *t,MiniEvent_t &ev);
void attachToMiniEventTree(TTree *t, MiniEvent_t &ev);

#endif

