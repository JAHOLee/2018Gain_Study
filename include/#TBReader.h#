#ifndef TBReader_h
#define TBReader_h

#include "setup_config.h"
#include "MakePlots.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TApplication.h"
#include "TChain.h"
#include <vector>


class TBReader{
 public:
  // Constructor & Destructor
  TBReader();
  TBReader( TChain *c1,TChain *c2,string filename ); 
  ~TBReader();

  // Outside member
  bool Injection_Run;
  string beam_str; // "Ele","Pi","Mu"
  string fname;
  string dirpath;
  int    given_config;
  TCanvas *c1;

  // Working functions called by main.cc
  void Ntuple_Maker(setup_config *SC);
  void TProfile_Maker(setup_config *SC, MakePlots *M);
  bool Check_Config(int given_config);  // Check if config conflict with main
  
 private:
  
  //member  
  TTree        *T_Rechit;
  TTree        *T_DWC;
  int nevents;
  
  // Mainframe functions
  void Init(bool quite = false); // Run all initial stuff 
  void Init_Pointers(); // give all pointer = 0 value
  void SetRootBranch(); // Set input root branch
  void Init_Beaminfo(bool quite = false); // Read Run information from root
  bool Check_run_filled(TTree* tree); // Check if a run was already filled by
                                      // look into the Run_history tree
  
  // Ntuple members
  int  RunN;
  int  beamE; // In GeV, -1 for injection
  int  PID;   // -1 for injection, 0 for electron, 1 for pion, 2 for muon


  
  ///////////////////////////////
  // Declaration of leaf types //
  ///////////////////////////////

  /*Data*/
  //For Rechit
   UInt_t          event;
   UInt_t          run;
   Int_t           pdgID;
   Float_t         beamEnergy;
   Float_t         trueBeamEnergy;
   Int_t           NRechits;
   vector<unsigned int> *rechit_detid;
   vector<unsigned int> *rechit_module;
   vector<unsigned int> *rechit_layer;
   vector<unsigned int> *rechit_chip;
   vector<unsigned int> *rechit_channel;
   vector<unsigned int> *rechit_type;
   vector<float>   *rechit_x;
   vector<float>   *rechit_y;
   vector<float>   *rechit_z;
   vector<int>     *rechit_iu;
   vector<int>     *rechit_iv;
   vector<float>   *rechit_energy;
   vector<float>   *rechit_energy_noHG;
   vector<float>   *rechit_amplitudeHigh;
   vector<float>   *rechit_amplitudeLow;
   vector<bool>     *rechit_hg_goodFit;
   vector<bool>     *rechit_lg_goodFit;
   vector<bool>     *rechit_hg_saturated;
   vector<bool>     *rechit_lg_saturated;
   vector<bool>     *rechit_fully_calibrated;
   vector<float>   *rechit_TS2High;
   vector<float>   *rechit_TS2Low;
   vector<float>   *rechit_TS3High;
   vector<float>   *rechit_TS3Low;
   vector<unsigned short>   *rechit_Tot;
   vector<float>   *rechit_time;   
   vector<float>   *rechit_timeMaxHG;
   vector<float>   *rechit_timeMaxLG;
   vector<unsigned short>   *rechit_toaRise;
   vector<unsigned short>   *rechit_toaFall;

   // For ImpactPoints (from Delayed wire chamber)
   Int_t           ntracks;
   // ignore the layers currently
   Float_t         trackChi2_X;
   Float_t         trackChi2_Y;
   Int_t           dwcReferenceType;
   Double_t        m_x;
   Double_t        m_y;
   Double_t        b_x;
   Double_t        b_y;
  
};

#endif
