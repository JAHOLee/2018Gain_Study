/////////////////////////////////////////////////////////
// Arthor: Chia-hung Chien  chchien521@gmail.com       
// Just use the same class name as we used to.
// Date : 5-Nov-2018
/////////////////////////////////////////////////////////


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
#include <string>
#include <vector>
#include <map>
#include <fstream>

using namespace std;


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
  TCanvas *c1;

  // Working functions called by main.cc
  void Make_dir();
  void Read_Module_List();
  void Ntuple_Maker();
  void TProfile_Maker(MakePlots *M);

 private:
  
  TTree        *T_Rechit;
  TTree        *T_DWC;
  int nevents;
  int Module_List[MAXBOARDS];
  std::map<int,int> moduleID2BDorder;
  
  // Mainframe functions
  void Init(); // Run all initial stuff 
  void Init_Pointers(); // give all pointer = 0 value
  void SetRootBranch(); // Set input root branch
  void Init_Beaminfo(); // Read Run information from root
  bool DirectoryExists( const char* pzPath ); // Check if a directory exist
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
