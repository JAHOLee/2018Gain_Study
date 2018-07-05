/////////////////////////////////////////////////////////
// Arthor: Chia-hung Chien  chchien521@gmail.com       
// Just use the same class name as we used to.
// Date : 20-June-2018
/////////////////////////////////////////////////////////

#ifndef makePlots_h
#define makePlots_h

#include "TTree.h"
#include "TROOT.h"
#include "TH2Poly.h"
#include "TApplication.h"
#include "TChain.h"
#include <string>
#include <vector>

using namespace std;

class makePlots{
 public:
  makePlots( TChain *c1 , string filename);
  ~makePlots();
  
  bool DBG();
  void Loop();
  
  //member
  int  beamE;
  int  PID; // 0 for electron, 1 for pion, 2 for muon
  string beam_str; // "Ele","Pi","Mu"
  string fname;
  TCanvas *c1;
  
 private:
  //TApplication *app;
  
  TTree        *T_Rawhit;
  TH1D         *h_tprLGUS;
  TH1D         *h_LGundershoot;
  TFile        *root_out;
  int          nevents;
  
  // Mainframe functions
  void Init();
  void Init_BeamE();

  // Tool functions
  TGraph* Draw_HG_LG(int BD);
  void InitTH2Poly(TH2Poly& poly); //Give frame to TH2Poly
  void root_logon();

  ///////////////////////////////
  // Declaration of leaf types //
  ///////////////////////////////

  Int_t           eventID;
  Int_t           skirocID;
  Int_t           boardID;
  Int_t           moduleID;
  Int_t           channelID;
  Float_t         HighGainADC;
  Float_t         HighGainTmax;
  Float_t         HighGainChi2;
  Float_t         HighGainErrorADC;
  Float_t         HighGainErrorTmax;
  Int_t           HighGainStatus;
  Int_t           HighGainNCalls;
  Float_t         LowGainADC;
  Float_t         LowGainTmax;
  Float_t         LowGainChi2;
  Float_t         LowGainErrorADC;
  Float_t         LowGainErrorTmax;
  Int_t           LowGainStatus;
  Int_t           LowGainNCalls;
  Int_t           TotSlow;
  Int_t           ToaRise;
  Int_t           ToaFall;

};

#endif
