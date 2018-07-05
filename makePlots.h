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
  void Draw_HG_LG();
  void Draw_HG_LG(int BD,bool Draw_SCAT);
  void InitTH2Poly(TH2Poly& poly); //Give frame to TH2Poly
  void root_logon();

  ///////////////////////////////
  // Declaration of leaf types //
  ///////////////////////////////

   // Declaration of leaf types
   Int_t           eventID;
   vector<int>     *skirocID;
   vector<int>     *boardID;
   vector<int>     *moduleID;
   vector<int>     *channelID;
   vector<float>   *HighGainADC;
   vector<float>   *HighGainTmax;
   vector<float>   *HighGainChi2;
   vector<float>   *HighGainErrorADC;
   vector<float>   *HighGainErrorTmax;
   vector<int>     *HighGainStatus;
   vector<int>     *HighGainNCalls;
   vector<float>   *LowGainADC;
   vector<float>   *LowGainTmax;
   vector<float>   *LowGainChi2;
   vector<float>   *LowGainErrorADC;
   vector<float>   *LowGainErrorTmax;
   vector<int>     *LowGainStatus;
   vector<int>     *LowGainNCalls;
   vector<int>     *TotSlow;
   vector<int>     *ToaRise;
   vector<int>     *ToaFall;
   
};

#endif
