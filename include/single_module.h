/////////////////////////////////////////////////////////
// Arthor: Chia-hung Chien  chchien521@gmail.com       
// Just use the same class name as we used to.
// Date : 20-June-2018
/////////////////////////////////////////////////////////

#ifndef single_module_h
#define single_module_h

#include "TTree.h"
#include "TROOT.h"
#include "TH2Poly.h"
#include "TApplication.h"
#include "TChain.h"
#include <string>
#include <vector>

using namespace std;

class single_module{
 public:
  single_module( TChain *chain, string filename );
  ~single_module();
  

  void Loop();
  
  //member
  string fname;
  string dirpath;
  string inj_CH_str;
  bool   inj_sweep;
  int    inj_CH;
  int    inj_event;
  
 private:
  
  TCanvas *c1;
  TTree        *T_Rawhit;
  TFile        *root_out;
  int          nevents;

  // Mainframe functions
  void Init();
  void Root_logon();
  void Read_yaml(string yaml);
  void Setname();
  void Fill_Tprofile();
  
  //member
  string moduleID_str;
  string labelID ;
  string filepath;

  
  ///////////////////////////////
  // Declaration of leaf types //
  ///////////////////////////////

   Int_t           eventID;
   vector<int>     *skirocID;
   vector<int>     *boardID;
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
