#ifndef MakePlots_h
#define MakePlots_h

#include "setup_config.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
using namespace std;
class MakePlots{
 public:
  MakePlots();
  ~MakePlots();
  
  //TProfile members
  TProfile *HG_LG [MAXBOARDS][MAXSKI][MAXCH];
  TProfile *LG_TOT[MAXBOARDS][MAXSKI][MAXCH];
 
  // Function
  void Init_TFile(string TPro_outputname);
  void Write_TProfile();  // Takes time, fill every few runs
  bool Check_Run(int RunN);
  
  // member
  bool file_exist;
  int  history_Run;
  string TPro_output;
  
  
 private:
  // Function
  void Init_Pointers();

  TFile *TPro_root;
  TTree *TPro_history;
  
};

#endif
