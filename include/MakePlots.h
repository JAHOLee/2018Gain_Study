#ifndef MakePlots_h
#define MakePlots_h

#include "setup_config.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
class MakePlots{
 public:
  MakePlots(setup_config *SC);
  ~MakePlots();
  
  //TProfile members
  TProfile *HG_LG [MAXBOARDS][MAXSKI][MAXCH];
  TProfile *LG_TOT[MAXBOARDS][MAXSKI][MAXCH];
 
  // Function
  bool Init_TFile(string TPro_outputname);
  void Write_TProfile();  // Takes time, fill every few runs
  bool Check_Run(int RunN);

  // Drawing Function
  void root_logon();
  
  // member
  bool file_exist;
  int  history_Run;
  string TPro_output;
  
  
 private:
  // Function
  void Init_Pointers();

  TFile *TPro_root;
  TTree *TPro_history;
  setup_config *mysetup;
};

#endif
