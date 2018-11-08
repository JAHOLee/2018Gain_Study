#ifndef MakePlots_h
#define MakePlots_h

#include "setup_config.h"
#include "TProfile.h"


class MakePlots{
 public:
  MakePlots();
  ~MakePlots();
  
  //TProfile members
  TProfile *HG_LG [MAXBOARDS][MAXSKI][MAXCH];
  TProfile *LG_TOT[MAXBOARDS][MAXSKI][MAXCH];

  //Function
  void Init_TProfile(int BD); //Init a module (4*32plots)
  void Init_Pointers();
};

#endif
