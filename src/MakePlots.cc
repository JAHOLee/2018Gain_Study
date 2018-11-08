#include "MakePlots.h"
#include <iostream>

MakePlots::MakePlots(){}
MakePlots::~MakePlots(){}

void MakePlots::Init_TProfile(int BD){
  int HGLGBIN  = 400;
  int LGTOTBIN = 200;
  
  char p_name[200];
  
  for(int chip = 0 ; chip < MAXSKI ; ++chip){
    for(int ch = 0 ; ch < MAXCH ; ++ch){
      
      sprintf(p_name,"HG_LG_BD%d_chip%d_ch%d",BD,chip,ch*2);
      HG_LG[BD][chip][ch] = new TProfile(p_name,"",HGLGBIN,0,800,0,4000);
      sprintf(p_name,"LG_TOT_BD%d_chip%d_ch%d",BD,chip,ch*2);
      LG_TOT[BD][chip][ch] = new TProfile(p_name,"",LGTOTBIN,0,800,0,2000);
    }
  }

}
void MakePlots::Init_Pointers(){
  for(int BD = 0 ; BD < MAXBOARDS ; ++BD){
    for(int chip = 0 ; chip < MAXSKI ; ++chip){
      for(int ch = 0 ; ch < MAXCH ; ++ch){
	HG_LG [BD][chip][ch] = 0;
        LG_TOT[BD][chip][ch] = 0;
      }
    }
  }
}
