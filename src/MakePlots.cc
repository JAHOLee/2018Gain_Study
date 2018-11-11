#include "MakePlots.h"
#include <iostream>
#include <fstream>
#include "TCanvas.h"

MakePlots::MakePlots(){
  cout << "Constructor of MakePlots ..." << endl;
}
MakePlots::~MakePlots(){
  cout << "Destructor of MakePlots ..." << endl;
}

void MakePlots::Init_TFile(string TPro_outputname){
  Init_Pointers();
  ifstream f_check(TPro_outputname.c_str());
  file_exist = f_check.good();
  TPro_output = TPro_outputname;
  char p_name[200];
  if(file_exist){
    TPro_root = new TFile(TPro_outputname.c_str(),"update");
    TPro_history = (TTree*)TPro_root->Get("history");
    TPro_history -> SetBranchAddress("history_Run",&history_Run);

    // Get the TProfiles
    for(int BD = 0 ; BD < MAXBOARDS ; ++BD){
      for(int chip = 0 ; chip < MAXSKI ; ++chip){
	for(int ch = 0 ; ch < MAXCH ; ++ch){	  
	  sprintf(p_name,"Board%d/HGLG_chip%d_ch%d",BD,chip,ch*2);
	  HG_LG[BD][chip][ch] = (TProfile*)TPro_root->Get(p_name);
	  sprintf(p_name,"Board%d/LGTOT_chip%d_ch%d",BD,chip,ch*2);
	  LG_TOT[BD][chip][ch] = (TProfile*)TPro_root->Get(p_name);
	}
      }
    }

  }
  else{
    TPro_root = new TFile(TPro_outputname.c_str(),"recreate");
    TPro_history = new TTree("history","history");
    TPro_history-> Branch("history_Run",&history_Run);
    
    int HGLGBIN  = 400;
    int LGTOTBIN = 200;
  
    for(int BD = 0 ; BD < MAXBOARDS ; ++BD){
      for(int chip = 0 ; chip < MAXSKI ; ++chip){
	for(int ch = 0 ; ch < MAXCH ; ++ch){	  
	  sprintf(p_name,"HG_LG_BD%d_chip%d_ch%d",BD,chip,ch*2);
	  HG_LG[BD][chip][ch] = new TProfile(p_name,"",HGLGBIN,0,800,0,4000);
	  sprintf(p_name,"LG_TOT_BD%d_chip%d_ch%d",BD,chip,ch*2);
	  LG_TOT[BD][chip][ch] = new TProfile(p_name,"",LGTOTBIN,0,800,0,2000);
	  HG_LG[BD][chip][ch]->SetMarkerStyle(22);
	  HG_LG[BD][chip][ch]->SetMarkerColor(chip);
	  HG_LG[BD][chip][ch]->SetMarkerSize(1);
	  LG_TOT[BD][chip][ch]->SetMarkerStyle(21);
	  LG_TOT[BD][chip][ch]->SetMarkerColor(chip);
	  LG_TOT[BD][chip][ch]->SetMarkerSize(1);
	}
      }
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

bool MakePlots::Check_Run(int RunN){
  
  int history_run = TPro_history->GetEntries();
  bool doublefill = false;
  for(int i = 0 ; i < history_run ; ++i){
    TPro_history->GetEntry(i);
    if(history_Run == RunN){ doublefill = true; }
  }
  if(!doublefill){
    history_Run = RunN;
    TPro_history->Fill();}
  
  return doublefill;
}

void MakePlots::Write_TProfile(){
  TPro_root->cd();
  cout << "Filling TProfiles ..." << endl;
  TDirectory *dir;
  char title[50];
  
  for(int BD = 0 ; BD < MAXBOARDS ; BD++){
    sprintf(title,"Board%i",BD);
    if(!file_exist){
      dir = new TDirectory();
      dir = TPro_root->mkdir(title); }
    else{ dir = (TDirectory*)TPro_root->Get(title); }
    dir->cd();
    
    for(int chip = 0 ; chip < MAXSKI ; ++chip){
      for(int ch = 0 ; ch < MAXCH ; ++ch){
	sprintf(title,"HGLG_chip%i_ch%i",chip,ch*2);
	HG_LG[BD][chip][ch]->SetTitle(title);
	HG_LG[BD][chip][ch]->Write(title,TObject::kOverwrite);
	sprintf(title,"LGTOT_chip%i_ch%i",chip,ch*2);
	LG_TOT[BD][chip][ch]->SetTitle(title);
	LG_TOT[BD][chip][ch]->Write(title,TObject::kOverwrite);
      }
    }
  }
  
  TPro_root->cd();
  TPro_history->Write("history",TObject::kOverwrite);
  
  file_exist = true;
}
