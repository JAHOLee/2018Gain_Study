#include "TBReader.h"
#include "MakePlots.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <utility>
#include <sstream>
#include <dirent.h>
TProfile *HGLG [MAXBOARDS][MAXSKI][MAXCH];
TProfile *LGTOT[MAXBOARDS][MAXSKI][MAXCH];
int RUNNUM;
TFile *f;
TTree *history;
bool dirty_way_fexist; 

TBReader::TBReader(){
  cout << "Constructor of TBReader, blank constructor... " << endl;
}

TBReader::TBReader( TChain *c1, TChain *c2, string filename):T_Rechit(c1),T_DWC(c2){
  cout << "Constructor of TBReader, Test Beam Run ... \n\n" << endl;
  fname = filename;
}

TBReader::~TBReader(){
  cout << "\n\n";
  cout << "Destructor of TBReader ... " << endl;
}

void TBReader::Init_Pointers(){
  rechit_detid = 0;
  rechit_module = 0;
  rechit_layer = 0;
  rechit_chip = 0;
  rechit_channel = 0;
  rechit_type = 0;
  rechit_x = 0;
  rechit_y = 0;
  rechit_z = 0;
  rechit_iu = 0;
  rechit_iv = 0;
  rechit_energy = 0;
  rechit_energy_noHG = 0;
  rechit_amplitudeHigh = 0;
  rechit_amplitudeLow = 0;
  rechit_hg_goodFit = 0;
  rechit_lg_goodFit = 0;
  rechit_hg_saturated = 0;
  rechit_lg_saturated = 0;
  rechit_fully_calibrated = 0;
  rechit_TS2High = 0;
  rechit_TS2Low = 0;
  rechit_TS3High = 0;
  rechit_TS3Low = 0;
  rechit_Tot = 0;
  rechit_time = 0;
  rechit_timeMaxHG = 0;
  rechit_timeMaxLG = 0;
  rechit_toaRise = 0;
  rechit_toaFall = 0;
   
}

void TBReader::SetRootBranch(){
  nevents = T_Rechit->GetEntries();

  T_Rechit->SetBranchAddress("event", &event);
  T_Rechit->SetBranchAddress("run", &run);
  T_Rechit->SetBranchAddress("pdgID", &pdgID);
  T_Rechit->SetBranchAddress("beamEnergy", &beamEnergy);
  T_Rechit->SetBranchAddress("trueBeamEnergy", &trueBeamEnergy);
  T_Rechit->SetBranchAddress("NRechits", &NRechits);
  T_Rechit->SetBranchAddress("rechit_detid", &rechit_detid);
  T_Rechit->SetBranchAddress("rechit_module", &rechit_module);
  T_Rechit->SetBranchAddress("rechit_layer", &rechit_layer);
  T_Rechit->SetBranchAddress("rechit_chip", &rechit_chip);
  T_Rechit->SetBranchAddress("rechit_channel", &rechit_channel);
  T_Rechit->SetBranchAddress("rechit_type", &rechit_type);

  T_Rechit->SetBranchAddress("rechit_x", &rechit_x);
  T_Rechit->SetBranchAddress("rechit_y", &rechit_y);
  T_Rechit->SetBranchAddress("rechit_z", &rechit_z);
  T_Rechit->SetBranchAddress("rechit_iu", &rechit_iu);
  T_Rechit->SetBranchAddress("rechit_iv", &rechit_iv);
  T_Rechit->SetBranchAddress("rechit_energy", &rechit_energy);
  T_Rechit->SetBranchAddress("rechit_energy_noHG", &rechit_energy_noHG);

  T_Rechit->SetBranchAddress("rechit_amplitudeHigh", &rechit_amplitudeHigh);
  T_Rechit->SetBranchAddress("rechit_amplitudeLow", &rechit_amplitudeLow);
  T_Rechit->SetBranchAddress("rechit_hg_goodFit", &rechit_hg_goodFit);
  T_Rechit->SetBranchAddress("rechit_lg_goodFit", &rechit_lg_goodFit);
  T_Rechit->SetBranchAddress("rechit_hg_saturated", &rechit_hg_saturated);
  T_Rechit->SetBranchAddress("rechit_lg_saturated", &rechit_lg_saturated);
  T_Rechit->SetBranchAddress("rechit_fully_calibrated", &rechit_fully_calibrated);
  T_Rechit->SetBranchAddress("rechit_TS2High", &rechit_TS2High);
  T_Rechit->SetBranchAddress("rechit_TS2Low", &rechit_TS2Low);
  T_Rechit->SetBranchAddress("rechit_TS3High", &rechit_TS3High);
  T_Rechit->SetBranchAddress("rechit_TS3Low", &rechit_TS3Low);
    
  T_Rechit->SetBranchAddress("rechit_Tot", &rechit_Tot);
  T_Rechit->SetBranchAddress("rechit_time", &rechit_time);
  T_Rechit->SetBranchAddress("rechit_timeMaxHG", &rechit_timeMaxHG);
  T_Rechit->SetBranchAddress("rechit_timeMaxLG", &rechit_timeMaxLG);
  T_Rechit->SetBranchAddress("rechit_toaRise", &rechit_toaRise);
  T_Rechit->SetBranchAddress("rechit_toaFall", &rechit_toaFall);

  T_DWC->SetBranchAddress("ntracks", &ntracks);
  T_DWC->SetBranchAddress("trackChi2_X", &trackChi2_X);
  T_DWC->SetBranchAddress("trackChi2_Y", &trackChi2_Y);
  T_DWC->SetBranchAddress("dwcReferenceType", &dwcReferenceType);
  T_DWC->SetBranchAddress("m_x", &m_x);
  T_DWC->SetBranchAddress("m_y", &m_y);
  T_DWC->SetBranchAddress("b_x", &b_x);
  T_DWC->SetBranchAddress("b_y", &b_y);

}

void TBReader::Init_Beaminfo(){
  T_Rechit->GetEntry(0);
  beamE = beamEnergy;
  if( pdgID == 11 ){
    beam_str = "Electron";
    PID = 0;}
  else if( pdgID == 13){
    beam_str = "Muon";
    PID = 2;}
  else if( pdgID == 211){
    beam_str = "Pion";
    PID = 1;}
  else{
    cout << "unknown PDGID QQ" << endl;
    beam_str = "??";
    PID = -1;}
  RunN = run;
  cout << beam_str.c_str() << " Run, "<< beamE << "GeV , with " << nevents
       << " events." << endl;
}

void TBReader::Init(){
  Init_Pointers();
  SetRootBranch();
  Init_Beaminfo();
}

void TBReader::Read_Module_List(){
  // Using the first configuration so all boards are included
  string Module_Layout = "./configs/all_config.csv";
  // 6-Nov-2018 copy from
  // https://docs.google.com/spreadsheets/d/1KcvFr3JG69plQeVy4AA8nI63ZuTMc6hGDWGwo-F5xPE/edit#gid=337823033
  
  ifstream infile(Module_Layout.c_str());
  string line;
  int line_count = 0;
  int members = 6;
  string line_contents[members];
  
  // Get the headers
  getline(infile,line); 
  getline(infile,line);
  
  while(true){
    getline(infile,line);
    if( infile.eof() ) {break;};
    std::istringstream iss(line);
    for(int i = 0 ; i < members ; ++i){
      getline(iss,line_contents[i], ',' );}
    int ModuleID = std::stoi( line_contents[0] );
    Module_List[line_count] = ModuleID;
    moduleID2BDorder.insert( std::pair<int,int>(ModuleID,line_count) );
    //cout << "Module "<< ModuleID << " correspond to BD " << line_count << endl;
    line_count++;
  }
  infile.close();
}

void TBReader::Make_dir(){

  cout << "Creating output directories..., will skip if exist..." << endl;
  Read_Module_List();
  
  string outpath = string( dirpath + string("/Module_Ntuple") );

  if( DirectoryExists(outpath.c_str()) ){
    cout << outpath << " has already exist, give up creating other dirs..."
	 << endl;
    return;}

  char command[150];
  int file_check = 0;
  sprintf(command,"mkdir -p %s",outpath.c_str());
  file_check += system(command);

  char dirname[150];
  sprintf(dirname,"%s/TB",outpath.c_str() );
  sprintf(command,"mkdir -p %s",dirname);
  file_check += system(command);

  sprintf(dirname,"%s/Inj",outpath.c_str() );
  sprintf(command,"mkdir -p %s",dirname);
  file_check += system(command);

  string tpro_outpath = string( dirpath + string("/Module_TProfile") );
  sprintf(dirname,"%s",tpro_outpath.c_str() );
  sprintf(command,"mkdir -p %s",dirname);
  file_check += system(command);

  
  cout << "Output directories has been created ... " << endl;
}

void TBReader::Ntuple_Maker(){
  Init();
  
  string outpath = string( dirpath + string("/Module_Ntuple") );
  TFile *outNtuple[MAXBOARDS];
  TTree *outTree[MAXBOARDS];
  TTree *outTree_history[MAXBOARDS];

  Read_Module_List();
  // Variable for output trees
  // vector<unsigned int> *TB_chip;
  // vector<unsigned int> *TB_channel;
  // vector<unsigned int> *TB_type;
  
  for(int ifile = 0 ; ifile < 1 ; ifile++){
    char fpath[200];
    
    // Check root file exist
    int moduleID = Module_List[ifile];
    sprintf(fpath,"%s/TB/Module%d_Oct18.root",outpath.c_str(),moduleID);
    
    ifstream f_check(fpath);
    if( !f_check.good() ){
      outNtuple[ifile] = new TFile(fpath,"recreate");
      outTree[ifile]   = new TTree("TB","TB");

      //TODO: Create outtree branch

      //Fill Run history
      outTree_history[ifile] = new TTree("Run_history","Run_history");
      int Run_hist;
      outTree_history[ifile]->Branch("RunNumber",&Run_hist);
      Run_hist = RunN; 
      outTree_history[ifile]->Fill();
    }
    else{
      outNtuple[ifile] = new TFile(fpath,"update");
      outTree[ifile]   = (TTree*)outNtuple[ifile]->Get("TB");
      //TODO: Set outtree branch Address
      
      // Check and update run_history      
      outTree_history[ifile] = (TTree*)outNtuple[ifile]->Get("Run_history");
      bool already_filled = Check_run_filled(outTree_history[ifile]);
      if(already_filled){
	cout << "Run " << RunN << " has already filled in "
	     << fpath << ", skip it!"<< endl;
	outTree[ifile] = NULL;}
    }
  }
  
  for(int ev = 0 ; ev < nevents ; ++ev){
    T_Rechit->GetEntry(ev);
    T_DWC->GetEntry(ev);
    
    //Event selection
    if(dwcReferenceType != 15) continue;
      
    for(int ifile = 0 ; ifile < 1 ; ifile++){
      if(outTree[ifile] != NULL){	
	// InitoutTreeBranch();
	// SetoutTreeBranchAddress();
	
	outTree[ifile]->Fill();}
    }
  }
  
  for(int ifile = 0 ; ifile < 1 ; ifile++){
    if(outTree[ifile] != NULL){
      outNtuple[ifile]->cd();
      outTree[ifile]->Write("TB",TObject::kOverwrite);
      outTree_history[ifile]->Write("Run_history",TObject::kOverwrite);
    }
    outNtuple[ifile]->Close();
  }
  

}

void TBReader::TProfile_Maker(){
  // TODO: Last Board will has large amount of TProfile, don't know why yet.
  Init();
  //Run selection
  if( PID != 0 ){ cout << "Not e- runs, skip for now." << endl; return; }
    
  string outpath = string( dirpath + string("/Module_TProfile") );
  TFile *outFile[MAXBOARDS];
  TTree *outTree_history[MAXBOARDS];
  int startBD = 0;
  int endBD   = MAXBOARDS;

  Read_Module_List();

  for(int ifile = startBD ; ifile < endBD ; ifile++){
    char fpath[150];
    
    // Check root file exist
    int moduleID = Module_List[ifile];
    sprintf(fpath,"%s/Module%d_Oct18.root",outpath.c_str(),moduleID);
    
    ifstream f_check(fpath);

    if( !f_check.good() ){
      outFile[ifile] = new TFile(fpath,"recreate");
      outTree_history[ifile] = new TTree("Run_history","Run_history");

      //Fill run history
      int Run_hist;
      outTree_history[ifile]->Branch("RunNumber",&Run_hist);
      Run_hist = RunN; 
      outTree_history[ifile]->Fill();
    }
    else{
      outFile[ifile] = new TFile(fpath,"update");

      // Check and update run_history      
      outTree_history[ifile] = (TTree*)outFile[ifile]->Get("Run_history");
      bool already_filled = Check_run_filled(outTree_history[ifile]);
      if(already_filled){
	cout << "Run " << RunN << " has already filled in "
	     << fpath << ", skip it!"<< endl;
	outFile[ifile] = NULL;
      }
    }
  }

  // Set TProfile
  MakePlots M;
  M.Init_Pointers();
  char p_name[200];
  for(int ifile = startBD ; ifile < endBD ; ifile++){
    if(outFile[ifile] == NULL){ continue; }
    // Create if not exist
    sprintf(p_name,"HG_LG_BD%d_chip0_ch0",Module_List[ifile]);
    bool exist = outFile[ifile]->GetListOfKeys()->Contains(p_name);
    int moduleID = Module_List[ifile];
    if(exist){
      for(int chip = 0 ; chip < MAXSKI ; ++chip){
	for(int ch = 0 ; ch < MAXCH ; ++ch){
	  sprintf(p_name,"HG_LG_MID%d_chip%d_ch%d",moduleID,chip,ch*2);
	  M.HG_LG[ifile][chip][ch] = (TProfile*)outFile[ifile]->Get(p_name);
	  sprintf(p_name,"LG_TOT_MID%d_chip%d_ch%d",moduleID,chip,ch*2);
	  M.LG_TOT[ifile][chip][ch] = (TProfile*)outFile[ifile]->Get(p_name);
	}
      }
    }
    //40,54,
    else{
      M.Init_TProfile(ifile);
      // Force write so all channels will exist
      for(int chip = 0 ; chip < MAXSKI ; ++chip){
	for(int ch = 0 ; ch < MAXCH ; ++ch){
	  sprintf(p_name,"HG_LG_chip%d_ch%d",chip,ch*2);
	  M.HG_LG[ifile][chip][ch]->Write(p_name);
	  sprintf(p_name,"LG_TOT_chip%d_ch%d",chip,ch*2);
	  M.LG_TOT[ifile][chip][ch]->Write(p_name);
	}
      }
    }
  }  
  cout << "Looping evts" << endl;
  for(int ev = 0 ; ev < nevents ; ++ev){
    T_Rechit->GetEntry(ev);
    T_DWC->GetEntry(ev);
    
    //Event selection
    if(dwcReferenceType != 15) continue;    

    for(int ihit = 0 ; ihit < NRechits ; ++ihit){
      
      double HG,LG,TOT;
      int chip,ch,BD_order,moduleID;
      HG   = rechit_amplitudeHigh->at(ihit);
      LG   = rechit_amplitudeLow->at(ihit);
      TOT  = rechit_Tot->at(ihit);
      chip = (int)rechit_chip->at(ihit);
      ch   = (int)rechit_channel->at(ihit);
      ch   /= 2;
      moduleID = rechit_module->at(ihit);
      BD_order = moduleID2BDorder.find(moduleID)->second;
      if(outFile[BD_order] == NULL){ continue; }
      if( LG < 5 ) continue;
      M.HG_LG[BD_order][chip][ch]->Fill(LG,HG,1);
      if( TOT < 100 ) continue;
      M.LG_TOT[BD_order][chip][ch]->Fill(TOT,LG,1);      
    }
    
  }
  cout << "End of looping evts" << endl;

  for(int ifile = startBD ; ifile < endBD ; ifile++){
    if( outFile[ifile] == NULL ){ continue; }
    outFile[ifile]->cd();
    
    //cout << outFile[ifile]->GetName() << endl;
    for(int chip = 0 ; chip < MAXSKI ; ++chip){
      for(int ch = 0 ; ch < MAXCH ; ++ch){
	sprintf( p_name,"%s",M.HG_LG[ifile][chip][ch]->GetName() );
	M.HG_LG[ifile][chip][ch]->Write(p_name,TObject::kOverwrite);
	sprintf( p_name,"%s",M.LG_TOT[ifile][chip][ch]->GetName() );
	M.LG_TOT[ifile][chip][ch]->Write(p_name,TObject::kOverwrite);
	delete M.HG_LG [ifile][chip][ch];
	delete M.LG_TOT[ifile][chip][ch];
      }
    }
    outTree_history[ifile]->Write("Run_history",TObject::kOverwrite);
    
    outFile[ifile]->Close();
    delete outFile[ifile];
  }

}


void TBReader::dirty_way(){
  Init();
  //Run selection
  if( PID != 0 ){ cout << "Not e- runs, skip for now." << endl; return; }
  
  Read_Module_List(); 
  int history_run = history->GetEntries();
  bool doublefill = false;
  for(int i = 0 ; i < history_run ; ++i){
    history->GetEntry(i);
    if(RUNNUM == RunN){ doublefill = true; }
  }
  if(doublefill){
    cout << "Run " << RunN << " double filled!" << endl;
    return; }
  else{
    RUNNUM = RunN;
    history->Fill();
    f->cd();
    history->Write("history",TObject::kOverwrite);}

  cout << "Looping evts" << endl;
  for(int ev = 0 ; ev < nevents ; ++ev){
    T_Rechit->GetEntry(ev);
    T_DWC->GetEntry(ev);
    
    //Event selection
    if(dwcReferenceType != 15) continue;    

    for(int ihit = 0 ; ihit < NRechits ; ++ihit){
      
      double HG,LG,TOT;
      int chip,ch,BD_order,moduleID;
      HG   = rechit_amplitudeHigh->at(ihit);
      LG   = rechit_amplitudeLow->at(ihit);
      TOT  = rechit_Tot->at(ihit);
      chip = (int)rechit_chip->at(ihit);
      ch   = (int)rechit_channel->at(ihit);
      ch   /= 2;
      moduleID = rechit_module->at(ihit);
      BD_order = moduleID2BDorder.find(moduleID)->second;
      if( LG < 5 ) continue;
      HGLG[BD_order][chip][ch]->Fill(LG,HG,1);
      if( TOT < 100 ) continue;
      LGTOT[BD_order][chip][ch]->Fill(TOT,LG,1);      
    }    
  }
  cout << "End of looping evts" << endl;

  cout << "Writing TProfiles ..." << endl;
  TDirectory *dir;
  char title[50];
  for(int BD = 0 ; BD < MAXBOARDS ; BD++){
    sprintf(title,"Board_%i",BD);
    dir = (TDirectory*)f->Get(title);
    dir->cd();
    for(int chip = 0 ; chip < MAXSKI ; ++chip){
      for(int ch = 0 ; ch < MAXCH ; ++ch){
	sprintf(title,"HGLG_chip%i_ch%i",chip,ch*2);
	HGLG[BD][chip][ch]->SetTitle(title);
	HGLG[BD][chip][ch]->Write(title,TObject::kOverwrite);
	sprintf(title,"LGTOT_chip%i_ch%i",chip,ch*2);
	LGTOT[BD][chip][ch]->SetTitle(title);
	LGTOT[BD][chip][ch]->Write(title,TObject::kOverwrite);
      }
    }
  }
  cout << "Done!" << endl;

}


bool TBReader::DirectoryExists( const char* pzPath ){
    if ( pzPath == NULL) return false;
 
    DIR *pDir;
    bool bExists = false;
 
    pDir = opendir (pzPath);
 
    if (pDir != NULL)
    {
        bExists = true;    
        (void) closedir (pDir);
    }
 
    return bExists;
}

bool TBReader::Check_run_filled(TTree* tree){
  //Check if this run have already been filled
  int Filled_Runs = tree-> GetEntries();
  int Run_hist;
  tree->SetBranchAddress("RunNumber",&Run_hist);
  bool already_filled = false;
  for(int irun = 0 ; irun < Filled_Runs ; ++ irun){
    tree->GetEntry(irun);
    if(Run_hist == RunN){
      already_filled = true;
    }
  }
  if(!already_filled){
    Run_hist = RunN;
    tree->Fill();
  }
  return already_filled;
}

void dirty_way(){
  char p_name[150];
  ifstream f_check("TPro.root");
  if( !f_check.good() ){
    dirty_way_fexist = true;
    f = new TFile("TPro.root","recreate");
    history = new TTree("Run_history","Run_history");
    history->Branch("Runnum",&RUNNUM);
    int HGLGBIN  = 400;
    int LGTOTBIN = 200;  
    char p_name[200],title[100];
    TDirectory *dir;
    for(int BD = 0 ; BD < MAXBOARDS ; ++BD){
      dir = new TDirectory();
      sprintf(title,"Board_%i",BD);
      dir = f->mkdir(title);
      for(int chip = 0 ; chip < MAXSKI ; ++chip){
	for(int ch = 0 ; ch < MAXCH ; ++ch){	  
	  sprintf(p_name,"HG_LG_BD%d_chip%d_ch%d",BD,chip,ch*2);
	  HGLG[BD][chip][ch] = new TProfile(p_name,"",HGLGBIN,0,800,0,4000);
	  sprintf(p_name,"LG_TOT_BD%d_chip%d_ch%d",BD,chip,ch*2);
	  LGTOT[BD][chip][ch] = new TProfile(p_name,"",LGTOTBIN,0,800,0,2000);
	}
      }
    }
  } 
  else{
    dirty_way_fexist = false;
    f = new TFile("TPro.root","update");
    history = (TTree*) f->Get("Run_history");
    history->SetBranchAddress("Runnum",&RUNNUM);
    for(int BD = 0 ; BD < MAXBOARDS ; ++BD){
      for(int chip = 0 ; chip < MAXSKI ; ++chip){
	for(int ch = 0 ; ch < MAXCH ; ++ch){	  
	  sprintf(p_name,"BD%d/HG_LG_chip%d_ch%d",BD,chip,ch*2);
	  HGLG[BD][chip][ch] = (TProfile*)f->Get(p_name);
	  sprintf(p_name,"BD%d/LG_TOT_chip%d_ch%d",BD,chip,ch*2);
	  LGTOT[BD][chip][ch] = (TProfile*)f->Get(p_name);
	}
      }
    }
  }  
}
