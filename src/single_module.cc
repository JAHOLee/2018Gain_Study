#include "single_module.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include <utility>
#include "TProfile.h"
#include <sstream>

single_module::single_module( TChain *chain, string filename ):T_Rawhit(chain)
{
  cout << "Constructor of makePlot ... \n\n" << endl;

  //TFile f(filename.c_str());
  //T_Rawhit = (TTree*)f.Get("pulseshapeplotter/tree");
   
  root_out = new TFile("inj.root","update");
  if(root_out->IsZombie())
    root_out = new TFile("Inj.root","recreate");    
  c1 = new TCanvas();
  fname = filename;
}

//Destructor
single_module::~single_module()
{
  root_out->Close();
  delete c1;
  cout << "\n\n";
  cout << "Destructor of makePlot ... " << endl;
}

void single_module::Init(){
  
   skirocID = 0;
   boardID = 0;
   channelID = 0;
   HighGainADC = 0;
   HighGainTmax = 0;
   HighGainChi2 = 0;
   HighGainErrorADC = 0;
   HighGainErrorTmax = 0;
   HighGainStatus = 0;
   HighGainNCalls = 0;
   LowGainADC = 0;
   LowGainTmax = 0;
   LowGainChi2 = 0;
   LowGainErrorADC = 0;
   LowGainErrorTmax = 0;
   LowGainStatus = 0;
   LowGainNCalls = 0;
   TotSlow = 0;
   ToaRise = 0;
   ToaFall = 0;

  
   T_Rawhit->SetBranchAddress("eventID", &eventID);
   T_Rawhit->SetBranchAddress("skirocID", &skirocID);
   T_Rawhit->SetBranchAddress("boardID", &boardID);
   T_Rawhit->SetBranchAddress("channelID", &channelID);
   T_Rawhit->SetBranchAddress("HighGainADC", &HighGainADC);
   T_Rawhit->SetBranchAddress("HighGainTmax", &HighGainTmax);
   T_Rawhit->SetBranchAddress("HighGainChi2", &HighGainChi2);
   T_Rawhit->SetBranchAddress("HighGainErrorADC", &HighGainErrorADC);
   T_Rawhit->SetBranchAddress("HighGainErrorTmax", &HighGainErrorTmax);
   T_Rawhit->SetBranchAddress("HighGainStatus", &HighGainStatus);
   T_Rawhit->SetBranchAddress("HighGainNCalls", &HighGainNCalls);
   T_Rawhit->SetBranchAddress("LowGainADC", &LowGainADC);
   T_Rawhit->SetBranchAddress("LowGainTmax", &LowGainTmax);
   T_Rawhit->SetBranchAddress("LowGainChi2", &LowGainChi2);
   T_Rawhit->SetBranchAddress("LowGainErrorADC", &LowGainErrorADC);
   T_Rawhit->SetBranchAddress("LowGainErrorTmax", &LowGainErrorTmax);
   T_Rawhit->SetBranchAddress("LowGainStatus", &LowGainStatus);
   T_Rawhit->SetBranchAddress("LowGainNCalls", &LowGainNCalls);
   T_Rawhit->SetBranchAddress("TotSlow", &TotSlow);
   T_Rawhit->SetBranchAddress("ToaRise", &ToaRise);
   T_Rawhit->SetBranchAddress("ToaFall", &ToaFall);

}

void single_module::Loop(){
  Setname();
  if(!inj_sweep) {
    cout << "single_module::Loop only deal with sweep injection run!" << endl;
    return;}
  
  Init();
  //gROOT->SetBatch(kTRUE);
  nevents = T_Rawhit->GetEntries();
  if( nevents != inj_event ){
    cout << nevents << " , " << inj_event << endl;
    cout << "yaml events not match! skip!" << endl;
    return;}
  Fill_Tprofile();
}

void single_module::Fill_Tprofile(){
  //Assume single channel injection
  int MAXCHIP = 4;
  char title[50];

  int layer_to_moduleID[28] = { 78, 90, 89, 88, 77,
				85, 84, 32, 69, 79,
				67, 65, 76, 83, 35,
				36, 70, 73, 44, 51,
				86, 87, 54, 62, 64,
				55, 59, 71 };
  

  string  moduleID = moduleID_str.substr(6);
  int moduleID_int = atoi( moduleID.c_str() );
  int BD_layer = -1;

  for( int i = 0 ; i < 28; ++i ){
    if(layer_to_moduleID[i] ==  moduleID_int)
      BD_layer = i;  }
  if(BD_layer == -1){
    cout << moduleID_str << " not used in June TB!" << endl;
    return;  }
  
}

void single_module::Setname(){
  int module_start = fname.find("module");
  int module_end   = fname.find("/",module_start+1);
  int lastslash    = fname.find_last_of("/");
  int findroot     = fname.find(".root");
  moduleID_str = fname.substr(module_start,module_end-module_start);
  labelID   = fname.substr(lastslash+1,findroot-lastslash-1);
  filepath  = fname.substr(0,module_end+1);
  if((int)labelID.find("_pedestal") != -1){
    labelID = labelID.substr(0,(int)labelID.length()-9);}
  cout << "module : " << moduleID_str << "\nlabel : " << labelID << endl;
  cout << "path : " << filepath << endl;
  string yaml;
  yaml.append(filepath);
  yaml.append("yaml/");
  yaml.append(labelID);
  yaml.append(".yaml");
  Read_yaml(yaml);

}

void single_module::Read_yaml(string yaml){
  cout << "yaml file: " << yaml << endl;

  ifstream yaml_in(yaml);
  if(!yaml_in.is_open()){
    cout << "can't find yaml file " << yaml << endl;
    return;}

  string line;
  int line_label = 0;

  while(true){
    getline(yaml_in,line);
    if(yaml_in.eof()) break;

    if(line_label == 1){
      if((int)line.find("sweep") != -1)
	inj_sweep = true;
      else{
	inj_sweep = false;}    }

    if(line_label == 2){
      string before_str = "channelIds: [";
      int start = line.find(before_str);
      int end   = line.find("]");
      inj_CH_str = line.substr(start+before_str.length(),end - before_str.length() - start);
      if(inj_CH_str.length() != 0){
	istringstream iss(inj_CH_str);
	string token;
	while(getline(iss,token,',')){
	  inj_CH.push_back(stoi(token));
	}
      }
    }

    if(line_label == 9){
      string before_str = "nEvent: ";
      int start = line.find(before_str);
      string tmp_str = line.substr(start+before_str.length());
      inj_event = atoi( tmp_str.c_str() );
    }
      
    //cout << " Line: " << line_label << ", " << line << endl;
    line_label++;
  }

  cout << "type: " << inj_sweep << ", CH: ";
  for(int i = 0 ; i < (int)inj_CH.size() ; ++i){
    cout << inj_CH[i] << " ";
  }
  cout << ", evt = " << inj_event << endl;

}

