#include "makePlots.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "TApplication.h"
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
#include "fitter.h"
#include <utility>
#include "TProfile.h"

int Hit_counter[28][4][32];
//Constructor
makePlots::makePlots( TChain *c1 ,string filename):T_Rawhit(c1)
{
  cout << "Constructor of makePlot ... \n\n" << endl;
  fname = filename;
}
makePlots::makePlots( TChain *c1,TChain *c2,string filename ):T_Rechit(c1),T_DWC(c2)
{
  cout << "Constructor of makePlot ... \n\n" << endl;
  fname = filename;
  // Set object pointer(Data)
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

//Destructor
makePlots::~makePlots()
{
  root_out->Close();
  cout << "\n\n";
  cout << "Destructor of makePlot ... " << endl;
}

void makePlots::begin(){
  ofstream of(string(dirpath+string("LG_US.txt")).c_str());
  of << "BD\tchip\tCH\n";
  of.close();
  root_out = new TFile(string(dirpath+string("TPro.root")).c_str(),"recreate");

}


bool makePlots::DBG(){
  return false;
}

void makePlots::Init(){

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
void makePlots::Init_BeamE(){
  /*
    string beamE_str;
    int end = fname.find("GeV");
    int par;
    int par_e  = fname.find("Ele");
    int par_pi = fname.find("Pi");
    if(par_e != -1){
      PID = 0;
      par = par_e;
      beam_str = "Ele";
      beamE_str = fname.substr(par+3,end-par-3);}
    else if( par_pi != -1 ){
      PID = 1;
      par = par_pi;
      beam_str = "Pi";      
      beamE_str = fname.substr(par+2,end-par-2);}
    beamE = stoi(beamE_str);}
  cout << beam_str.c_str()  << " , "<< beamE << "GeV\n" << endl;
  */
}
void makePlots::Loop(){
  root_logon();
  begin();
  Init();
  Draw_HG_LG();
  
}
void makePlots::Draw_HG_LG(){
  
  int MAXBD   = 28;
  int MAXCHIP = 4;
  int MAXCH   = 32;
  char title[50];
  TProfile *tpr_HGLG[MAXBD][MAXCHIP][MAXCH];
  TProfile *tpr_LGTOT[MAXBD][MAXCHIP][MAXCH];
  int tpr_LS[MAXBD][MAXCHIP][MAXCH];
  h_tprLGUS = h_tprLGUS = new TH1D("LG_US_tpr","LG_US_tpr",50,0,10);
  
  for(int BD = 0; BD < MAXBD ; ++BD){
    for(int chip = 0 ; chip < MAXCHIP ; ++chip){
      for(int ch = 0 ; ch < MAXCH ;++ch){
	sprintf(title,"HGLG_BD%i_chip%i_ch%i",BD,chip,ch*2);
	tpr_HGLG[BD][chip][ch] = new TProfile(title,title,400,0,800,0,4000);
	sprintf(title,"LGTOT_BD%i_chip%i_ch%i",BD,chip,ch*2);
	tpr_LGTOT[BD][chip][ch] = new TProfile(title,title,300,0,800,0,3000);
	tpr_LS[BD][chip][ch] = 0;      }}}
  
  for(int ev = 0 ; ev < nevents ; ++ev){
    T_Rechit->GetEntry(ev);
    if(ev % 10000 == 0)
      cout << "Processing " << ev << " / " << nevents << " ..." << endl;
    //T_DWC   ->GetEntry(ev);
    
    for(int hit = 0 ; hit < (int) rechit_amplitudeHigh->size() ; ++hit){
      
      double HG,LG,TOT;
      int chip,ch,BD;
      HG   = rechit_amplitudeHigh->at(hit);
      LG   = rechit_amplitudeLow->at(hit);
      TOT  = rechit_Tot->at(hit);
      chip = (int)rechit_chip->at(hit);
      ch   = (int)rechit_channel->at(hit);
      ch   /= 2;
      BD   = rechit_layer->at(hit);
      BD   -= 1;
      if( HG > 200 && LG < 20) tpr_LS[BD][chip][ch]++;
      if( LG < 5 ) continue;
      tpr_HGLG[BD][chip][ch]->Fill(LG,HG,1);
      tpr_LGTOT[BD][chip][ch]->Fill(TOT,LG,1);
      
    }
  }

  TDirectory *dir;
  for(int BD = 0; BD < MAXBD ; ++BD){
    dir = new TDirectory();
    sprintf(title,"Board_%i",BD);
    dir = root_out->mkdir(title);
    dir->cd();
    for(int chip = 0 ; chip < MAXCHIP ; ++chip){
      for(int ch = 0 ; ch < MAXCH ;++ch){
	if(tpr_HGLG[BD][chip][ch]->GetEntries() == 0){
	  continue;}	
	sprintf(title,"HGLG_chip%i_ch%i",chip,ch*2);
	tpr_HGLG[BD][chip][ch]->SetTitle(title);
	tpr_HGLG[BD][chip][ch]->SetName(title);
	tpr_HGLG[BD][chip][ch]->SetMarkerStyle(20);
	tpr_HGLG[BD][chip][ch]->SetMarkerSize(1.2);
	tpr_HGLG[BD][chip][ch]->SetMarkerColor(chip+1);
	tpr_HGLG[BD][chip][ch]->Write(title);

	if(tpr_LGTOT[BD][chip][ch]->GetEntries() == 0){
	  continue;}
	sprintf(title,"LGTOT_chip%i_ch%i",chip,ch*2);
	tpr_LGTOT[BD][chip][ch]->SetTitle(title);
	tpr_LGTOT[BD][chip][ch]->SetName(title);
	tpr_LGTOT[BD][chip][ch]->SetMarkerStyle(20);
	tpr_LGTOT[BD][chip][ch]->SetMarkerSize(1.2);
	tpr_LGTOT[BD][chip][ch]->SetMarkerColor(chip+1);
	tpr_LGTOT[BD][chip][ch]->Write(title);

	
	if(tpr_LS[BD][chip][ch] != 0){
	  double US_LG_percernt = tpr_LS[BD][chip][ch]*100./tpr_HGLG[BD][chip][ch]->GetEntries();
	  h_tprLGUS->Fill(US_LG_percernt);
	  ofstream of(string(dirpath+string("LG_US.txt")).c_str(), std::ios::app);
	  of << BD << "\t" << chip << "\t" << ch*2 << "\t" << US_LG_percernt
	     << endl;
	     of.close();}	
      }
    }	
  }
  
  
  for(int BD = 0; BD < MAXBD ; ++BD){
    for(int chip = 0 ; chip < MAXCHIP ; ++chip){
      for(int ch = 0 ; ch < MAXCH ;++ch){
	delete tpr_HGLG[BD][chip][ch];
	delete tpr_LGTOT[BD][chip][ch];}}}
  
} 

void makePlots::InitTH2Poly(TH2Poly& poly)
{
  int MAXVERTICES = 6;
  double HexX[MAXVERTICES];
  double HexY[MAXVERTICES];
  int iu,iv,CellXYsize;
  ifstream file("src_txtfile/poly_frame.txt");
  string line;

  
  for(int header = 0; header < 4; ++header )     getline(file,line);
  
  while(true){
    getline(file,line);
    if( file.eof() ) break;
    file >> iu >> iv >> CellXYsize;    
    for(int i = 0; i < CellXYsize ; ++i){
      getline(file,line);
      file >> HexX[i] >> HexY[i];
    }
    
    poly.AddBin(CellXYsize, HexX, HexY);
  }
  file.close();

}

void makePlots::root_logon(){

cout << endl << "Welcome to the ATLAS rootlogon.C" << endl;
//
// based on a style file from BaBar
//

//..BABAR style from RooLogon.C in workdir
TStyle *atlasStyle= new TStyle("ATLAS","Atlas style");

// use plain black on white colors
 Int_t icol=0;
atlasStyle->SetFrameBorderMode(icol);
atlasStyle->SetCanvasBorderMode(icol);
atlasStyle->SetPadBorderMode(icol);
atlasStyle->SetPadColor(icol);
atlasStyle->SetCanvasColor(icol);
atlasStyle->SetStatColor(icol);
//atlasStyle->SetFillColor(icol);

// set the paper & margin sizes
atlasStyle->SetPaperSize(20,26);
atlasStyle->SetPadTopMargin(0.1);
//atlasStyle->SetPadRightMargin(0.05);
atlasStyle->SetPadRightMargin(0.12);
atlasStyle->SetPadBottomMargin(0.16);
atlasStyle->SetPadLeftMargin(0.12);

// use large fonts
//Int_t font=72;
Int_t font=32;
Double_t tsize=0.05;
atlasStyle->SetTextFont(font);


atlasStyle->SetTextSize(tsize);
atlasStyle->SetLabelFont(font,"x");
atlasStyle->SetTitleFont(font,"x");
atlasStyle->SetLabelFont(font,"y");
atlasStyle->SetTitleFont(font,"y");
atlasStyle->SetLabelFont(font,"z");
atlasStyle->SetTitleFont(font,"z");

atlasStyle->SetLabelSize(tsize,"x");
atlasStyle->SetTitleSize(tsize,"x");
atlasStyle->SetLabelSize(tsize,"y");
atlasStyle->SetTitleSize(tsize,"y");
atlasStyle->SetLabelSize(tsize,"z");
atlasStyle->SetTitleSize(tsize,"z");
//atlasStyle->SetTitleOffset(1.2,"y");

//use bold lines and markers
atlasStyle->SetMarkerStyle(20);
atlasStyle->SetMarkerSize(1.2);
atlasStyle->SetHistLineWidth(2.);
atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

//get rid of X error bars and y error bar caps
//atlasStyle->SetErrorX(0.001);

//do not display any of the standard histogram decorations
//atlasStyle->SetOptTitle(0);
//atlasStyle->SetOptStat(1111);
atlasStyle->SetOptStat(0);
//atlasStyle->SetOptFit(1111);
atlasStyle->SetOptFit(0);

// put tick marks on top and RHS of plots
atlasStyle->SetPadTickX(1);
atlasStyle->SetPadTickY(1);
 

gROOT->SetStyle("Plain");

//gStyle->SetPadTickX(1);
//gStyle->SetPadTickY(1);

}
