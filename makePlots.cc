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
  ofstream of("LG_US.txt");
  of << "BD\tchip\tCH\n";
  of.close();
  root_out = new TFile("TPro.root","recreate");
}

//Destructor
makePlots::~makePlots()
{
  root_out->Close();
  cout << "\n\n";
  cout << "Destructor of makePlot ... " << endl;
}

bool makePlots::DBG(){
  return false;
}

void makePlots::Init(){
  T_Rawhit->SetBranchAddress("eventID", &eventID );
  T_Rawhit->SetBranchAddress("skirocID", &skirocID);
  T_Rawhit->SetBranchAddress("boardID", &boardID);
  T_Rawhit->SetBranchAddress("moduleID", &moduleID);
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
  int NLAYER = 28;
  root_logon();
  Init();
  gROOT->SetBatch(kTRUE);
  nevents = T_Rawhit->GetEntries();
  c1 = new TCanvas();
  
  
  for(int l = 0; l < NLAYER ; ++l){
    for(int chip = 0 ; chip < 4 ; ++chip){
      for(int Nch = 0 ; Nch < 32 ; ++Nch){
	Hit_counter[l][chip][Nch] = 0;
      }
    }
  }
  cout << "Looping over all hits..." << endl;
  int count_0 = 0;
  int count_1 = 0;
  
  for(int ev = 0; ev < nevents; ++ev){
    T_Rawhit->GetEntry(ev);
    Hit_counter[boardID][skirocID][channelID/2]++;
    if(HighGainStatus == 1) count_1++;
    else count_0++;
  }
  
  cout << "count 0: " << count_0 << ", count1: " << count_1 << endl;
  //getchar();
  h_LGundershoot = new TH1D("LG_undershoot","LG_undershoot",40,0,20);
  h_tprLGUS = new TH1D("LG_US_tpr","LG_US_tpr",50,0,10);
  
  for(int i = 0 ; i < NLAYER; ++i)
    Draw_HG_LG(i);

  h_LGundershoot->Draw();
  c1->Update();
  c1->SaveAs("LG_US.png");
  //getchar();
  TGraph* gr = 0;
  if(gr != NULL){
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(2);
    gr->GetYaxis()->SetTitle("HG(ADC)");
    gr->GetYaxis()->SetTitleOffset(1.2);
    gr->GetXaxis()->SetTitle("LG(ADC)");
    gr->Draw("AP");
    c1->Update();
    //getchar();
  }
  else
    cout << "haha" << endl;
  root_out->cd();
  h_tprLGUS->Write("tprLGUS");
}

TGraph* makePlots::Draw_HG_LG(int BD = 0){
  bool save_png = false;
  //if(Hit_counter[BD][chip][ch] < 1000) return NULL;
  int NCHIP = 4;
  int NCH   = 32;
  vector< vector< vector< double > > > HG_vec,LG_vec,TOT_vec;
  for(int chip = 0 ; chip < NCHIP ; ++chip){
    HG_vec.resize(NCHIP);
    LG_vec.resize(NCHIP);
    TOT_vec.resize(NCHIP);
    for(int ch = 0; ch < NCH ; ++ch){
      HG_vec[chip].resize(NCH);
      LG_vec[chip].resize(NCH);
      TOT_vec[chip].resize(NCH);
    }
  }
  for(int ev = 0; ev < nevents; ++ev){
    if(ev % 100000 == 0 && DBG()) cout << "processing hit " << ev << endl;
    T_Rawhit->GetEntry(ev);
    if(boardID != BD) continue;
      HG_vec[skirocID][channelID/2].push_back(HighGainADC);
      LG_vec[skirocID][channelID/2].push_back(LowGainADC);
      TOT_vec[skirocID][channelID/2].push_back(TotSlow);}
  //}
  
  TGraph  *gr;
  TMultiGraph  *mgr;
  TLegend *leg;

  char title[50];
  char title_sub[20];
  /*
  for(int ch = 0; ch < NCH ; ++ch){
    mgr = new TMultiGraph();
    leg = new TLegend(0.65,0.13,0.9,0.4);
    leg->SetBorderSize(0);

    if( Hit_counter[BD][0][ch] < 500 ||  Hit_counter[BD][1][ch] < 500) continue;
    if( Hit_counter[BD][2][ch] < 500 ||  Hit_counter[BD][3][ch] < 500) continue;

    for(int chip = 0 ; chip < NCHIP ; ++chip){
 
      gr = new TGraph(HG_vec[chip][ch].size(),&LG_vec[chip][ch][0],&HG_vec[chip][ch][0]);
      //cout << HG_vec[chip][ch].size() << endl;
      //cout << LG_vec[chip][ch].size() << endl;
      
      fitter f(gr,HG_vec[chip][ch],LG_vec[chip][ch],TOT_vec[chip][ch]);
      f.fit_Graph();
      h_LGundershoot->Fill(f.undershoot_percent);
      gr->Draw("AP");
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(0.2);
      gr->SetMarkerColor(chip+1);
      
      sprintf(title_sub,"chip%i",chip);
      leg->AddEntry(gr,title_sub,"P");
      mgr->Add(gr);
      
      if(f.undershoot_percent > 2){
	ofstream of("LG_US.txt", std::ios::app);
	of << "BD ,chip, CH " << BD << "," << chip << "," << ch*2
	   << " LG_US = " << f.undershoot_percent << endl;
	of.close();}
      
    }

    mgr->Draw("AP");
    mgr->GetYaxis()->SetTitle("HG(ADC)");
    mgr->GetYaxis()->SetTitleOffset(1.2);
    mgr->GetXaxis()->SetTitle("LG(ADC)");
    
    sprintf(title,"Board_%iCH%i",BD,ch*2);
    leg->SetHeader(title);
    leg->Draw("same");
    c1->Update();
    sprintf(title,"plot_out/BD_%iCH%i.png",BD,ch*2);
    if(save_png)
      c1->SaveAs(title);
  }
  */    
  TProfile *tpr_HGLG;
  TProfile *tpr_LGTOT;
  
  //root_out->cd();
  TDirectory *dir = new TDirectory();
  sprintf(title,"Board_%i",BD);
  dir = root_out->mkdir(title);
  dir->cd();
  for(int chip = 0 ; chip < NCHIP ; ++chip){
    for(int ch = 0 ; ch < NCH ;++ch){
      if(HG_vec[chip][ch].size() < 1000) continue;
      tpr_HGLG  = new TProfile("Tpro","Tpro",40,0,800,0,4000);
      tpr_LGTOT = new TProfile("Tpro2","Tpro2",40,0,800,0,4000);
      
      int tpr_LS = 0;
      for(size_t i = 0; i < HG_vec[chip][ch].size(); ++i){
	if(LG_vec[chip][ch].at(i) < 20 && HG_vec[chip][ch].at(i) > 200)
	  tpr_LS++;
	if(LG_vec[chip][ch].at(i) < 20) continue;
	
        tpr_HGLG->Fill(LG_vec[chip][ch].at(i),HG_vec[chip][ch].at(i),1);
	tpr_LGTOT->Fill(TOT_vec[chip][ch].at(i),LG_vec[chip][ch].at(i));
      }
      
      tpr_HGLG->SetMarkerStyle(20);
      tpr_HGLG->SetMarkerSize(1.2);
      tpr_HGLG->SetMarkerColor(chip+1);
      
      tpr_HGLG->Draw("e");
      c1->Update();
      
      tpr_HGLG->Draw("esame");
      c1->Update();
      sprintf(title,"HGLG_chip%i,ch%i",chip,ch*2);
      tpr_HGLG->Write(title);
      sprintf(title,"LGTOT_chip%i,ch%i",chip,ch*2);
      tpr_LGTOT->Write(title);
      
      delete tpr_HGLG;
      delete tpr_LGTOT;
      //getchar();
      if(tpr_LS != 0){
	//cout << tpr_LS*100./LG_vec[chip][ch].size() << endl;
	h_tprLGUS->Fill(tpr_LS*100./LG_vec[chip][ch].size());
	ofstream of("LG_US.txt", std::ios::app);
	of << BD << "\t" << chip << "\t" << ch*2
	   << "\t" <<  tpr_LS*100./LG_vec[chip][ch].size() << endl;
	of.close();
      }
    }
  }
    
  return NULL;
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
