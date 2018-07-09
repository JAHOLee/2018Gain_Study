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
  
   skirocID = 0;
   boardID = 0;
   moduleID = 0;
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
  begin();
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
  cout << "Looping over all hits...(for counting)" << endl;

  int count_0 = 0;
  int count_1 = 0;
  int m_BDID,m_SKIID,m_CHID;
  
  for(int ev = 0; ev < nevents; ++ev){
    if(ev % 10000 == 0) cout << "processing evt " << ev << endl;
    T_Rawhit->GetEntry(ev);
    for(int hit = 0 ; hit < (int)boardID->size() ; ++hit){
      m_BDID  = boardID->at(hit);
      m_SKIID = skirocID->at(hit);
      m_CHID  = channelID->at(hit)/2;
      Hit_counter[m_BDID][m_SKIID][m_CHID]++;
      if(HighGainStatus->at(hit) == 1) count_1++;
      else count_0++;}
  }
  
  cout << "count 0: " << count_0 << ", count1: " << count_1 << endl;
  //getchar();
  h_LGundershoot = new TH1D("LG_undershoot","LG_undershoot",40,0,20);
  h_tprLGUS = new TH1D("LG_US_tpr","LG_US_tpr",50,0,10);
  
  //for(int i = 0 ; i < NLAYER; ++i)
    //Draw_HG_LG(i,1);
  Draw_HG_LG();
  h_LGundershoot->Draw();
  c1->Update();
  c1->SaveAs(string(dirpath+string("LG_US.png")).c_str());
  
  root_out->cd();
  h_tprLGUS->Write("tprLGUS");
}
void makePlots::Draw_HG_LG(){
  int MAXBD   = 28;
  int MAXCHIP = 4;
  int MAXCH   = 32;
  char title[50];
  TProfile *tpr_HGLG[MAXBD][MAXCHIP][MAXCH];
  TProfile *tpr_LGTOT[MAXBD][MAXCHIP][MAXCH];
  int tpr_LS[MAXBD][MAXCHIP][MAXCH];
  
  
  for(int BD = 0; BD < MAXBD ; ++BD){
    for(int chip = 0 ; chip < MAXCHIP ; ++chip){
      for(int ch = 0 ; ch < MAXCH ;++ch){
	sprintf(title,"HGLG_BD%i_chip%i_ch%i",BD,chip,ch*2);
	tpr_HGLG[BD][chip][ch] = new TProfile(title,title,400,0,800,0,4000);
	sprintf(title,"LGTOT_BD%i_chip%i_ch%i",BD,chip,ch*2);
	tpr_LGTOT[BD][chip][ch] = new TProfile(title,title,300,0,800,0,3000);
	tpr_LS[BD][chip][ch] = 0;      }}}
  
  for(int ev = 0 ; ev < nevents ; ++ev){
    T_Rawhit->GetEntry(ev);
    for(int hit = 0 ; hit < (int) HighGainADC->size() ; ++hit){
      double HG,LG,TOT;
      int chip,ch,BD;
      HG   = HighGainADC->at(hit);
      LG   = LowGainADC->at(hit);
      TOT  = TotSlow->at(hit);
      chip = skirocID->at(hit);
      ch   = channelID->at(hit);
      ch   /= 2;
      BD   = boardID->at(hit);
      if( HG > 200 && LG < 20) tpr_LS[BD][chip][ch]++;
      if( LG < 20 ) continue;
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
	sprintf(title,"HGLG_chip%i_ch%i",chip,ch*2);
	tpr_HGLG[BD][chip][ch]->SetTitle(title);
	tpr_HGLG[BD][chip][ch]->SetName(title);
	tpr_HGLG[BD][chip][ch]->SetMarkerStyle(20);
	tpr_HGLG[BD][chip][ch]->SetMarkerSize(1.2);
	tpr_HGLG[BD][chip][ch]->SetMarkerColor(chip+1);
	tpr_HGLG[BD][chip][ch]->Write(title);
	sprintf(title,"LGTOT_chip%i_ch%i",chip,ch*2);
	tpr_LGTOT[BD][chip][ch]->SetTitle(title);
	tpr_LGTOT[BD][chip][ch]->SetName(title);
	tpr_LGTOT[BD][chip][ch]->SetMarkerStyle(20);
	tpr_LGTOT[BD][chip][ch]->SetMarkerSize(1.2);
	tpr_LGTOT[BD][chip][ch]->SetMarkerColor(chip+1);
	tpr_LGTOT[BD][chip][ch]->Write(title);
  
	if(tpr_LS[BD][chip][ch] != 0){
	  double US_LG_percernt = tpr_LS[BD][chip][ch]*100./Hit_counter[BD][chip][ch];
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

void makePlots::Draw_HG_LG(int BD = 0,bool Draw_SCAT = 0){
  bool save_png = false;
  //if(Hit_counter[BD][chip][ch] < 1000) return NULL;
  int NCHIP = 4;
  int NCH   = 32;
  char title[50];
  char title_sub[20];

  
   vector< vector< vector< double > > > HG_vec,LG_vec,TOT_vec;
  // for(int chip = 0 ; chip < NCHIP ; ++chip){
  //   HG_vec.resize(NCHIP);
  //   LG_vec.resize(NCHIP);
  //   TOT_vec.resize(NCHIP);
  //   for(int ch = 0; ch < NCH ; ++ch){
  //     HG_vec[chip].resize(NCH);
  //     LG_vec[chip].resize(NCH);
  //     TOT_vec[chip].resize(NCH);
  //   }
  // }
   cout << "staring Board " << BD << endl;

  for(int ev = 0; ev < nevents; ++ev){
    //if(ev % 10000 == 0) cout << "processing evt " << ev << endl;
    cout << "processing evt " << ev << endl;
    T_Rawhit->GetEntry(ev);
    //for(int hit = 0 ;hit < (int)channelID->size(); ++hit ){
      //if(boardID->at(hit) != BD) continue;
      //cout << skirocID->at(hit) << "," << channelID->at(hit)/2 << "," << HighGainADC->at(hit);
      //HG_vec[skirocID->at(hit)][channelID->at(hit)/2].push_back(HighGainADC->at(hit));
      //LG_vec[skirocID->at(hit)][channelID->at(hit)/2].push_back(LowGainADC->at(hit));
      //TOT_vec[skirocID->at(hit)][channelID->at(hit)/2].push_back(TotSlow->at(hit));
      cout << "finish pushing " << ev << endl;
  }

    //}
  if(Draw_SCAT){  

    TGraph  *gr;
    TMultiGraph  *mgr;
    TLegend *leg;
    for(int ch = 0; ch < NCH ; ++ch){
      mgr = new TMultiGraph();
      leg = new TLegend(0.65,0.13,0.9,0.4);
      leg->SetBorderSize(0);

      for(int chip = 0 ; chip < NCHIP ; ++chip){
	if( Hit_counter[BD][chip][ch] < 1000 ) continue;

	
	gr = new TGraph(HG_vec[chip][ch].size(),&LG_vec[chip][ch][0],&HG_vec[chip][ch][0]);
	//cout << HG_vec[chip][ch].size() << endl;
	//cout << LG_vec[chip][ch].size() << endl;
      
	//fitter f(gr,HG_vec[chip][ch],LG_vec[chip][ch],TOT_vec[chip][ch]);
	//f.fit_Graph();
	//h_LGundershoot->Fill(f.undershoot_percent);
	gr->Draw("AP");
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.2);
	gr->SetMarkerColor(chip+1);
      
	sprintf(title_sub,"chip%i",chip);
	leg->AddEntry(gr,title_sub,"P");
	mgr->Add(gr);
      } 
      mgr->Draw("AP");
      mgr->GetYaxis()->SetTitle("HG(ADC)");
      mgr->GetYaxis()->SetTitleOffset(1.2);
      mgr->GetXaxis()->SetTitle("LG(ADC)");
    
      sprintf(title,"Board_%iCH%i",BD,ch*2);
      leg->SetHeader(title);
      leg->Draw("same");
      c1->Update();
      sprintf(title,"~/HG_LG/plot_out/BD_%iCH%i.png",BD,ch*2);
      if(save_png)
	c1->SaveAs(title);
    }
  }
  

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
