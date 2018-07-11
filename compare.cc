#include "compare.h"
#include "TLegend.h"
#include "TTree.h"
#include "TROOT.h"
#include <vector>
#include <array>


compare::compare(){
  c1 = new TCanvas();
}
compare::~compare(){}

void compare::compare_Ene(int method){

  if(method != 1) return;
  int Ene[6] = {10 , 30 ,50 ,80 ,100 ,150};
  for(int i = 0 ; i < 6 ; ++i){
    sprintf(title,"Linear_result/HGLG_sat_%iGeV.root",Ene[i]);
    store_GR(string(title),i);}
  
  for(int BD = 0; BD < MAXBD ; ++BD){
    TMultiGraph *mgr = new TMultiGraph();
    TLegend *leg = new TLegend(0.65,0.65,0.87,0.87);
    leg->SetBorderSize(0);
    for(int lab = 0 ; lab < MAXLABEL ; ++lab){
      if(gr_p1[lab][BD] == NULL) continue;
      gr_p1[lab][BD]->SetMarkerStyle(20);
      gr_p1[lab][BD]->SetMarkerSize(0.8);
      gr_p1[lab][BD]->SetMarkerColor(lab+1);
      mgr->Add(gr_p1[lab][BD],"P");
      sprintf(title,"Energy_%iGeV",Ene[lab%6]);
      leg->AddEntry(gr_p1[lab][BD],title,"P");
    }
    
    string titlex("skiroc*32+channel/2");
    string titley("p1");
    Set_mgr(*mgr,titlex,titley);

    mgr->SetMaximum(10);
    mgr->Draw("sameAP");
    
    sprintf(title,"Board%i",BD);
    leg->SetHeader(title,"C");
    leg->Draw("same");
    c1->Update();
    getchar();  
  }
}

void compare::store_GR(string fname,int label){
  if( label >= MAXLABEL || label < 0 ){
    cout << "store_GR() label exceed!" << endl;
    return ;}
  
  TFile f(fname.c_str());
  if( f.IsZombie() ) {
    cout << fname << " : no such file exist!" << endl;
    return;  }
  else{
    cout << "open file: " << fname << endl;
  }
  TTree *tt = (TTree*) f.Get("tree");
  int layerID          = 0;
  int skirocID         = 0;
  int channelID        = 0;
  double p0            = 0;
  double p1            = 0;
  double p_saturation  = 0;
  bool good_saturation = 0;
  
  tt->SetBranchAddress("layerID",&layerID);
  tt->SetBranchAddress("skirocID",&skirocID);
  tt->SetBranchAddress("channelID",&channelID);
  tt->SetBranchAddress("p0",&p0);
  tt->SetBranchAddress("p1",&p1);
  tt->SetBranchAddress("p_saturation",&p_saturation);
  tt->SetBranchAddress("good_saturation",&good_saturation);
  
  int nentry = tt->GetEntries();

  array< vector<double>, MAXBD> p0_vec,p1_vec,sat_vec,X_vec;
  for(int data = 0 ; data < nentry ; ++data){
    tt->GetEntry(data);
    if(!good_saturation) continue;
    p0_vec[layerID].push_back(p0);
    p1_vec[layerID].push_back(p1);
    sat_vec[layerID].push_back(p_saturation);
    X_vec[layerID].push_back(skirocID*32+channelID/2);    }

  for(int BD = 0 ; BD < MAXBD ; ++BD){
    gr_p0[label][BD]  = NULL;
    gr_p1[label][BD]  = NULL;
    gr_sat[label][BD] = NULL;

    if( X_vec[BD].size() != 0 ){
      gr_p0[label][BD] =
	new TGraph(X_vec[BD].size(),&X_vec[BD][0],&p0_vec[BD][0]);
      gr_p1[label][BD] =
	new TGraph(X_vec[BD].size(),&X_vec[BD][0],&p1_vec[BD][0]);
      gr_sat[label][BD] =
	new TGraph(X_vec[BD].size(),&X_vec[BD][0],&sat_vec[BD][0]);
    }
  }
  
}

void compare::Set_mgr(TMultiGraph& mgr,string xtitle,string ytitle){
  mgr.Draw("AP");
  TAxis *ax = mgr.GetXaxis();
  TAxis *ay = mgr.GetYaxis();
  
  ax->SetTitle(xtitle.c_str());
  ay->SetTitle(ytitle.c_str());
  double tsize = 0.05;
  ax->SetLabelSize(tsize);
  ax->SetTitleSize(tsize);    
  ax->SetLabelFont(32);
  ax->SetTitleFont(32);
  ax->SetTitleOffset(0.8);
  ay->SetLabelSize(tsize);
  ay->SetTitleSize(tsize);
  ay->SetLabelFont(32);
  ay->SetTitleFont(32);
  ay->SetTitleOffset(0.8);
}
