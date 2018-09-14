#include "compare.h"
#include "TLegend.h"
#include "TTree.h"
#include "TROOT.h"
#include "TStyle.h"
#include <vector>
#include <array>


compare::compare(){
  c1 = new TCanvas();
}
compare::~compare(){}

void compare::compare_Ene(int method){

  //if(method != 1) return;
  int Ene[6] = {10 , 30 ,50 ,80 ,100 ,150};
  for(int i = 0 ; i < 12 ; ++i){
    if(Ene[i%6] == 10 || Ene[i%6] == 30 || Ene[i%6] == 150) continue;
    // if(i == 6){
    //   sprintf(title,"Linear_result/HGLG_sat_all.root");
    //   store_GR(string(title),i);}
    if(i < 6){
      sprintf(title,"Linear_result/HGLG_sat_%iGeV.root",Ene[i]);
      store_GR(string(title),i);}
    else{
      sprintf(title,"Spline_result/HGLG_sat_Spline_%iGeV.root",Ene[i%6]);
      store_GR(string(title),i);}
  }


  bool savepng = true;

  for(int BD = 0; BD < MAXBD ; ++BD){
    TMultiGraph *mgr = new TMultiGraph();
    TLegend *leg = new TLegend(0.65,0.65,0.87,0.87);
    leg->SetBorderSize(0);
    for(int lab = 0 ; lab < MAXLABEL ; ++lab){
      if(gr_p0[lab][BD] == NULL) continue;
      if(Ene[lab%6] == 10 || Ene[lab%6] == 30 || Ene[lab%6] == 150) continue;
      if(lab < 6)
	gr_p0[lab][BD]->SetMarkerStyle(20);
      else
	gr_p0[lab][BD]->SetMarkerStyle(22);
      gr_p0[lab][BD]->SetMarkerSize(0.8);
      gr_p0[lab][BD]->SetMarkerColor((lab+1)%6);
      mgr->Add(gr_p0[lab][BD],"P");
      if( lab < 6 )
	sprintf(title,"Linear_%iGeV",Ene[lab%6]);
      else
	sprintf(title,"Spline_%iGeV",Ene[lab%6]);
      leg->AddEntry(gr_p0[lab][BD],title,"P");
    }
    
    string titlex("skiroc*32+channel/2");
    string titley("p0");
    Set_mgr(*mgr,titlex,titley);

    mgr->SetMaximum(50);
    mgr->SetMinimum(-50);

    mgr->Draw("sameAP");
    sprintf(title,"Board%i",BD);
    leg->SetHeader(title,"C");
    leg->Draw("same");
    c1->Update();
    sprintf(title,"plot_out/EnergyComp/BD%i_p0.png",BD);
    if(savepng)
      c1->SaveAs(title);
  }

  for(int BD = 0; BD < MAXBD ; ++BD){
    TMultiGraph *mgr = new TMultiGraph();
    TLegend *leg = new TLegend(0.65,0.65,0.87,0.87);
    leg->SetBorderSize(0);
    for(int lab = 0 ; lab < MAXLABEL ; ++lab){
      if(gr_p1[lab][BD] == NULL) continue;
      if(Ene[lab%6] == 10 || Ene[lab%6] == 30 || Ene[lab%6] == 150) continue;
      if(lab < 6)
	gr_p1[lab][BD]->SetMarkerStyle(20);
      else
	gr_p1[lab][BD]->SetMarkerStyle(22);
      gr_p1[lab][BD]->SetMarkerSize(0.8);
      gr_p1[lab][BD]->SetMarkerColor((lab+1)%6);
      mgr->Add(gr_p1[lab][BD],"P");
      if( lab < 6 )
	sprintf(title,"Linear_%iGeV",Ene[lab%6]);
      else
	sprintf(title,"Spline_%iGeV",Ene[lab%6]);
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
    sprintf(title,"plot_out/EnergyComp/BD%i_p1.png",BD);
    if(savepng)
      c1->SaveAs(title);
  }

  for(int BD = 0; BD < MAXBD ; ++BD){
    TMultiGraph *mgr = new TMultiGraph();
    TLegend *leg = new TLegend(0.65,0.65,0.87,0.87);
    leg->SetBorderSize(0);
    for(int lab = 0 ; lab < MAXLABEL ; ++lab){
      if(gr_sat[lab][BD] == NULL) continue;
      if(Ene[lab%6] == 10 || Ene[lab%6] == 30 || Ene[lab%6] == 150) continue;
      if(lab < 6)
	gr_sat[lab][BD]->SetMarkerStyle(20);
      else
	gr_sat[lab][BD]->SetMarkerStyle(22);
      gr_sat[lab][BD]->SetMarkerSize(0.8);
      gr_sat[lab][BD]->SetMarkerColor((lab+1)%6);
      mgr->Add(gr_sat[lab][BD],"P");
      if( lab < 6 )
	sprintf(title,"Linear_%iGeV",Ene[lab%6]);
      else
	sprintf(title,"Spline_%iGeV",Ene[lab%6]);
      leg->AddEntry(gr_sat[lab][BD],title,"P");
    }
    
    string titlex("skiroc*32+channel/2");
    string titley("sat");
    Set_mgr(*mgr,titlex,titley);

    mgr->SetMaximum(2500);
    mgr->SetMinimum(1500);

    mgr->Draw("sameAP");
    sprintf(title,"Board%i",BD);
    leg->SetHeader(title,"C");
    leg->Draw("same");
    c1->Update();
    sprintf(title,"plot_out/EnergyComp/BD%i_sat.png",BD);
    if(savepng)
      c1->SaveAs(title);
  }
  

  /*  
  for(int BD = 0; BD < MAXBD ; ++BD){
    TMultiGraph *mgr = new TMultiGraph();
    TLegend *leg = new TLegend(0.65,0.65,0.87,0.87);
    leg->SetBorderSize(0);
    for(int lab = 0 ; lab < MAXLABEL ; ++lab){
      if(gr_p0[lab][BD] == NULL) continue;
      gr_p0[lab][BD]->SetMarkerStyle(20);
      gr_p0[lab][BD]->SetMarkerSize(0.8);
      gr_p0[lab][BD]->SetMarkerColor(lab+1);
      mgr->Add(gr_p0[lab][BD],"P");
      if( lab == 6 )
	sprintf(title,"Energy_all");
      else
	sprintf(title,"Energy_%iGeV",Ene[lab%6]);
      leg->AddEntry(gr_p0[lab][BD],title,"P");
    }
    
    string titlex("skiroc*32+channel/2");
    string titley("p0");
    Set_mgr(*mgr,titlex,titley);

    mgr->SetMaximum(50);
    mgr->SetMinimum(-50);

    mgr->Draw("sameAP");
    sprintf(title,"Board%i",BD);
    leg->SetHeader(title,"C");
    leg->Draw("same");
    c1->Update();
    sprintf(title,"plot_out/EnergyComp/BD%i_p0.png",BD);
    if(savepng)
      c1->SaveAs(title);
  }

  
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
      if( lab == 6)
	sprintf(title,"Energy_all");
      else
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
    sprintf(title,"plot_out/EnergyComp/BD%i_p1.png",BD);
    if(savepng)
      c1->SaveAs(title);
  }


    for(int BD = 0; BD < MAXBD ; ++BD){
    TMultiGraph *mgr = new TMultiGraph();
    TLegend *leg = new TLegend(0.65,0.65,0.87,0.87);
    leg->SetBorderSize(0);
    for(int lab = 0 ; lab < MAXLABEL ; ++lab){
      if(gr_sat[lab][BD] == NULL) continue;
      gr_sat[lab][BD]->SetMarkerStyle(20);
      gr_sat[lab][BD]->SetMarkerSize(0.8);
      gr_sat[lab][BD]->SetMarkerColor(lab+1);
      mgr->Add(gr_sat[lab][BD],"P");
      if( lab == 6)
	sprintf(title,"Energy_all");
      else
	sprintf(title,"Energy_%iGeV",Ene[lab%6]);
      leg->AddEntry(gr_sat[lab][BD],title,"P");
    }
    
    string titlex("skiroc*32+channel/2");
    string titley("sat");
    Set_mgr(*mgr,titlex,titley);

    mgr->SetMaximum(2500);
    mgr->SetMinimum(1000);
    mgr->Draw("sameAP");
    sprintf(title,"Board%i",BD);
    leg->SetHeader(title,"C");
    leg->Draw("same");
    c1->Update();
    sprintf(title,"plot_out/EnergyComp/BD%i_sat.png",BD);
    if(savepng)
      c1->SaveAs(title);
  }
  */
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

  std::array< vector<double>, MAXBD> p0_vec,p1_vec,sat_vec,X_vec;
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
void compare::compare_method(int Ene){
  
  fitter::root_logon();
  sprintf(title,"HGLG_sat_%iGeV.root",Ene);
  TFile f_linear(title);
  sprintf(title,"HGLG_sat_Spline_%iGeV.root",Ene);
  TFile f_spline(title);

  TTree *T_linear = (TTree*)f_linear.Get("tree");
  TTree *T_spline = (TTree*)f_spline.Get("tree");

  TH1D *h[2][3];
  for(int i = 0 ; i < 2 ; ++i){
      h[i][0] = new TH1D("","",100,-50,50);
      h[i][1] = new TH1D("","",40,6,10);
      h[i][2] = new TH1D("","",120,1200,2400);
  }
  
  double p0,p1,sat;
  bool   goodsat;
  T_linear->SetBranchAddress("p0",&p0);
  T_linear->SetBranchAddress("p1",&p1);
  T_linear->SetBranchAddress("p_saturation",&sat);
  T_linear->SetBranchAddress("good_saturation",&goodsat);

  T_spline->SetBranchAddress("p0",&p0);
  T_spline->SetBranchAddress("p1",&p1);
  T_spline->SetBranchAddress("p_saturation",&sat);
  T_spline->SetBranchAddress("good_saturation",&goodsat);


  for(int i = 0 ;i < T_linear->GetEntries();++i ){
    T_linear->GetEntry(i);
    if(goodsat){
      h[0][0]->Fill(p0);
      h[0][1]->Fill(p1);
      h[0][2]->Fill(sat);
    }
  }
  
  for(int i = 0 ;i < T_spline->GetEntries();++i ){
    T_spline->GetEntry(i);
    if(goodsat){
      h[1][0]->Fill(p0);
      h[1][1]->Fill(p1);
      h[1][2]->Fill(sat);
    }
  }
  
  TLegend *leg = new TLegend(0.6,0.65,0.87,0.87);
  leg->SetBorderSize(0);
  gStyle->SetOptStat(0);
  for(int i = 0 ; i < 3 ; ++i){
    
    if( i == 0){
      leg->SetHeader("P0","C");
      leg->AddEntry(h[0][i],"Linear","L");
      leg->AddEntry(h[1][i],"Spline","L");
      h[0][i]->SetXTitle("P0");}
    
    else if( i == 1){
      leg->SetHeader("P1","C");
      h[0][i]->SetXTitle("P1");}
    else if( i == 2){
      leg->SetHeader("saturation_point","C");
      h[0][i]->SetXTitle("saturationP(HGADC)");}

    h[1][i]->SetLineColor(2);

    h[0][i]->SetLineWidth(2.5);
    h[1][i]->SetLineWidth(2.5);
    
    h[0][i]->SetMaximum( h[0][i]->GetMaximum()*1.2 );
    h[0][i]->Draw();
    h[1][i]->Draw("same");
    leg->Draw("same");
    c1->Update();
    if(i == 0) c1->SaveAs("compare_p0.png");
    if(i == 1) c1->SaveAs("compare_p1.png");
    if(i == 2) c1->SaveAs("compare_sat.png");

    //c1->WaitPrimitive();
  }


  
  // Correlation

  vector<double> linear_p0,linear_p1,linear_sat;
  vector<double> spline_p0,spline_p1,spline_sat;

  for(int i = 0 ;i < T_linear->GetEntries();++i ){
    T_linear->GetEntry(i);
    if(goodsat){
      h[0][0]->Fill(p0);
      h[0][1]->Fill(p1);
      h[0][2]->Fill(sat);
      linear_p0.push_back(p0);
      linear_p1.push_back(p1);
      linear_sat.push_back(sat);
      
      T_spline->GetEntry(i);
      if(goodsat){
	spline_p0.push_back(p0);
	spline_p1.push_back(p1);
	spline_sat.push_back(sat);
      }
      else{
	linear_p0.pop_back();
	linear_p1.pop_back();
	linear_sat.pop_back();	
      }
    }
  }

  TGraph *gr[3];
  for(int i = 0 ; i < 3 ; ++i){
    if( i == 0 ){
      gr[i] = new TGraph(linear_p0.size(),&linear_p0[0],&spline_p0[0]);
      gr[i]->GetXaxis()->SetTitle("linear_p0");
      gr[i]->GetYaxis()->SetTitle("spline_p0");}
    if( i == 1 ){
      gr[i] = new TGraph(linear_p0.size(),&linear_p1[0],&spline_p1[0]);
      gr[i]->GetXaxis()->SetTitle("linear_p1");
      gr[i]->GetYaxis()->SetTitle("spline_p1");}
    if( i == 2 ){
      gr[i] = new TGraph(linear_p0.size(),&linear_sat[0],&spline_sat[0]);
      gr[i]->GetXaxis()->SetTitle("linear_sat");
      gr[i]->GetYaxis()->SetTitle("spline_sat");}
    gr[i]->SetMarkerSize(1.2);
    gr[i]->SetMarkerStyle(20);
    gr[i]->SetMarkerColor(2);
    gr[i]->Draw("AP");
    gr[i]->SetTitle("");
    c1->Update();
    

    if(i == 0){
      c1->SaveAs("correlation_p0.png");}
    if(i == 1){
      c1->SaveAs("correlation_p1.png");}
    if(i == 2){
      gr[i]->GetYaxis()->SetTitleOffset(1.2);
      c1->SaveAs("correlation_sat.png");}


  }
  
  
}

void compare::compare_Data_Inj(){
  int MAXBD = 28;
  int MAXCHIP = 4;
  int MAXCH = 32;
  double HLcoeff[2][MAXBD][MAXCHIP][MAXCH];
  double LTcoeff[2][MAXBD][MAXCHIP][MAXCH];

}
