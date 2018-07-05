#include "fitter.h"
#include "TF1.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TProfile.h"

fitter::fitter(TGraph *in_gr,vector<double> HG,vector<double> LG,vector<double> TOT){
  gr = in_gr;
  HG_vec = HG;
  LG_vec = LG;
  TOT_vec = TOT;
  undershoot_percent = -1;
}

fitter::~fitter(){
}

void fitter::fit_Graph(){
  int npoint = gr->GetN();
  //cout << "Fitting size "<< npoint << endl;
  if(npoint < 1000){
    status = false;
    p0 = -1;
    p1 = -1;
    sat_point = -1;
    return;  }
  status = true;
  
  //TF1 *fit = new TF1("fit","pol1");
  //gr->Fit(fit,"0EMR");
  //p0 = fit->GetParameter(0);
  //p1 = fit->GetParameter(1);

  //fit_Draw();
}


void fitter::fit_Draw(){
  
  TCanvas *c1 = new TCanvas();
  TF1 *fit = new TF1("fit","pol1",0,500);
  TMultiGraph *mgr = new TMultiGraph();

  gr->Fit(fit,"0EMR");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.2);
  
  /*
  TGraph *g_res = new TGraph(fit_res.size(),&LG_vec[0],&fit_res[0]);
  g_res->Draw("AP");
  g_res->SetMaximum(0.5);
  g_res->SetMinimum(-0.5);
  c1->Update();
  cout << "res!" << endl;
  getchar();
  */
  int loop_time = 0;
  vector<double> HG_tmp,LG_tmp,HG_fit,LG_fit,fit_res,undershoot_tmp,othershoot_tmp,US_LG,OS_LG;
  TMultiGraph *mgr_fit = new TMultiGraph();
  mgr_fit->Add(gr,"P");
  fit_remove.resize(HG_vec.size());
  undershoot_tmp.resize(HG_vec.size());
  othershoot_tmp.resize(HG_vec.size());
  US_LG.resize(HG_vec.size());
  OS_LG.resize(HG_vec.size());
  
  vector<TF1 *> fit_S_vec;
  int LG_under_shoot = 0;
  while(true){
    TF1 *fit_S = new TF1("fit_S","[0]*x",0,500);
    fit_S_vec.push_back(fit_S);
    fit_res.clear();
    HG_tmp.clear();
    LG_tmp.clear();
    HG_fit.clear();
    LG_fit.clear();
    US_LG.clear();
    OS_LG.clear();
    
    for(int i = 0 ; i < (int)HG_vec.size() ; ++i){
      if(fit_remove.at(i)){
	HG_tmp.push_back(HG_vec.at(i));
	LG_tmp.push_back(LG_vec.at(i));}
      else{
	//cout << HG_vec.at(i) << " , " << LG_vec.at(i) << endl;
	HG_tmp.push_back(0);
	LG_tmp.push_back(0);
	HG_fit.push_back(HG_vec.at(i));
	LG_fit.push_back(LG_vec.at(i));}}

    
    for(int i = 0 ; i < (int)HG_vec.size() ; ++i){
      double res = HG_vec.at(i) -  fit->Eval(LG_vec.at(i));
      fit_res.push_back(res/HG_vec.at(i));
      bool remove = true ? abs(res/HG_vec.at(i)) > 0.2 : false ;
      if(loop_time==0 && LG_vec.at(i) < 500 && LG_vec.at(i) > 20)
	remove = false;
      if(HG_vec.at(i) > 2500) remove = true;
      fit_remove[i] = remove;

      if( loop_time==1 ){
	US_LG[i] = 0;
	OS_LG[i] = 0;
	undershoot_tmp[i] = 0;
	othershoot_tmp[i] = 0; 
	if( res/HG_vec.at(i) > 0.2 ) {
	  LG_under_shoot++;
	  undershoot_tmp[i] = HG_vec.at(i);
	  US_LG[i] = LG_vec.at(i);      }
	if( res/HG_vec.at(i) < 0.2 ) {
	  othershoot_tmp[i] = HG_vec.at(i);
	  OS_LG[i] = LG_vec.at(i);}}
      
    }

    TGraph *g_rm = new TGraph(HG_vec.size(),&LG_tmp[0],&HG_tmp[0]);
    TGraph *g_fit= new TGraph(HG_fit.size(),&LG_fit[0],&HG_fit[0]);
    TGraph *g_US = new TGraph(HG_fit.size(),&US_LG[0],&undershoot_tmp[0]);
    TGraph *g_OS = new TGraph(HG_fit.size(),&OS_LG[0],&othershoot_tmp[0]);
    
    if(loop_time == 1){
    g_US->SetMarkerColor(3);
    g_US->SetMarkerStyle(20);
    g_US->SetMarkerSize(0.2);
    mgr_fit ->Add(g_US,"P");
    g_OS->SetMarkerColor(2);
    g_OS->SetMarkerStyle(20);
    g_OS->SetMarkerSize(0.2);
    mgr_fit ->Add(g_OS,"P"); }
      
    g_fit->Fit(fit_S,"0EMR");
    fit_S->SetLineColor(loop_time+2);
    fit_S->SetLineWidth(2);
    
    g_rm->SetMarkerColor(4);
    g_rm->SetMarkerStyle(20);
    g_rm->SetMarkerSize(0.2);
      
    //mgr_fit -> Add(g_rm,"P");
    if(loop_time == 0)
      mgr_fit -> Draw("AP");
    else{
      mgr_fit -> Draw("APsame");
      for(size_t i = 0 ; i < fit_S_vec.size() ; ++i)
	fit_S_vec.at(i)->Draw("same");
    }
    mgr_fit -> SetMinimum(0);
    mgr_fit -> SetMaximum(4096);
    c1->Update();
    //getchar();
    
    double new_p0 = 0;
    double new_p1 = fit_S->GetParameter(0);
    if(LG_under_shoot*100./HG_fit.size() > 5) {
      mgr_fit -> Draw("same");
      c1->Update();
      getchar();}
      
    if(abs(p1 - new_p1) < 0.01 && loop_time > 1) {
      undershoot_percent = LG_under_shoot*100./HG_fit.size();
      p0 = new_p0;
      p1 = new_p1;
      cout << "happy break! loop time = " << loop_time << endl;
      break;}
    p0 = new_p0;
    p1 = new_p1;
    loop_time++;
  }

  /*
  mgr -> Add(gr);
  mgr -> Add(g_rm);
  mgr -> Draw("AP");
  mgr -> SetMinimum(0);
  mgr -> SetMaximum(4096);

  //mgr->Fit(fit,"EMR");
  fit->Draw("same");
  c1->Update();
  getchar();
  //fit->Draw("same");
  c1->Update();
  getchar();
  */
  delete c1;

}
