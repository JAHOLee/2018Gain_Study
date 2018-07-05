#include "fitter.h"
#include "TFile.h"
#include "TROOT.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TStyle.h"
fitter::fitter(TGraph *in_gr,vector<double> HG,vector<double> LG,vector<double> TOT){
  c1 = new TCanvas();
  gr = in_gr;
  HG_vec = HG;
  LG_vec = LG;
  TOT_vec = TOT;
  undershoot_percent = -1;
}
fitter::fitter(){
  c1 = new TCanvas();
}
fitter::~fitter(){
}
void fitter::fit(int labelE = 10){
  root_logon();
  if(!(labelE == 10 || labelE == 30 || labelE == 50 || labelE == 80
       || labelE == 100 || labelE == 150 )) {
    cout << "invalid energy!" << endl;
    return;}
  char title[50];
  sprintf(title,"root_result/%iGeV.root",labelE);
  TFile f(title);
  int MAXBD  = 28;
  int MAXSKI = 4;
  int MAXCH  = 32;

  // TF1 *sat_fit = new TF1("1st_try"," [1]* (TMath::Exp(x/[0]) / (TMath::Exp(x/[0]) + 1) -0.5 )",0,500);
  // sat_fit->SetParLimits(0,50,120);
  // sat_fit->SetParLimits(1,2000,6000);

  //TF1 *sat_fit = new TF1("2nd_try"," [0]*tanh(x*[1])",0,500);
  //TF1 *sat_fit = new TF1("3rd_try"," [0]*([2]*x/(1+abs(x))+tanh(x*[1]))",0,500);
  TF1 *sat_fit   = new TF1("4th_try","[0]*x",0,500);
  double p0[MAXBD][MAXSKI][MAXCH];
  //TF1 *sat_fit_2 = new TF1("4th_try_2","[0]*x+[1]",350,500);
  gStyle->SetOptStat(0);
  TProfile *tpr = new TProfile("","",400,0,800,0,4000);
  TH1D *h1 = new TH1D("","",tpr->GetNbinsX(),0,800);
  for(int BD = 0 ;BD < MAXBD ; ++BD){
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	p0[BD][SKI][CH] = -999;
	int true_ch = CH*2;
	sprintf(title,"Board_%i/HGLG_chip%i_ch%i",BD,SKI,true_ch);
	tpr = (TProfile *)f.Get(title);
	if(tpr == NULL) continue;
	int nentry = tpr->GetEntries();
	if( nentry < 1000 ) continue;
	tpr->SetName(title);
	//tpr->Draw();
	//c1->Update();
	//getchar();
	
	tpr->Fit(sat_fit,"EMR");
	//tpr->Fit(sat_fit_2,"EMR");
	sat_fit->Draw("same");
	sat_fit->SetLineColor(6);
	p0[BD][SKI][CH] = sat_fit->GetParameter(0);
	
	h1->Reset();
	//h1->Sumw2();
	for(int i = 0 ; i < tpr->GetNbinsX () ; ++i){
	  double x = tpr->GetBinCenter(i);
	  double y = tpr->GetBinContent(i);
	  if(x == 0 || y == 0) continue;
	  //cout << "x = "<< x << ", y = " << y << endl;
	  double res = ( y - x*p0[BD][SKI][CH]) / (x*p0[BD][SKI][CH]);
	  h1->SetBinContent(i,res);
	}
	tpr->SetMarkerColor(SKI+1);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(1.2);
	h1->SetMarkerColor(SKI+1);

	//h1->Draw("e");
	//sat_fit_2->Draw("same");
	//sat_fit_2->SetLineColor(7);
	//tpr->Draw();
	ratio_plot(tpr,sat_fit,h1);
	c1->Update();
	getchar();
      }
    }
  }

  TMultiGraph *mgr;
  TGraph *gr;
  //  for(int BD = 0 ;BD < MAXBD ; ++BD){
  int BD = 7;
    mgr = new TMultiGraph();
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      vector<double> ch_arr;
      vector<double> p0_arr;
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	if(p0[BD][SKI][CH] != -999){
	  ch_arr.push_back(CH*2);
	  p0_arr.push_back(p0[BD][SKI][CH]);
	}}
      gr = new TGraph(ch_arr.size(),&ch_arr[0],&p0_arr[0]);
      gr->SetMinimum(0);
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(1.2);
      gr->SetMarkerColor(SKI+1);
      mgr->Add(gr);
    }
    //}
  mgr->Draw("AP");
  c1->Update();
  getchar();
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


}

void fitter::ratio_plot(TProfile *tpr,TF1 *fit,TH1D *hratio){
  
  //Ratio plot
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.02); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  tpr->Draw();              // Draw h1
  fit->Draw("same");
  
  // Do not draw the Y axis label on the upper plot and redraw a small
  // axis instead, in order to avoid the first label (0) to be clipped.
  //h1->GetYaxis()->SetLabelSize(0.);
  
  TAxis *axis = tpr->GetYaxis();
  axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis->SetLabelSize(15);
  axis->SetTitleSize(0.06);
  axis->SetTitleOffset(0.7);
  axis->Draw();

  //h1->SetXTitle("DAC");
  //h1->GetXaxis()->SetTitleOffset(4);
  axis = tpr->GetXaxis();
  axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis->SetLabelSize(0);
  axis->Draw();

  // lower plot will be in pad
  c1->cd();          // Go back to the main canvas before defining pad2

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  pad2->SetTopMargin(0.1);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  hratio->SetTitle("");
  //hratio->SetXTitle("DAC");
  //hratio->SetYTitle("ratio");
  
  hratio->SetMaximum(0.1);
  hratio->SetMinimum(-0.1); 
  //hratio->Draw();
  
  hratio->GetYaxis()->SetTitle("#frac{Fit - HG}{Fit}");
  hratio->GetYaxis()->SetTitleOffset(0.35);
  hratio->GetYaxis()->SetLabelFont(43); 
  hratio->GetYaxis()->SetLabelSize(15);

  axis = hratio->GetYaxis();
  axis->SetNdivisions(404);
  axis->Draw();

  axis = hratio->GetXaxis();
  axis->SetTitle("LG");
  axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis->SetLabelSize(15);
  
  

  hratio->GetXaxis()->SetTitleSize(0.12);
  hratio->GetXaxis()->SetTitleOffset(0.76);
  hratio->GetYaxis()->SetTitleSize(0.12);
  //axis->Draw();

  
  hratio->Draw("e");
  //Add TLine to show mean value
  TLine *Gline = new TLine(0,0,tpr->GetXaxis()->GetXmax(),0);
  Gline->SetLineColor(1);
  Gline->SetLineWidth(4.8);
  Gline->SetLineStyle(7);
  Gline->Draw();  
  c1->Update();
  getchar();
}

void fitter::root_logon(){

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
