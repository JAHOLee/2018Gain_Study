#include "fitter.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include <algorithm>

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
  if(!(labelE == 10 || labelE == 30 || labelE == 50 || labelE == 80
       || labelE == 100 || labelE == 150 )) {
    cout << "invalid energy!" << endl;
    return;}
  char title[50];
  sprintf(title,"root_result/400Bin/%iGeV.root",labelE);
  TFile f(title);

   // TF1 *sat_fit = new TF1("1st_try"," [1]* (TMath::Exp(x/[0]) / (TMath::Exp(x/[0]) + 1) -0.5 )",0,800);
   // sat_fit->SetParLimits(0,50,120);
   // sat_fit->SetParLimits(1,2000,6000);

  //TF1 *sat_fit = new TF1("2nd_try"," [0]*tanh(x*[1])",0,500);
  //TF1 *sat_fit = new TF1("3rd_try"," [0]*([2]*x/(1+abs(x))+tanh(x*[1]))",0,500);
  TF1 *sat_fit   = new TF1("4th_try","[0]*x+[1]",50,150);
  sat_fit->SetParLimits(1,-50,50);
  double p0_ARR[MAXBD][MAXSKI][MAXCH];
  double p1_ARR[MAXBD][MAXSKI][MAXCH];
  double sat_ARR[MAXBD][MAXSKI][MAXCH];
  bool   sat_good[MAXBD][MAXSKI][MAXCH];

  gStyle->SetOptStat(0);
  TProfile *tpr = (TProfile *)f.Get("Board_7/HGLG_chip2_ch44");
  //  TProfile *tpr = new TProfile("","",200,0,800,0,4000);
  TH1D *h1 = new TH1D("","",tpr->GetNbinsX(),0,800);
  int rebinN = 2;
  h1->Rebin(rebinN);
  bool savepng = 0;
  int stop_and_look = -1;
  gROOT->SetBatch(kTRUE);

  for(int BD = 0 ;BD < MAXBD ; ++BD){
    cout << "BD "<< BD << endl;
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      cout << "SKI "<< SKI << endl;
      //if( BD == 8 && (SKI ==2 || SKI ==3))
      //stop_and_look = 1;
      //else
      //stop_and_look = -1;
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	p0_ARR[BD][SKI][CH] = -999;
	p1_ARR[BD][SKI][CH] = -999;
	sat_ARR[BD][SKI][CH] = -999;
	sat_fit->SetRange(50,150);
        
	int true_ch = CH*2;
	sprintf(title,"Board_%i/HGLG_chip%i_ch%i",BD,SKI,true_ch);
	tpr = (TProfile *)f.Get(title);
	if(tpr == NULL) continue;

	tpr->Rebin(rebinN);
	
	int nentry = tpr->GetEntries();
	if( nentry < 1000 ) continue;
	tpr->SetName(title);
	//tpr->Draw();
	//c1->Update();
	//getchar();
 
	tpr->Fit(sat_fit,"QEMR");
	//tpr->Fit(sat_fit2,"QEMR");
	sat_fit->Draw("same");
	sat_fit->SetLineColor(6);
	
	h1->Reset();
	//h1->Sumw2();
	double fit_start = 0,fit_end = 0;
	bool right_most = false,left_most = false;
	for(int i = 0 ; i < tpr->GetNbinsX () ; ++i){
	  double x = tpr->GetBinCenter(i);
	  double y = tpr->GetBinContent(i);
	  if(x == 0 || y == 0) continue;
	  //cout << "x = "<< x << ", y = " << y << endl;
	  double res = ( y - sat_fit->Eval(x)) / sat_fit->Eval(x);
	  double error = tpr->GetBinError(i)/sat_fit->Eval(x);
	  h1->SetBinContent(i,res);
	  h1->SetBinError(i,error);
	  if( abs(res) < 0.03 && x > 200 && !right_most) {
	    right_most = true;
	    fit_end = x;
	  }
	  if( abs(res) < 0.03 && x < 70 && !left_most) {
	    left_most = true;
	    fit_start = x;
	  }

	}

	tpr->SetMarkerColor(SKI+1);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(1.2);
	h1->SetMarkerColor(SKI+1);
	ratio_plot(tpr,sat_fit,h1);
	if(stop_and_look == 0)
	  getchar();
	if(savepng)
	  c1->SaveAs("step1.png");
	
	sat_fit->SetRange(fit_start,fit_end);
        tpr->Fit(sat_fit,"QEMR");
	sat_fit->Draw("same");
	h1->Reset();
	double sat_point = 0;
	double sat_point_x = 0;
	
	for(int i = 0 ; i < tpr->GetNbinsX () ; ++i){
	  double x = tpr->GetBinCenter(i);
	  double y = tpr->GetBinContent(i);
	  if(x == 0 || y == 0) continue;
	  //cout << "x = "<< x << ", y = " << y << endl;
	  double res = ( y - sat_fit->Eval(x)) / sat_fit->Eval(x);
	  double error = tpr->GetBinError(i)/(x*sat_fit->Eval(x));
	  h1->SetBinContent(i,res);
	  h1->SetBinError(i,error);	}
	
	ratio_plot(tpr,sat_fit,h1);
	if(stop_and_look == 0)
	  getchar();
	
	if(savepng)
	  c1->SaveAs("step2.png");

	double lowx,lowy,highx,highy;
	Find_low(h1,&lowx,&lowy);
	Find_high(h1,&highx,&highy);
	//cout << lowx << " , " << sat_fit->Eval(lowx) << endl;
	//cout << highx << " , " << sat_fit->Eval(highx) << endl;
	sat_fit->SetRange(lowx,fit_end);
	tpr->Fit(sat_fit,"QEMR");
	ratio_plot(tpr,sat_fit,h1);
	sat_point_x = highx;
	sat_point   = sat_fit->Eval(highx);

	double thres = Calc_avg(h1,lowx,highx);
	
	if(stop_and_look == 0)
	  getchar();

	sat_good[BD][SKI][CH] = Find_sat(tpr,h1,&sat_point,&sat_point_x,thres);
	
	h1->Reset();
	for(int i = 0 ; i < tpr->GetNbinsX () ; ++i){
	  double x = tpr->GetBinCenter(i);
	  double y = tpr->GetBinContent(i);
	  if(x == 0 || y == 0) continue;
	  //cout << "x = "<< x << ", y = " << y << endl;
	  double res = ( y - sat_fit->Eval(x)) / sat_fit->Eval(x);
	  double error = tpr->GetBinError(i)/(x*sat_fit->Eval(x));
	  h1->SetBinContent(i,res);
	  h1->SetBinError(i,error);	 
	}

	//sat_good[BD][SKI][CH] = Find_sat(tpr,h1,&sat_point,&sat_point_x);
	/*
	if( sat_fit->Eval(highx) < sat_point || sat_point < 500 ){	  
	  sat_point_x = highx;
	  sat_point   = sat_fit->Eval(highx);	}
	*/

	//ratio_plot(tpr,sat_fit,h1);
	//cout << "coooooooooooool" << endl;
	//getchar();
	
	// TF1 *sat_fit2   = new TF1("4th_try2","[0]*x",0,fit_end);
	// tpr->Fit(sat_fit2);

	// for(int i = 0 ; i < tpr->GetNbinsX () ; ++i){
	//   double x = tpr->GetBinCenter(i);
	//   double y = tpr->GetBinContent(i);
	//   if(x == 0 || y == 0) continue;
	//   //cout << "x = "<< x << ", y = " << y << endl;
	//   double res = ( y - x*p0[BD][SKI][CH]) / (x*p0[BD][SKI][CH]);
	//   double error = tpr->GetBinError(i)/(x*p0[BD][SKI][CH]);
	//   h1->SetBinContent(i,res);
	//   h1->SetBinError(i,error);
	// }
		
	

	//h1->Draw("e");
	//sat_fit_2->Draw("same");
	//sat_fit_2->SetLineColor(7);
	//tpr->Draw();
	// sat_fit2->SetLineColor(7);
	// sat_fit2->SetMarkerStyle(20);
	// sat_fit2->SetMarkerSize(1.2);
	// sat_fit2->Draw("same");
	ratio_plot(tpr,sat_fit,h1);
	TLine *Gline = new TLine(0,sat_point,tpr->GetXaxis()->GetXmax(),sat_point);
	Gline->SetLineColor(1);
	Gline->SetLineWidth(4.8);
	//Gline->SetLineStyle(7);
	pad1->cd();
	Gline->Draw("same");  
	c1->cd();
	c1->Update();
	
	if(stop_and_look == 0)
	  getchar();
	if(savepng)
	  c1->SaveAs("step3.png");
	
	sat_fit->SetRange(lowx,sat_point_x);
	tpr->Fit(sat_fit,"QEMR");
	ratio_plot(tpr,sat_fit,h1);
	
	pad1->cd();
	Gline->Draw("same");  
	c1->cd();
	c1->Update();
	
	if(stop_and_look == 1)
	  getchar();
	
	if(savepng)
	  c1->SaveAs("step4.png");
	
	if(savepng)
	  getchar();
	
	//if( sat_point > 2000 || sat_point < 500) getchar();
	if( sat_point < 500 || sat_fit->GetParameter(0) < 6 ||sat_point > 2400){
	  p0_ARR[BD][SKI][CH] = -999;
	  p1_ARR[BD][SKI][CH] = -999;
	  sat_ARR[BD][SKI][CH] = -999;
	  sat_good[BD][SKI][CH] = false;	}
	
	else{
	  p0_ARR[BD][SKI][CH] = sat_fit->GetParameter(0);
	  p1_ARR[BD][SKI][CH] = sat_fit->GetParameter(1);
	  sat_ARR[BD][SKI][CH] = sat_point;}
      }
    }
  }

  TMultiGraph *mgr[3];
  TGraph *gr;
  TGraph *gr_sat;
  TLegend *leg;
  
  for(int BD = 0 ;BD < MAXBD ; ++BD){
    mgr[0] = new TMultiGraph();
    mgr[1] = new TMultiGraph();
    mgr[2] = new TMultiGraph();
    leg = new TLegend(0.7,0.7,0.87,0.87);
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      vector<double> ch_arr;
      vector<double> p0_arr;
      vector<double> p1_arr;
      vector<double> sat_arr;
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	if(p0_ARR[BD][SKI][CH] != -999){
	  ch_arr.push_back(CH*2);
	  p0_arr.push_back(p0_ARR[BD][SKI][CH]);
	  p1_arr.push_back(p1_ARR[BD][SKI][CH]);
	  sat_arr.push_back(sat_ARR[BD][SKI][CH]);}	  
      }
      
      gr = new TGraph(ch_arr.size(),&ch_arr[0],&p0_arr[0]);
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(1.2);
      gr->SetMarkerColor(SKI+1);
      mgr[0]->Add(gr);

      gr_sat = new TGraph(ch_arr.size(),&ch_arr[0],&sat_arr[0]);
      gr_sat->SetMarkerStyle(20);
      gr_sat->SetMarkerSize(1.2);
      gr_sat->SetMarkerColor(SKI+1);
      mgr[1]->Add(gr_sat);

      gr_sat = new TGraph(ch_arr.size(),&ch_arr[0],&p1_arr[0]);
      gr_sat->SetMarkerStyle(20);
      gr_sat->SetMarkerSize(1.2);
      gr_sat->SetMarkerColor(SKI+1);
      mgr[2]->Add(gr_sat);

      sprintf(title,"chip %i",SKI);
      leg->AddEntry(gr,title,"P");
      
    }
    mgr[0]->Draw("AP");
    mgr[0]->GetXaxis()->SetTitle("CH");
    mgr[0]->GetYaxis()->SetTitle("p1");
    mgr[0]->SetMaximum(13);
    leg->Draw("same");
    
    c1->Update();
    sprintf(title,"plot_out/%iGeV/BD%i_p1.png",labelE,BD);
    c1->SaveAs(title);

    mgr[2]->Draw("AP");
    mgr[2]->GetXaxis()->SetTitle("CH");
    mgr[2]->GetYaxis()->SetTitle("p0");
    mgr[2]->SetMaximum(40);
    mgr[2]->SetMinimum(-40);
    leg->Draw("same");
    
    c1->Update();
    sprintf(title,"plot_out/%iGeV/BD%i_p0.png",labelE,BD);
    c1->SaveAs(title);

    
    mgr[1]->Draw("AP");
    mgr[1]->GetXaxis()->SetTitle("CH");
    mgr[1]->GetXaxis()->SetTitle("sat_point");
    mgr[1]->SetMaximum(3000);
      
    leg->Draw("same");
    c1->Update();
    sprintf(title,"plot_out/%iGeV/BD%i_sat.png",labelE,BD);
    c1->SaveAs(title);
    
    //getchar();
  }


  //Modify the wrong ones
  double p0_avg[MAXBD][MAXSKI];
  double p1_avg[MAXBD][MAXSKI];
  double sat_avg[MAXBD][MAXSKI];
  int    count_avg[MAXBD][MAXSKI];
  
  for(int BD = 0 ;BD < MAXBD ; ++BD){
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      p0_avg[BD][SKI] = 0;
      p1_avg[BD][SKI] = 0;
      sat_avg[BD][SKI] = 0;
      count_avg[BD][SKI] = 0;    }}

  for(int BD = 0 ;BD < MAXBD ; ++BD){
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	if(sat_ARR[BD][SKI][CH] == -999) continue;
	p0_avg[BD][SKI]    += p0_ARR[BD][SKI][CH];
	p1_avg[BD][SKI]    += p1_ARR[BD][SKI][CH];
	sat_avg[BD][SKI]   += sat_ARR[BD][SKI][CH];
	count_avg[BD][SKI] += 1;
      }
      p0_avg[BD][SKI]  /= count_avg[BD][SKI];
      p1_avg[BD][SKI]  /= count_avg[BD][SKI];
      sat_avg[BD][SKI] /= count_avg[BD][SKI];
    }
  }
  sprintf(title,"HGLG_sat_%iGeV.root",labelE);
  TFile outf(title,"recreate");
  TTree *outtree = new TTree("tree","tree");
  int layerID,skirocID,channelID;
  double m_p0,m_p1,m_sat;
  bool   m_goodsat;

  outtree->Branch("layerID",&layerID,"layerID/I");
  outtree->Branch("skirocID",&skirocID,"skirocID/I");
  outtree->Branch("channelID",&channelID,"channelID/I");
  outtree->Branch("p0",&m_p0,"p0/D");
  outtree->Branch("p1",&m_p1,"p1/D");
  outtree->Branch("p_saturation",&m_sat,"p_saturation/D");
  outtree->Branch("good_saturation",&m_goodsat,"good_saturation/O");
  
  
  for(int BD = 0 ;BD < MAXBD ; ++BD){
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	if( sat_ARR[BD][SKI][CH] == -999 ){
	  layerID  = BD;
	  skirocID = SKI;
	  channelID = CH*2;
	  m_p0  = p1_avg[BD][SKI];
	  m_p1  = p0_avg[BD][SKI];
	  m_sat = sat_avg[BD][SKI];
	  m_goodsat = false;
	  outtree->Fill();	}
	else{
	  layerID  = BD;
	  skirocID = SKI;
	  channelID = CH*2;
	  m_p0  = p1_ARR[BD][SKI][CH];
	  m_p1  = p0_ARR[BD][SKI][CH];
	  m_sat = sat_ARR[BD][SKI][CH];
	  m_goodsat = sat_good[BD][SKI][CH];
	  outtree->Fill();}
      }
    }
  }

  TH1D* hp0 = qualify(p1_ARR,2); //Actually p0
  TH1D* hp1 = qualify(p0_ARR,0); //Actually p1
  TH1D* hsat = qualify(sat_ARR,1);
  hp0->SetTitle("p0_qualify");
  hp0->SetName("p0_qualify");
  hp1->SetTitle("p1_qualify");
  hp1->SetName("p1_qualify");
  hsat->SetTitle("sat_qualify");
  hsat->SetName("sat_qualify");

  
  sprintf(title,"plot_out/%iGeV/p0_hist.png",labelE);
  hp0 ->Draw(); c1->Update(); c1->SaveAs(title);
  sprintf(title,"plot_out/%iGeV/p1_hist.png",labelE);
  hp1 ->Draw(); c1->Update(); c1->SaveAs(title);
  sprintf(title,"plot_out/%iGeV/sat_hist.png",labelE);
  hsat->Draw(); c1->Update(); c1->SaveAs(title);

  outf.Write();
  outf.Close();
  
}
  
void fitter::Find_low(TH1D* h1,double* lowx,double* lowy){
  double lowerbound = 100;
  int Nbin = h1->GetNbinsX ();
  double minx = 1,miny = 1;
  for(int i = 0 ; i < Nbin ; ++i){
    double x = h1->GetBinCenter(i);
    double y = h1->GetBinContent(i);
    if(y < miny && x < lowerbound){
      minx = x;
      miny = y;    }
  }
  *lowx = minx;
  *lowy = miny;
}
void fitter::Find_high(TH1D* h1,double* highx,double* highy){
  double lowerbound = 200;
  double upperbound = 300;
  
  int Nbin = h1->GetNbinsX ();
  double maxx = -1,maxy = -1;
  for(int i = 0 ; i < Nbin ; ++i){
    double x = h1->GetBinCenter(i);
    double y = h1->GetBinContent(i);
    if(y > maxy && x > lowerbound && x < upperbound){
      maxx = x;
      maxy = y;    }
  }
  *highx = maxx;
  *highy = maxy;
}

bool fitter::Find_sat(TProfile *tpr, TH1D* h1 , double *sat,double *sat_x,double thres){
  int Nbin = tpr->GetNbinsX ();
  int keeplow;
  bool ret_first = false;
  
  for(int i = 0 ; i < Nbin; ++i){
    double x = h1->GetBinCenter(i);
    double y = h1->GetBinContent(i);
    keeplow = 0;
    for( int check = 0 ; check < 5 ; ++check){
      if( check+i == Nbin ) break;
      if( h1->GetBinContent(i+check) < thres )
	keeplow++;    }
    if( x > 200 && y < thres - 0.01 && keeplow >= 3 && !ret_first){
      *sat   = tpr->GetBinContent(i);
      *sat_x = tpr->GetBinCenter(i);
      ret_first = true;
      if( *sat != 0 ) return true;
    }
  }
  /*
  if( *sat == 0 )
  for(int i = 0 ; i < Nbin; ++i){
    double x = h1->GetBinCenter(i);
    double y = h1->GetBinContent(i);
    if( x > 150 && y < - 0.01 && !ret_first){
      *sat   = tpr->GetBinContent(i);
      *sat_x = tpr->GetBinCenter(i);
      ret_first = true;
      if( *sat != 0 ) return true;
    }
  }
  */
  *sat   = -999;
  *sat_x = -999;
  cout << "fail!" << endl;
  return false;
    
}

double fitter::Calc_avg(TH1D* h1 , double min,double max){
  int Nbin = h1->GetNbinsX ();
  int    member = 0;
  double sum    = 0;

  for(int i = 0 ; i < Nbin; ++i){
    double x = h1->GetBinCenter(i);
    double y = h1->GetBinContent(i);
    if(x >= min && x <= max){
      member++;
      sum += y;    }
  }
  sum /= member;
  return sum;
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

TH1D* fitter::qualify( double ARR[MAXBD][MAXSKI][MAXCH],int option){
  TH1D *h_qua;
  if(option == 0)
    h_qua = new TH1D("","",20,0,0.1);
  else if(option == 1)
    h_qua = new TH1D("","",80,0,0.4);
  else if(option == 2)
    h_qua = new TH1D("","",200,0,1);
  else{
    cout << "wrong input option, return NULL!" << endl;
    return NULL;  }
  cout << "==============================" << endl;
  for(int BD = 0 ;BD < MAXBD ; ++BD){
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      vector<double> data;
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	data.push_back(ARR[BD][SKI][CH]);}
      std::sort( data.begin(), data.end() );
      double sigma_first = 0.16;
      double sigma_last  = 0.84;
      double sigma = data[data.size()*sigma_last]-data[data.size()*sigma_first];
      double median= data[data.size()*0.5];
      if (median < 0) median = -median;
      double quality = sigma/median;
      h_qua->Fill(quality);
    }
  }
  return h_qua;
}

void fitter::fit_Draw(){
  

  TF1 *fit = new TF1("fit","pol1",0,500);

  gr->Fit(fit,"0EMR");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.2);
  
  /*
  TMultiGraph *mgr = new TMultiGraph();
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
  pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
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

  pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
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
  
  hratio->GetYaxis()->SetTitle("#frac{HG - Fit}{Fit}");
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
  //getchar();
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
