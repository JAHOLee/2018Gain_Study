#include "fitter.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include <algorithm>
#include "TSpline.h"

const int MINPOINT = 5;
const int MAXPOINT = 40;

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
void fitter::fit(int labelE ){
  if(!(labelE == 10 || labelE == 30 || labelE == 50 || labelE == 80
       || labelE == 100 || labelE == 150 )) {
    cout << "invalid energy!" << endl;
    return;}
  char title[50];
  sprintf(title,"root_result/400Bin/update/%iGeV.root",labelE);
  TFile f(title);

   // TF1 *sat_fit = new TF1("1st_try"," [1]* (TMath::Exp(x/[0]) / (TMath::Exp(x/[0]) + 1) -0.5 )",0,800);
   // sat_fit->SetParLimits(0,50,120);
   // sat_fit->SetParLimits(1,2000,6000);

  //TF1 *sat_fit = new TF1("2nd_try"," [0]*tanh(x*[1])",0,500);
  //TF1 *sat_fit = new TF1("3rd_try"," [0]*([2]*x/(1+abs(x))+tanh(x*[1]))",0,500);  
  TF1 *sat_fit   = new TF1("4th_try","[0]*x+[1]",MINPOINT,150);
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
  int stop_and_look = 0;
  //gROOT->SetBatch(kTRUE);

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
	sat_good[BD][SKI][CH] = false;
	sat_fit->SetRange(MINPOINT,MAXPOINT);
        
	int true_ch = CH*2;
	sprintf(title,"Board_%i/HGLG_chip%i_ch%i",BD,SKI,true_ch);
	tpr = (TProfile *)f.Get(title);
	if(tpr == NULL) continue;
	if(tpr->GetEntries() < 1500) continue;
	if(BD == 9 && SKI == 3 && CH*2 == 36 ) continue;
	if(BD == 18 && SKI == 2 && CH*2 == 12) continue;
	if(BD == 19 && SKI == 0 && CH*2 == 46) continue;

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
	  //c1->WaitPrimitive();
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
	  //c1->WaitPrimitive();
	
	if(savepng)
	  c1->SaveAs("step2.png");

	double lowx,lowy,highx,highy;
	Find_low(h1,&lowx,&lowy);
	Find_high(h1,&highx,&highy);
	//cout << lowx << " , " << sat_fit->Eval(lowx) << endl;
	//cout << highx << " , " << sat_fit->Eval(highx) << endl;
	sat_fit->SetRange(MINPOINT,fit_end);
	tpr->Fit(sat_fit,"QEMR");
	ratio_plot(tpr,sat_fit,h1);
	sat_point_x = highx;
	sat_point   = sat_fit->Eval(highx);

	double thres = Calc_avg(h1,lowx,highx);
	
	if(stop_and_look == 0)
	  //c1->WaitPrimitive();

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
	  //c1->WaitPrimitive();
	if(savepng)
	  c1->SaveAs("step3.png");
	
	sat_fit->SetRange(MINPOINT,sat_point_x);
	tpr->Fit(sat_fit,"QEMR");
	ratio_plot(tpr,sat_fit,h1);
	
	pad1->cd();
	Gline->Draw("same");  
	c1->cd();
	c1->Update();
	
	if(stop_and_look == 1)
	  //c1->WaitPrimitive();
	
	if(savepng)
	  c1->SaveAs("step4.png");
	
	if(savepng)
	  c1->WaitPrimitive();
	
	//if( sat_point > 2000 || sat_point < 500) getchar();
	if( sat_fit->GetParameter(0) < 6 || sat_fit->GetParameter(0) > 10 ||sat_point > 2400 || sat_point < 500 || sat_point_x < 150 || sat_point_x > 300 ){
	  p0_ARR[BD][SKI][CH] = -999;
	  p1_ARR[BD][SKI][CH] = -999;
	  sat_ARR[BD][SKI][CH] = -999;
	  sat_good[BD][SKI][CH] = false;	}
	
	else{
	  p0_ARR[BD][SKI][CH] = sat_fit->GetParameter(0);
	  p1_ARR[BD][SKI][CH] = sat_fit->GetParameter(1);
	  sat_ARR[BD][SKI][CH] = sat_point;
	  sat_good[BD][SKI][CH] = true;
	}
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

      sprintf(title,"chip %i",SKI);
      leg->AddEntry(gr,title,"P");
      
    }
    /*
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
    */    
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
	if( sat_ARR[BD][SKI][CH] == -999 || sat_good[BD][SKI][CH] == false){
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

void fitter::look_detail(){


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
      //c1->WaitPrimitive();
    }
      
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

void fitter::fit_spline(int labelE){
  if(!(labelE == 10 || labelE == 30 || labelE == 50 || labelE == 80
       || labelE == 100 || labelE == 150 )) {
    cout << "invalid energy!" << endl;
    return;}
  //gROOT->SetBatch(kTRUE);
  char title[50];
  sprintf(title,"root_result/400Bin/update/%iGeV.root",labelE);
  TFile f(title);

  double p0_ARR[MAXBD][MAXSKI][MAXCH];
  double p1_ARR[MAXBD][MAXSKI][MAXCH];
  double sat_ARR[MAXBD][MAXSKI][MAXCH];
  bool   sat_good[MAXBD][MAXSKI][MAXCH];
  int    tpr_entry[MAXBD][MAXSKI][MAXCH];
  

  gStyle->SetOptStat(0);
  TProfile *tpr = (TProfile *)f.Get("Board_7/HGLG_chip2_ch44");
  //TH1D *h1 = new TH1D("","",tpr->GetNbinsX(),0,800);
  int rebinN = 2;
  //h1->Rebin(rebinN);
  bool savepng = 0;
  //int stop_and_look = -1;
  //gROOT->SetBatch(kTRUE);
  
  TF1 *linear = new TF1("","[0]+x*[i]",MINPOINT,150);
  linear->SetParLimits(0,-50,50);
  TPad pad1_sp("pad1_sp","",0,0,1,1);
  TPad pad2_sp("pad2_sp","",0,0,1,1);
  
  for(int BD = 0 ;BD < MAXBD ; ++BD){
    cout << "BD "<< BD << endl;
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      cout << "SKI "<< SKI << endl;
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	p0_ARR[BD][SKI][CH] = -999;
	p1_ARR[BD][SKI][CH] = -999;
	sat_ARR[BD][SKI][CH] = -999;
	sat_good[BD][SKI][CH] = false;
	tpr_entry[BD][SKI][CH] = 0;
	int true_ch = CH*2;
	sprintf(title,"Board_%i/HGLG_chip%i_ch%i",BD,SKI,true_ch);
	tpr = (TProfile *)f.Get(title);
	if(tpr == NULL) continue;
	tpr_entry[BD][SKI][CH] = tpr->GetEntries();
	if(tpr->GetEntries() < 1500) continue;
	if(BD == 9 && SKI == 3 && CH*2 == 36 ) continue;
	if(BD == 18 && SKI == 2 && CH*2 == 12) continue;
	if(BD == 19 && SKI == 0 && CH*2 == 46) continue;
	tpr->Rebin(rebinN);
	tpr->SetName(title);

	int Nbin = tpr->GetNbinsX();
	vector<double> HG,LG;
	HG.clear();
	LG.clear();
	for(int i = 0 ; i < Nbin ; ++i){
	  double x = tpr->GetBinCenter(i);
	  double y = tpr->GetBinContent(i);
	  
	  if( x > 0 && x < 300 && y != 0 && i%2 == 0){
	    LG.push_back(x);
	    HG.push_back(y);}
	  if( x > 300 && y != 0 && i%5 == 0){
	    LG.push_back(x);
	    HG.push_back(y);}
	}
	
	int np = LG.size();
	TSpline3 *s = new TSpline3("grs",&LG[0],&HG[0],np);

	
	double diff = 0.1;
	double h_start = 0;
        double h_end   = 800;
	TH1D *h_derivative = new TH1D("","",(h_end - h_start)/diff,h_start,h_end);
	
	double xx = MINPOINT;
	for(int i = 0 ; i < h_derivative->GetNbinsX() ; ++i){
	  double x = h_derivative->GetBinCenter(i);
	  if(x <= MINPOINT || x >=750) continue;
	  double y = s->Derivative(xx);
	  h_derivative->SetBinContent(i,y);
	  xx += diff;}

	double start_val = s->Derivative(MINPOINT);
	double tmp_sat = -999,tmp_sat_x,thres;
	thres = start_val - 0.1;
	bool shift = false;
	for(int i = 0 ; i < h_derivative->GetNbinsX() ; ++i){
	  double x = h_derivative->GetBinCenter(i);
	  if( x <= 150 ) continue;
	  double y = h_derivative->GetBinContent(i);
	  if (y < thres){
	    //h_derivative->GetXaxis()->SetRange(x,x+30);
	    double local_max = h_derivative->GetMaximum();
	    if( local_max > thres && !shift){
	      thres -= 0.1;
	      shift = true;
	    }
	    else{
	      tmp_sat   = s->Eval(x);
	      tmp_sat_x = x;
	      break;	    }
	  }}

	//*************************************************************
	//overwrite the sat_x and threshold calculate by 1st derivative
	//*************************************************************
	TH1D *derivative_2nd = TSpline_2nd_deri(h_derivative,&tmp_sat_x);
	tmp_sat   = s->Eval(tmp_sat_x);

	int pad1max  = 3500;
	int pad1xmax = 800;      		

	pad2_sp.SetFillColor(0);
	pad2_sp.SetFillStyle(4000);
	pad2_sp.SetFrameFillStyle(0);
	
	pad1_sp.Draw();
	pad1_sp.cd();
	s->SetLineColor(SKI+2);
	s->SetLineWidth(1.5);
	tpr->Draw();
	tpr->SetXTitle("LG");
	tpr->SetYTitle("HG");
	tpr->GetYaxis()->SetTitleOffset(1.2);
	tpr->SetMaximum(pad1max);
	s->Draw("same");

	
	pad2_sp.Draw();
	pad2_sp.cd();
	
	h_derivative->SetMaximum(20);
	h_derivative->SetMinimum(0);
	h_derivative->SetLineColor(SKI+3);
	h_derivative->SetLineWidth(1.5);
	h_derivative->Draw("Y+");
      	h_derivative->GetXaxis()->SetNdivisions(0);
	h_derivative->GetYaxis()->SetTitle("1st_deri");
	h_derivative->GetYaxis()->SetTitleOffset(0.8);
	h_derivative->Draw("Y+");

	
	//Start to Draw
	pad1_sp.cd();
	//pad1_sp.Draw();
	TLine *Gline2 = new TLine(tmp_sat_x,0,tmp_sat_x,pad1max);
	Gline2->SetLineColor(1);
	Gline2->SetLineWidth(4.8);
	Gline2->SetLineStyle(7);
	Gline2->Draw();

	TLine *Gline3 = new TLine(0,tmp_sat,pad1xmax,tmp_sat);
	Gline3->SetLineColor(2);
	Gline3->SetLineWidth(4.8);
	Gline3->SetLineStyle(7);
	Gline3->Draw();
	
	pad1_sp.Draw();

	c1->Update();
	//if(tmp_sat < 1400 )
	//c1->WaitPrimitive();
	// getchar();
	// c1->SaveAs("data_plus_1stderi.png");

	//Fit with pol1 from 50 to sat point
	linear->SetRange(MINPOINT,tmp_sat_x);
        tpr->Fit(linear,"EMRQ");
	
	p0_ARR[BD][SKI][CH] = linear->GetParameter(0);
	p1_ARR[BD][SKI][CH] = linear->GetParameter(1);
	sat_ARR[BD][SKI][CH] = tmp_sat;
	sat_good[BD][SKI][CH] = true;

	if( sat_ARR[BD][SKI][CH] > 2500 || sat_ARR[BD][SKI][CH] < 1200)
	  sat_good[BD][SKI][CH] = false;
	if( p1_ARR[BD][SKI][CH] < 6 || p1_ARR[BD][SKI][CH] > 10 )
	  sat_good[BD][SKI][CH] = false;

	

	delete s;delete h_derivative;
        delete Gline2;delete Gline3;
	delete derivative_2nd;

      }
    }
  }


  // Calculate avg per chip
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
	if(sat_ARR[BD][SKI][CH] == -999 || sat_good[BD][SKI][CH] == false)
	  continue;
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

  
  sprintf(title,"HGLG_sat_Spline_%iGeV.root",labelE);
  TFile outf(title,"recreate");
  TTree *outtree = new TTree("tree","tree");
  int layerID,skirocID,channelID;
  double m_p0,m_p1,m_sat,m_tpr_entry;
  bool   m_goodsat;

  outtree->Branch("layerID",&layerID,"layerID/I");
  outtree->Branch("skirocID",&skirocID,"skirocID/I");
  outtree->Branch("channelID",&channelID,"channelID/I");
  outtree->Branch("p0",&m_p0,"p0/D");
  outtree->Branch("p1",&m_p1,"p1/D");
  outtree->Branch("p_saturation",&m_sat,"p_saturation/D");
  outtree->Branch("good_saturation",&m_goodsat,"good_saturation/O");
  outtree->Branch("tpr_entry",&m_tpr_entry,"tpr_entry/D");
  
  
  for(int BD = 0 ;BD < MAXBD ; ++BD){
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	layerID  = BD;
	skirocID = SKI;
	channelID = CH*2;
	if( sat_good[BD][SKI][CH] == true){
	  m_p0  = p0_ARR[BD][SKI][CH];
	  m_p1  = p1_ARR[BD][SKI][CH];
	  m_sat = sat_ARR[BD][SKI][CH];}
	else{
	  m_p0  = p0_avg[BD][SKI];
	  m_p1  = p1_avg[BD][SKI];
	  m_sat = sat_avg[BD][SKI];}
	m_tpr_entry = tpr_entry[BD][SKI][CH];
	m_goodsat = sat_good[BD][SKI][CH];
	outtree->Fill();	
      }
    }
  }
  
  TH1D* hp0 = qualify(p0_ARR,2);
  TH1D* hp1 = qualify(p1_ARR,0);
  TH1D* hsat = qualify(sat_ARR,1);
  hp0->SetTitle("p0_qualify");
  hp0->SetName("p0_qualify");
  hp1->SetTitle("p1_qualify");
  hp1->SetName("p1_qualify");
  hsat->SetTitle("sat_qualify");
  hsat->SetName("sat_qualify");

  
  sprintf(title,"plot_out/%iGeV/p0_hist_spline.png",labelE);
  hp0 ->Draw(); c1->Update(); c1->SaveAs(title);
  sprintf(title,"plot_out/%iGeV/p1_hist_spline.png",labelE);
  hp1 ->Draw(); c1->Update(); c1->SaveAs(title);
  sprintf(title,"plot_out/%iGeV/sat_hist_spline.png",labelE);
  hsat->Draw(); c1->Update(); c1->SaveAs(title);

  outf.Write();
  outf.Close();

  
}
TH1D* fitter::TSpline_2nd_deri(TH1D *h_deri, double *sat_x){
  
  int Nbin = h_deri->GetNbinsX();
  vector<double> HG,LG;
  for(int i = 0 ; i < Nbin ; ++i){
    double x = h_deri->GetBinCenter(i);
    double y = h_deri->GetBinContent(i);
    if( x > MINPOINT && x < 300 ){
      LG.push_back(x);
      HG.push_back(y);}  
}
  int np = LG.size();
  TSpline3 *s = new TSpline3("grs",&LG[0],&HG[0],np);
  
  TH1D *h_2ndderi = new TH1D("","",Nbin*10,0,800);
  for(int i = 0 ; i < Nbin*10 ; ++i){
    
    double x = h_2ndderi->GetBinCenter(i);
    if( x > MINPOINT && x < 350 ){
      double y = s->Derivative(x);
      h_2ndderi->SetBinContent(i,y);    }
    else
      h_2ndderi->SetBinContent(i,0);
  }
  
  double thres = -0.04;
  double sat_x_tmp = 0;
  
  while(true){
    sat_x_tmp = first_P_lower_thres(h_2ndderi,thres,180,350);
    if(thres >= 0 || (sat_x_tmp > 150 && sat_x_tmp < 300) )
      break;
    thres += 0.004;}
    
  TPad pad1_t("pad1_t","",0,0,1,1);
  TPad pad2_t("pad2_t","",0,0,1,1);

  pad2_t.SetFillColor(0);
  pad2_t.SetFillStyle(4000);
  pad2_t.SetFrameFillStyle(0);

  pad1_t.Draw();
  pad1_t.cd();
  h_deri->Draw();
  h_deri->SetXTitle("LG(ADC)");
  h_deri->SetYTitle("1st_deri");
    
  pad2_t.Draw();
  pad2_t.cd();

  double canvas_xmin = h_2ndderi->GetXaxis()->GetXmin();
  double canvas_xmax = h_2ndderi->GetXaxis()->GetXmax();
  double canvas_ymin = h_2ndderi->GetMinimum();
  double canvas_ymax = h_2ndderi->GetMaximum();
  h_2ndderi->Draw("Y+");
  h_2ndderi->SetYTitle("2nd_deri");
  h_2ndderi->GetYaxis()->SetTitleOffset(1.2);
  h_2ndderi->SetLineColor(1);

  //cout << "thres_method : thres = " << thres << ", x = "<< sat_x_tmp << endl;
  vector<double> zeros = findzeros(h_2ndderi);
  for(int i = 0 ;i < (int)zeros.size() ; ++i){
    if(zeros[i] > sat_x_tmp && zeros[i] < 250){
      sat_x_tmp = zeros[i];
    }
  }
  
  
  TLine *line2 = new TLine(canvas_xmin,thres,canvas_xmax,thres);
  line2->SetLineColor(7);
  line2->SetLineWidth(4.8);
  line2->SetLineStyle(7);
  line2->Draw("same");  

  TLine *line3 = new TLine(sat_x_tmp,canvas_ymin,sat_x_tmp,canvas_ymax);
  line3->SetLineColor(2);
  line3->SetLineWidth(4.8);
  line3->SetLineStyle(7);
  line3->Draw("same");
  
  
  //c1->Update();
  //c1->SaveAs("deri.png");
  //c1->WaitPrimitive();
  
  delete s;delete line2;delete line3;


  *sat_x = sat_x_tmp;  
  return h_2ndderi;
}

double fitter::first_P_lower_thres(TH1D* hist,double thres,double lowerbond,double upperbond){

  bool first_time = false;
  double sat_x = 0;
  int Nbin = hist->GetNbinsX();
  
  for(int i = 0 ; i < Nbin ; ++i){
    double x = hist->GetBinCenter(i);
    double y = hist->GetBinContent(i);
    if( y < thres && !first_time && x > lowerbond && x < upperbond){
	first_time = true;
	sat_x = x;
    }
  }
    
  return sat_x;
}

vector<double> fitter::findzeros(TH1D* hist){
  vector<double> zeros;
  vector<double> target_b,target_x,target_y;
  int Nbin = hist->GetNbinsX();
 
  bool first_flag = false;
  double lastx;
  int cut = 0;
  for(int i = 0 ; i < Nbin ; ++i){
    double x = hist->GetBinCenter(i);
    double y = hist->GetBinContent(i);
    if( y == 0 ) continue;
    
    if( abs(y) < 0.0001 ) {
      if(!first_flag){
	lastx = x;
	first_flag = true;}      
      if(x - lastx > 10) {
	cut++;
	if(cut >= 1){
	  if(target_x.size() != 0){
	    if(target_x[(target_x.size()+1)/2] > 10){
	      zeros.push_back(target_x[(target_x.size()+1)/2]);
	      target_x.clear();
	    }
	  }
	}
      }
      lastx = x;
      target_x.push_back(x);
    }
  }
  if(target_x.size() != 0){
    zeros.push_back(target_x[(target_x.size()+1)/2]);
    target_x.clear();
  }

  //for(int i = 0 ; i < (int)zeros.size() ; ++i){
  //  cout << "No." << i << " = "<< zeros.at(i) << endl;
  // }
  
  return zeros;
}

void fitter::pol4(TProfile *tpr){
  char title[50];
  int N_of_par  = 4;
  sprintf(title,"pol%i",N_of_par);
  TF1    *mytry = new TF1("",title,50,600);
  double pars[N_of_par];
  double pars_der[N_of_par-1];
  double pars_2der[N_of_par-2];
	
  //TSpline3 *s = new TSpline3("grs",gr);
  tpr->Draw();
  tpr->SetMaximum(4000);
  tpr->Fit(mytry,"EMR0");
  for(int i = 0 ; i < N_of_par; ++i){
    pars[i] = mytry->GetParameter(i);
    if( i == 0) continue;
    pars_der[i-1] = mytry->GetParameter(i)*(i);
    if( i <= 1) continue;
    pars_2der[i-2] = mytry->GetParameter(i)*(i)*(i-1);
  }
  sprintf(title,"%f+x*%f+x^2*%f",pars_der[0],pars_der[1],pars_der[2]);
  TF1    *mytry_deri = new TF1("",title,50,600);
  mytry_deri->SetTitle("");

  sprintf(title,"10*(%f+x*%f)",pars_2der[0],pars_2der[1]);
  TF1    *mytry_2deri = new TF1("",title,50,600);
  mytry_2deri->SetTitle("");


  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(4000);
  pad2->SetFrameFillStyle(0);
	
  pad1->Draw();
  pad1->cd();
  tpr->Draw();
  tpr->SetXTitle("LG");
  tpr->SetYTitle("HG");
  tpr->Draw();
  
  mytry->Draw("same");
  mytry->GetXaxis()->SetLabelOffset(999);
  mytry->GetXaxis()->SetLabelSize(0);
	
  pad2->Draw();
  pad2->cd();
	
  mytry_deri->SetLineWidth(1.2);
  mytry_deri->SetLineColor(8);
  mytry_deri->GetXaxis()->SetNdivisions(0);
  mytry_deri->GetYaxis()->SetTitle("1st_deri");
  mytry_deri->Draw("Y+");

  mytry_2deri->SetLineWidth(1.2);
  mytry_2deri->SetLineColor(4);
  mytry_2deri->GetXaxis()->SetNdivisions(0);
  mytry_2deri->Draw("same");

  double tmp_sat = mytry_deri->GetMaximumX(50,200);
  cout << tmp_sat << endl;
  TLine *Gline = new TLine(tmp_sat,0,tmp_sat,4000);
  Gline->SetLineColor(1);
  Gline->SetLineWidth(4.8);
  Gline->SetLineStyle(7);
	  

  pad1->cd();
  Gline->Draw("same");
  c1->Update();
  c1->WaitPrimitive();
}

double fitter::spline_4nodes(double *x, double *par){
   /*Fit parameters:
   par[0-3]=X of nodes (to be fixed in the fit!)
   par[4-7]=Y of nodes
   par[8-9]=first derivative at begin and end (to be fixed in the fit!)
   */
   Double_t xx = x[0];

   Double_t xn[4] = { par[0], par[1], par[2], par[3] };
   Double_t yn[4] = { par[4], par[5], par[6], par[7] };

   Double_t b1 = par[8];
   Double_t e1 = par[9];

   TSpline3 sp3("sp3", xn, yn, 4, "b1e1", b1, e1);

   return sp3.Eval(xx);
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
