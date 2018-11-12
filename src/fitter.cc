#include "fitter.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include <algorithm>
#include <fstream>

int MINPOINT = 5;
int MAXPOINT = 40;

fitter::fitter(){

}
fitter::~fitter(){
}
void fitter::DEBUG(){
  
  char title[50];
  for(int labelE = 0 ; labelE < 4 ; ++labelE){
    
  
    if(labelE % 2 == 0)
      sprintf(title,"root_result/0916update/inj.root");
    else
      sprintf(title,"root_result/0916update/%iGeV.root",100);

    TFile f(title);

    TProfile *tpr = (TProfile *)f.Get("Board_7/HGLG_chip2_ch44");
    sprintf(title,"Board_%i/HGLG_chip%i_ch%i",0,1,2);
    tpr = (TProfile *)f.Get(title);
    

    f.Close();
  }
}

void fitter::fit(int labelE){
  if(!(labelE == 10 || labelE == 30 || labelE == 50 || labelE == 80
       || labelE == 100 || labelE == 150 || labelE == -1)) {
    cout << "invalid energy!" << endl;
    return;}
  char title[50];
  if(labelE == -1)
    sprintf(title,"root_result/0916update/inj.root");
  else if(labelE == 100)
    sprintf(title,"TPro_ele.root");
  else
    sprintf(title,"root_result/0916update/%iGeV.root",labelE);
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
  root_logon();
  gStyle->SetOptStat(0);
  TProfile *tpr = (TProfile *)f.Get("Board_7/HGLG_chip2_ch44");
  //  TProfile *tpr = new TProfile("","",200,0,800,0,4000);
  TH1D *h1 = new TH1D("","",tpr->GetNbinsX(),0,800);
  int rebinN = 2;
  h1->Rebin(rebinN);
  bool savepng = 0;
  int stop_and_look = -1;
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
	
	tpr->Rebin(rebinN);
	
	int nentry = tpr->GetEntries();
	if( nentry < 1000 && labelE != -1 ) continue;
	
	sprintf(title,"Board%i_HGLG_chip%i_ch%i",BD,SKI,true_ch);
	tpr->SetName(title);
	tpr->SetTitle(title);
	
	//tpr->Draw();
	//c1->Update();
	//getchar();
 
	tpr->Fit(sat_fit,"QEMR0");
	tpr->Draw();
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
        tpr->Fit(sat_fit,"QEMR0");
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
	tpr->Fit(sat_fit,"QEMR0");
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
			
	ratio_plot(tpr,sat_fit,h1);
	TLine *Gline = new TLine(0,sat_point,tpr->GetXaxis()->GetXmax(),sat_point);
	Gline->SetLineColor(1);
	Gline->SetLineWidth(4.8);
	Gline->SetLineStyle(7);
	pad1->cd();
	Gline->Draw("same");
	
	if(stop_and_look == 0){
	  c1->cd();
	  c1->Update();
	  c1->WaitPrimitive();}
	
	if(savepng)
	  c1->SaveAs("step3.png");
	
	sat_fit->SetRange(MINPOINT,sat_point_x);
	tpr->Fit(sat_fit,"QEMR0");
	ratio_plot(tpr,sat_fit,h1);

	if(stop_and_look == 0){
	  pad1->cd();}
	
	pad1->cd();
	Gline->Draw("same");
	
	c1->cd();
	c1->Update();
	//c1->WaitPrimitive();

	if(stop_and_look == 0){
	  c1->cd();
	  c1->Update();
	  c1->WaitPrimitive();}

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

  sprintf(title,"HGLG_sat_linear.root");
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

  
  // sprintf(title,"plot_out/%iGeV/p0_hist.png",labelE);
  // hp0 ->Draw(); c1->Update(); c1->SaveAs(title);
  // sprintf(title,"plot_out/%iGeV/p1_hist.png",labelE);
  // hp1 ->Draw(); c1->Update(); c1->SaveAs(title);
  // sprintf(title,"plot_out/%iGeV/sat_hist.png",labelE);
  // hsat->Draw(); c1->Update(); c1->SaveAs(title);

  outf.Write();
  outf.Close();

  ofstream outcalib("linear_Calib.txt");
  outcalib << "Layer  ASIC_ID  Channel  p0  p1  sat  true"
	   <<"  totalcalibch" << endl;
  for(int BD = 0 ;BD < MAXBD ; ++BD){
    int count = 0;
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      for(int CH = 0 ; CH < MAXCH ; ++CH){ count+= sat_good[BD][SKI][CH]; }
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	outcalib << BD << "\t" << SKI << "\t" << CH*2 << "\t"
		 << p1_ARR[BD][SKI][CH] << "\t" << p0_ARR[BD][SKI][CH]
		 << "\t" << sat_ARR[BD][SKI][CH] << "\t"
		 << sat_good[BD][SKI][CH] << "\t" << count << endl;
      }
    }
  }
  outcalib.close();
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
void fitter::Find_high(TH1D* h1,double* highx,double* highy,double lowerbound ,double upperbound){
  
  int Nbin = h1->GetNbinsX ();
  double maxx = -1,maxy = -1;
  for(int i = 0 ; i < Nbin ; ++i){
    double x = h1->GetBinCenter(i);
    double y = h1->GetBinContent(i);
    if(y > maxy && x > lowerbound && x < upperbound){
      maxx = x;
      maxy = y;
    }
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
      bool remove = (abs(res/HG_vec.at(i)) > 0.2) ? true  : false ;
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
       || labelE == 100 || labelE == 150 || labelE == -1)) {
    cout << "invalid energy!" << endl;
    return;}
  
  char title[50];
  if(labelE == -1)
    sprintf(title,"root_result/0916update/inj.root");
  else if(labelE == 100){
    sprintf(title,"TPro_ele.root");}
  else
    sprintf(title,"root_result/0916update/%iGeV.root",labelE);

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
  //bool savepng = 0;
  //int stop_and_look = -1;
  //gROOT->SetBatch(kTRUE);
  
  TF1 *linear = new TF1("","[0]+x*[1]",MINPOINT,150);
  linear->SetParLimits(0,-50,50);
  TPad pad1_sp("pad1_sp","",0,0,1,1);
  TPad pad2_sp("pad2_sp","",0,0,1,1);

  //Loop over all channels (Get TProfile)
  for(int BD = 0 ;BD < MAXBD ; ++BD){
    cout << "HGLG BD "<< BD << endl;
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      cout << "SKI "<< SKI << endl;
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	cout << "CH " << CH << endl;
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
	if(labelE != -1 && tpr->GetEntries() < 1500) continue;

	sprintf(title,"Board%i_HGLG_chip%i_ch%i",BD,SKI,true_ch);
	tpr->Rebin(rebinN);
	tpr->SetName(title);
	tpr->SetTitle(title);

	//Get TProfile x,y to vectors for TSpline input
	int Nbin = tpr->GetNbinsX();
	vector<double> HG,LG;
	HG.clear();
	LG.clear();
	int get_pt_for_each = (labelE == -1) ? 2 : 8;
	for(int i = 0 ; i < Nbin ; ++i){
	  double x = tpr->GetBinCenter(i);
	  double y = tpr->GetBinContent(i);
	  
	  if( x > 0 && x < 300 && y != 0 && i % get_pt_for_each == 0){
	    LG.push_back(x);
	    HG.push_back(y);}
	  if( x > 300 && y != 0 && i%5 == 0){
	    LG.push_back(x);
	    HG.push_back(y);}
	}

	
	int np = LG.size();
	if(np == 0) continue;
	TSpline3 *s = new TSpline3("grs",&LG[0],&HG[0],np);

	//Calculate derivative for spline shown by TH1
	double diff = 0.1;
	double h_start = 0;
        double h_end   = 800;
	TH1D *h_derivative = new TH1D("","",(h_end - h_start)/diff,h_start,h_end);
	
	double xx = MINPOINT;
	for(int i = 0 ; i < h_derivative->GetNbinsX() ; ++i){
	  double x = h_derivative->GetBinCenter(i);
	  if(x <= MINPOINT || x >=350) continue;
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
	TH1D *derivative_2nd = new TH1D();
	TSpline_2nd_deri(*derivative_2nd,h_derivative,&tmp_sat_x);
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

	//Start to Draw	
	TLine *Gline2 = new TLine(tmp_sat_x,0,tmp_sat_x,pad1max);
	Gline2->SetLineColor(1);
	Gline2->SetLineWidth(4.8);
	Gline2->SetLineStyle(7);
	Gline2->Draw("same");

	TLine *Gline3 = new TLine(0,tmp_sat,pad1xmax,tmp_sat);
	Gline3->SetLineColor(2);
	Gline3->SetLineWidth(4.8);
	Gline3->SetLineStyle(7);
	Gline3->Draw("same");
	
	//if(tmp_sat < 1400 )
	// getchar();
	// c1->SaveAs("data_plus_1stderi.png");

	//Fit with pol1 from 50 to sat point

	double setsat = 1400;
	double diff_y_sat = 100;
	for(int i = 0 ; i < Nbin ; ++i){
	  double x = tpr->GetBinCenter(i);
	  double y = tpr->GetBinContent(i);
	  if( abs(y-setsat) < diff_y_sat ){
	    tmp_sat_x = x;
	    diff_y_sat = abs(y-setsat); }
	    
	}

	
	linear->SetRange(MINPOINT,tmp_sat_x);
	linear->SetLineColor(7);
        tpr->Fit(linear,"EMRQ");


	pad2_sp.Draw();
	pad2_sp.cd();
	
	h_derivative->SetMaximum(20);
	h_derivative->SetMinimum(0);
	h_derivative->SetLineColor(SKI+3);
	h_derivative->SetLineWidth(1.5);
      	h_derivative->GetXaxis()->SetNdivisions(0);
	h_derivative->GetYaxis()->SetTitle("1st_deri");
	h_derivative->GetYaxis()->SetTitleOffset(0.8);
	h_derivative->Draw("Y+");

	c1->cd();
	
	p0_ARR[BD][SKI][CH] = linear->GetParameter(0);
	p1_ARR[BD][SKI][CH] = linear->GetParameter(1);
	sat_ARR[BD][SKI][CH] = tmp_sat;
	sat_good[BD][SKI][CH] = true;

	if( sat_ARR[BD][SKI][CH] > 2500 || sat_ARR[BD][SKI][CH] < 1200)
	  sat_good[BD][SKI][CH] = false;
	if( p1_ARR[BD][SKI][CH] < 6 || p1_ARR[BD][SKI][CH] > 10 )
	  sat_good[BD][SKI][CH] = false;

	c1->Update();
	//c1->WaitPrimitive();
	
	
	delete s;delete h_derivative;
        delete Gline2;delete Gline3;
	//delete derivative_2nd;

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

  
  sprintf(title,"HGLG_sat_Spline.root");
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

	opt_val[(BD*4+SKI)*32+CH].L_ID = BD;
	opt_val[(BD*4+SKI)*32+CH].S_ID = SKI;
	opt_val[(BD*4+SKI)*32+CH].C_ID = CH*2;
	opt_val[(BD*4+SKI)*32+CH].L2HT = m_sat;
	opt_val[(BD*4+SKI)*32+CH].L2H  = m_p1;
	opt_val[(BD*4+SKI)*32+CH].HLTYPE = m_goodsat;

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

  
  // sprintf(title,"plot_out/%iGeV/p0_hist_spline.png",labelE);
  // hp0 ->Draw(); c1->Update(); c1->SaveAs(title);
  // sprintf(title,"plot_out/%iGeV/p1_hist_spline.png",labelE);
  // hp1 ->Draw(); c1->Update(); c1->SaveAs(title);
  // sprintf(title,"plot_out/%iGeV/sat_hist_spline.png",labelE);
  // hsat->Draw(); c1->Update(); c1->SaveAs(title);

  outf.Write();
  outf.Close();

  delete linear;
  f.Close();
}
void fitter::TSpline_2nd_deri(TH1D& h_2nd_deri,TH1D *h_deri, double *sat_x){
  
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
  
  TH1D h_2ndderi("","",Nbin*10,0,800);
  for(int i = 0 ; i < Nbin*10 ; ++i){
    
    double x = h_2ndderi.GetBinCenter(i);
    if( x > MINPOINT && x < 350 ){
      double y = s->Derivative(x);
      h_2ndderi.SetBinContent(i,y);    }
    else
      h_2ndderi.SetBinContent(i,0);
  }
  
  double thres = -0.04;
  double sat_x_tmp = 0;
  
  while(true){
    sat_x_tmp = first_P_lower_thres(&h_2ndderi,thres,180,350);
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

  double canvas_xmin = h_2ndderi.GetXaxis()->GetXmin();
  double canvas_xmax = h_2ndderi.GetXaxis()->GetXmax();
  double canvas_ymin = h_2ndderi.GetMinimum();
  double canvas_ymax = h_2ndderi.GetMaximum();
  h_2ndderi.Draw("Y+");
  h_2ndderi.SetYTitle("2nd_deri");
  h_2ndderi.GetYaxis()->SetTitleOffset(1.2);
  h_2ndderi.SetLineColor(1);

  //cout << "thres_method : thres = " << thres << ", x = "<< sat_x_tmp << endl;
  vector<double> zeros = findzeros(&h_2ndderi);
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
  h_2nd_deri = h_2ndderi;
  //return h_2ndderi;
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
  double pars_der[N_of_par-1];
  double pars_2der[N_of_par-2];
	
  //TSpline3 *s = new TSpline3("grs",gr);
  tpr->Draw();
  tpr->SetMaximum(4000);
  tpr->Fit(mytry,"EMR0");
  for(int i = 0 ; i < N_of_par; ++i){
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

void fitter::ratio_plot(TProfile *tpr,TF1 *fit,TH1D *hratio,string X_title,string Y_title){
  char label_char[50];
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

  sprintf(label_char,"#frac{%s - Fit}{Fit}",Y_title.c_str());
  hratio->GetYaxis()->SetTitle(label_char);
  hratio->GetYaxis()->SetTitleOffset(0.35);
  hratio->GetYaxis()->SetLabelFont(43); 
  hratio->GetYaxis()->SetLabelSize(15);

  axis = hratio->GetYaxis();
  axis->SetNdivisions(404);
  axis->Draw();

  axis = hratio->GetXaxis();
  axis->SetTitle(X_title.c_str());
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
  c1->cd();
  c1->Update();
  
  //getchar();
}


void fitter::fit_LGTOT(int labelE ){
  if(!(labelE == 10 || labelE == 30 || labelE == 50 || labelE == 80
       || labelE == 100 || labelE == 150 || labelE == -1)) {
    cout << "invalid energy!" << endl;
    return;}

  char title[50];
  if(labelE == -1)
    sprintf(title,"root_result/0916update/inj.root");
  else if(labelE == 100){
    sprintf(title,"TPro_ele.root",labelE);  }
  else
    sprintf(title,"root_result/0916update/%iGeV.root",labelE);
  TFile f(title);


  double Offset_avg[MAXBD][MAXSKI],Gain_avg[MAXBD][MAXSKI];
  double Offset[MAXBD][MAXSKI][MAXCH],Gain[MAXBD][MAXSKI][MAXCH];
  bool   Good_fit[MAXBD][MAXSKI][MAXCH];
  
  for(int BD = 0 ;BD < MAXBD ; ++BD){
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      Offset_avg[BD][SKI] = 0;
      Gain_avg  [BD][SKI] = 0;
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	Offset  [BD][SKI][CH] = -1;
	Gain    [BD][SKI][CH] = -1;
	Good_fit[BD][SKI][CH] = false;      }    }  }
  
  TF1 *sat_fit   = new TF1("","pol1");

  root_logon();
  gStyle->SetOptStat(0);
  TProfile *tpr = (TProfile *)f.Get("Board_7/LGTOT_chip2_ch44");
 
  TH1D *h_gain = new TH1D("","",100,0,10);
  TH1D *h1 = new TH1D("","",tpr->GetNbinsX(),0,800);
  int rebinN = 1;
  h1->Rebin(rebinN);
  //bool savepng = 0;
  //int stop_and_look = -1;
  //gROOT->SetBatch(kTRUE);

  MINPOINT = 300;
  MAXPOINT = 500;

  
  for(int BD = 0 ;BD < MAXBD ; ++BD){
    cout << "LGTOT BD "<< BD << endl;
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      cout << "SKI "<< SKI << endl;
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	//cout << "CH " << CH << endl;
	sat_fit->SetRange(MINPOINT,MAXPOINT);
        
	int true_ch = CH*2;
	sprintf(title,"Board_%i/LGTOT_chip%i_ch%i",BD,SKI,true_ch);
	tpr = (TProfile *)f.Get(title);
	if(tpr == NULL) continue;
	//tpr->Draw();
	if( BD == 16 && SKI == 2 && CH == 27 ){ continue; }
	if( BD == 23 && SKI == 0 && CH == 23 ){ continue; }

	tpr->Rebin(rebinN);
	 
	int nentry = tpr->GetEntries();
	if( nentry < 100 && labelE != -1 ) continue;
	//if(labelE != -1 && BD >= ) continue;

	//cout << BD <<" ," << SKI << ", " << CH << endl;
	sprintf(title,"Board%i_LGTOT_chip%i_ch%i",BD,SKI,true_ch);
	tpr->SetName(title);
	tpr->SetTitle(title);
 
	//if(BD >= 5 && BD <= 15)
	//  getchar();

	int Nbin = tpr->GetNbinsX();
	vector<double> LG,TOT;
	LG.clear();
	TOT.clear();
	
	// Fill spline(LG vs TOT) :separate injection and TB data

	bool AllTotlessthan100 = true;
	if(labelE == -1){
	  double continue_check_y = -1;
	  bool first_non_zero = false;
	  int point_counter = 0;
	  int required_value = 10;
	  for(int i = 1 ; i < Nbin ; ++i){
	    double x = tpr->GetBinCenter(i);
	    double y = tpr->GetBinContent(i);
	    //cout << "x = " << x << ", y = " << y << endl;
	    if(!first_non_zero){
	      if(y == 0) continue;
	      continue_check_y = y;
	      first_non_zero = true;	  }
	    else{
	      if( abs(y - continue_check_y) > continue_check_y*0.15 )
		continue;
	      continue_check_y = y;
	      point_counter++;
	      if(point_counter == required_value)
		point_counter = 0;
	  
	      if( x > 100 && point_counter == 0 ){
		AllTotlessthan100 = false;
		if( x < 300 && y > 1000){
		  point_counter-=1;
		  continue;}
		TOT.push_back(x);
		LG.push_back(y); }
	    }
	  }
	}
	
	else{
	  vector<int> ithbin;
	  for(int i = 1 ; i < Nbin ; ++i){
	    double x = tpr->GetBinCenter(i);
	    double y = tpr->GetBinContent(i);
	    if(y > 5 && x > 100){
	      AllTotlessthan100 = false;
	      ithbin.push_back(i);}
	  }
	  int spline_points = 10;
	  float factor = (float)ithbin.size()/spline_points;
	  for(int i = 0 ; i < spline_points ; ++i){
	    double x = tpr->GetBinCenter(ithbin[(int)i*factor]);
	    double y = tpr->GetBinContent(ithbin[(int)i*factor]);
	    TOT.push_back(x);
	    LG.push_back(y);
	    //cout << x << ", " << y << endl;
	  }
	}
	if(AllTotlessthan100) continue;	
	tpr->Draw();
	c1->Update();
	//getchar();
	for(int i = 1 ; i < (int)LG.size()-1 ; ++i){
	  if(TOT[i-1] > TOT[i] && TOT[i+1] > TOT[i]){
	    LG.erase(LG.begin() + i);
	    TOT.erase(TOT.begin() + i);
	  }
	}
	//cout << LG.size() << endl;
	//getchar();
	int np = LG.size();
	if(np == 0) continue;
	TSpline3 *s = new TSpline3("grs",&TOT[0],&LG[0],np);
	
	double diff = 0.1;
	double h_start = 0;
	double h_end   = 800;

	TH1D *h_derivative = new TH1D("","",(h_end - h_start)/diff,h_start,h_end);
	
	double xx = h_start;
	for(int i = 0 ; i < h_derivative->GetNbinsX() ; ++i){
	  double x = h_derivative->GetBinCenter(i);
	  double y = s->Derivative(xx);
	  if (y <=6 && x < s->GetXmax() && x > s->GetXmin())
	    h_derivative->SetBinContent(i,y);
	  xx += diff;	}
	h_derivative->SetMaximum(10);
	h_derivative->SetMinimum(0);
	Draw_Spline_and_1stderi(*tpr,*s,*h_derivative);
	//c1->Update();
	//c1->WaitPrimitive();
	
	double localmax_x,localmax_y;
	Find_high(h_derivative,&localmax_x,&localmax_y,s->GetXmin()+20,s->GetXmin()+300);

	
	//Prevention if fail
	if(localmax_x < 10) continue;
	
	sat_fit->SetRange(localmax_x*0.95,localmax_x*1.05);
	tpr->Fit(sat_fit,"QEMR0");
	//cout << 	localmax_x << endl;
	int fit_pt = sat_fit->GetNumberFitPoints();
	//cout << fit_pt << endl;
	if(fit_pt < 5)
	  continue;
	
	double tmp_gain = sat_fit->GetParameter(1);
	h_gain->Fill(tmp_gain);

	TH1D *h_res = new TH1D("","",tpr->GetNbinsX(),0,800);
	for(int i = 0 ; i < tpr->GetNbinsX () ; ++i){
	  double x = tpr->GetBinCenter(i);
	  double y = tpr->GetBinContent(i);
	  if(x == 0 || y == 0) continue;
	  //cout << "x = "<< x << ", y = " << y << endl;
	  double res = ( y - sat_fit->Eval(x)) / sat_fit->Eval(x);
	  double error = tpr->GetBinError(i)/(x*sat_fit->Eval(x));
	  h_res->SetBinContent(i,res);
	  h_res->SetBinError(i,error);	 }

	ratio_plot(tpr,sat_fit,h_res,"TOT","LG");
	
	//c1->WaitPrimitive();
	//c1->SaveAs("LGTOT_spline_fit.png");
	//getchar();

	//cout << localmax_x << " - " <<  localmax_y/tmp_gain << endl;
	Offset  [BD][SKI][CH] = localmax_x - sat_fit->Eval(localmax_x)/tmp_gain;
	Gain    [BD][SKI][CH] = tmp_gain;

	delete s; delete h_derivative; delete h_res;
	
	
      }
    }
  }

  for(int BD = 0 ;BD < MAXBD ; ++BD){
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      int counter = 0;
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	if(Offset  [BD][SKI][CH] < 100 || Offset  [BD][SKI][CH] > 600)
	  continue;
	if(Gain  [BD][SKI][CH] < 3 || Gain  [BD][SKI][CH] > 7)
	  continue;
	Good_fit [BD][SKI][CH] = true;
	counter++;
	Offset_avg[BD][SKI] += Offset  [BD][SKI][CH];
	Gain_avg  [BD][SKI] += Gain    [BD][SKI][CH];
      }
      Offset_avg[BD][SKI] /= counter;
      Gain_avg  [BD][SKI] /= counter;
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	opt_val[(BD*4+SKI)*32+CH].LTTYPE = Good_fit [BD][SKI][CH];
	if(Good_fit [BD][SKI][CH]){
	  opt_val[(BD*4+SKI)*32+CH].T2L  = Gain[BD][SKI][CH];
	  opt_val[(BD*4+SKI)*32+CH].TOFF = Offset[BD][SKI][CH];}
	else if(counter == 0){
	  opt_val[(BD*4+SKI)*32+CH].T2L  = 5.;
	  opt_val[(BD*4+SKI)*32+CH].TOFF = 180;
	}
	else{
	  opt_val[(BD*4+SKI)*32+CH].T2L  = Gain_avg[BD][SKI];
	  opt_val[(BD*4+SKI)*32+CH].TOFF = Offset_avg[BD][SKI];
	}
      }      
    }
  }
  
  //h_gain->Draw();
  //c1->WaitPrimitive();
    
  // sprintf(title,"Spline_%iGeV.root",labelE);
  // TFile outf(title,"recreate");
  // TTree *outtree = new TTree("tree","tree");
  // int layerID,skirocID,channelID;
  // double m_p0,m_p1,m_sat,m_tpr_entry;
  // bool   m_goodsat;

  // outtree->Branch("layerID",&layerID,"layerID/I");
  // outtree->Branch("skirocID",&skirocID,"skirocID/I");
  // outtree->Branch("channelID",&channelID,"channelID/I");
  // outtree->Branch("p0",&m_p0,"p0/D");
  // outtree->Branch("p1",&m_p1,"p1/D");
  // outtree->Branch("p_saturation",&m_sat,"p_saturation/D");
  // outtree->Branch("good_saturation",&m_goodsat,"good_saturation/O");
  // outtree->Branch("tpr_entry",&m_tpr_entry,"tpr_entry/D");
 
  delete sat_fit; delete h_gain; delete h1;  
  f.Close();  
}

void fitter::Draw_Spline_and_1stderi(TProfile& tpr, TSpline3 &s, TH1D& h_deri){
  TPad pad1_sp("pad1_sp","",0,0,1,1);
  TPad pad2_sp("pad2_sp","",0,0,1,1);

  pad2_sp.SetFillColor(0);
  pad2_sp.SetFillStyle(4000);
  pad2_sp.SetFrameFillStyle(0);
	
  pad1_sp.Draw();
  pad1_sp.cd();
  s.SetLineColor(2);
  s.SetLineWidth(1.5);
  tpr.Draw();
  tpr.SetXTitle("TOT");
  tpr.SetYTitle("LG");
  tpr.GetYaxis()->SetTitleOffset(1.2);
  s.Draw("same");

  pad2_sp.Draw();
  pad2_sp.cd();
  
  h_deri.SetLineColor(3);
  h_deri.SetLineWidth(1.5);
  h_deri.GetXaxis()->SetNdivisions(0);
  h_deri.GetYaxis()->SetTitle("1st_deri");
  h_deri.GetYaxis()->SetTitleOffset(0.8);
  h_deri.Draw("Y+");
  
  c1->cd();
  c1->Update();
  //c1->SaveAs("LGTOT_splineD1.png");
  //getchar();

  //getchar();
}

void fitter::fit_output(int labelE){

  gROOT->SetBatch(kTRUE);
  cout << "Starting HGLG " << labelE << " fitting... " << endl;
  fit_spline(labelE);

  cout << "Starting LGTOT " << labelE << " fitting... " << endl;
  fit_LGTOT(labelE);
  
  ofstream calib_result;
  if(labelE == -1)
    calib_result.open("Chia-hung_Inj_Calib.txt");
  else
    calib_result.open("Chia-hung_TB_Calib.txt");
  calib_result << "Layer  Module_ID  ASIC_ID  Channel  ADC_To_MIP  LowGain_To_HighGain_Transition  LowGain_To_HighGain_Conversion  TOT_To_LowGain_Transition  TOT_To_LowGain_Conversion  TOT_Offset  HLType LTType\n";
  int layer_to_moduleID[28] = {  78,  90,  89,  88,  77,
				 85,  84,  32,  69,  79,
				 76,  83,  70,  73,  86,
				 87,  82,  72,  67,  65,
				 35,  36,  44,  51, 142,
				143, 145, 144 };

  for(int BD = 0 ;BD < MAXBD ; ++BD){
    for(int SKI = 0 ; SKI < MAXSKI ; ++SKI){
      for(int CH = 0 ; CH < MAXCH ; ++CH){
	output O = opt_val[(BD*4+SKI)*32+CH];
	O.M_ID = layer_to_moduleID[BD];
	calib_result << O.L_ID << "\t" << O.M_ID << "\t" << O.S_ID << "\t"
		     << O.C_ID << "\t" << O.A2M  << "\t" << O.L2HT << "\t"
		     << O.L2H  << "\t" << O.T2LT << "\t" << O.T2L  << "\t"
		     << O.TOFF << "\t" << O.HLTYPE << "\t" << O.LTTYPE << endl;
      }
    }
  }
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


