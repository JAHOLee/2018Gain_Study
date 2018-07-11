#ifndef fitter_h
#define fitter_h

#include "TGraph.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH1.h"
#include "TF1.h"
#include <vector>

const int MAXBD  = 28;
const int MAXSKI = 4;
const int MAXCH  = 32;

using namespace std;
class fitter{
 public:
  fitter();
  fitter(TGraph *in_gr,vector<double> HG,vector<double> LG,vector<double> TOT);
 ~fitter();
  
 void fit_Graph();
 void fit(int labelE);
 
 bool   status;
 double p0;
 double p1;
 double sat_point;
 double undershoot_percent;

 
 private:
 
 TGraph *gr;
 TCanvas *c1;
 void fit_Draw();
 void ratio_plot(TProfile *tpr,TF1 *fit,TH1D *hratio);
 void root_logon();

 
 void Find_low(TH1D* h1,double* lowx,double* lowy);
 void Find_high(TH1D* h1,double* highx,double* highy);
 bool Find_sat(TProfile *tpr, TH1D* h1, double *sat,double *sat_x,double thres);
 double Calc_avg(TH1D* h1 , double min,double max);
 TH1D* qualify(double ARR[MAXBD][MAXSKI][MAXCH],int option);

 
 TPad *pad1,*pad2;
 vector<bool>   fit_remove;
 vector<double> HG_vec;
 vector<double> LG_vec;
 vector<double> TOT_vec;
 
};
#endif
