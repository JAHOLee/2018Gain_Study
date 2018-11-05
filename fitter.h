#ifndef fitter_h
#define fitter_h

#include "TGraph.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TSpline.h"
#include "TH1.h"
#include "TF1.h"
#include <vector>

const int MAXBD  = 1;
const int MAXSKI = 1;
const int MAXCH  = 2;

using namespace std;

class output{
public:
  int L_ID = -1;
  int M_ID = -1;
  int S_ID = -1;
  int C_ID = -1;
  double A2M  = -1;
  double L2HT = -1;
  double L2H  = -1;
  double T2L  = -1;
  double T2LT = -1;
  double TOFF = -1;
  bool   HLTYPE = false;
  bool   LTTYPE = false;
};


class fitter{
 public:
  fitter();
  fitter(TGraph *in_gr,vector<double> HG,vector<double> LG,vector<double> TOT);
 ~fitter();
  
 void fit_Graph();
 void fit(int labelE = 100);
 void fit_spline(int labelE = 100);
 void look_detail();
 void fit_LGTOT(int labelE = 100);
 void root_logon();
 void fit_output(int labelE = 100);
 void DEBUG(); 


 
 bool   status;
 double p0;
 double p1;
 double sat_point;
 double undershoot_percent;

 string fname;
 
 private:
 
 TGraph *gr;
 TCanvas *c1;
 output opt_val[28*4*32];
 void fit_Draw();
 void ratio_plot(TProfile *tpr,TF1 *fit,TH1D *hratio,string X_title = "LG",string Y_title = "HG");
 void Draw_Spline_and_1stderi(TProfile& tpr, TSpline3 &s, TH1D& h_deri);

 void Find_low(TH1D* h1,double* lowx,double* lowy);
 void Find_high(TH1D* h1,double* highx,double* highy, double lowerbound = 200,
 double upperbound = 300);
 bool Find_sat(TProfile *tpr, TH1D* h1, double *sat,double *sat_x,double thres);
 double Calc_avg(TH1D* h1 , double min,double max);
 TH1D* qualify(double ARR[MAXBD][MAXSKI][MAXCH],int option);
 double spline_4nodes(double *x, double *par);
 void pol4(TProfile *tpr);
 void TSpline_2nd_deri(TH1D& h_2nd_deri,TH1D *h_deri,double* sat_x);
 double first_P_lower_thres(TH1D* hist,double thres,double lowerbond,double upperbond);
 vector<double> findzeros(TH1D* hist);
 
 
 TPad *pad1,*pad2;

 vector<bool>   fit_remove;
 vector<double> HG_vec;
 vector<double> LG_vec;
 vector<double> TOT_vec;
 
};


#endif
