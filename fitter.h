#ifndef fitter_h
#define fitter_h

#include "TGraph.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TSpline.h"
#include "TH1.h"
#include "TF1.h"
#include <vector>

const int MAXBD  = 28;
const int MAXSKI = 4;
const int MAXCH  = 32;

using namespace std;

class output{
 public:
  output(){ Init(); };
  void Init(){
    L_ID = -1;
    M_ID = -1;
    S_ID = -1;
    C_ID = -1;
    A2M  = -1;
    L2HT = -1;
    L2H  = -1;
    T2L  = -1;
    T2LT = -1;
    TOFF = -1;
    HLTYPE = false;
    LTTYPE = false;
    HGLG_FitSKI  = 0;
    LGTOT_FitSKI = 0;}
  int L_ID ;
  int M_ID ;
  int S_ID ;
  int C_ID ;
  double A2M  ;
  double L2HT ;
  double L2H  ;
  double T2L  ;
  double T2LT ;
  double TOFF ;
  bool   HLTYPE;
  bool   LTTYPE;
  int    HGLG_FitSKI;
  int LGTOT_FitSKI;


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