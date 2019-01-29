#ifndef fitter_h
#define fitter_h

#include "setup_config.h"
#include "MakePlots.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TSpline.h"
#include "TLine.h"
#include "TH1.h"
#include "TF1.h"
#include <vector>

const int MAXBD  = 94;

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
    LGTOT_FitSKI = 0;
    TOT_THRES_LG = -1;
    THRES_TYPE   = 0;
  };
  
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
  int    LGTOT_FitSKI;
  double TOT_THRES_LG;
  bool   THRES_TYPE;

};


class fitter{
 public:
  fitter(setup_config *SC,MakePlots *M, string filename);
 ~fitter();
  
 void fit_Graph();
 void fit(int labelE = 100);
 void fit_spline();
 void look_detail();
 void fit_LGTOT();
 void fit_output();
 void DEBUG(); 

 
 bool   status;
 double p0;
 double p1;
 double sat_point;
 double undershoot_percent;

 
 private:

 //Option to disable canvas
 bool batch_mode;

 // members
 TGraph *gr;
 TCanvas *c1;
 output opt_val[MAXBOARDS*MAXSKI*MAXCH];
 void fit_Draw();
 void ratio_plot(TProfile *tpr,TF1 *fit,TH1D *hratio,string X_title = "LG",string Y_title = "HG");
 void Draw_Spline_and_1stderi(TProfile& tpr, TSpline3 &s, TH1D& h_deri);

 void Find_low(TH1D* h1,double* lowx,double* lowy);
 void Find_high(TH1D* h1,double* highx,double* highy, double lowerbound = 200,
		double upperbound = 300,bool prevent_below_zero = false);
 bool Find_sat(TProfile *tpr, TH1D* h1, double *sat,double *sat_x,double thres);
 double Calc_avg(TH1D* h1 , double min,double max);
 TH1D* qualify(double ARR[MAXBOARDS][MAXSKI][MAXCH],int option);
 double spline_4nodes(double *x, double *par);
 void pol4(TProfile *tpr);
 void TSpline_2nd_deri(TH1D& h_2nd_deri,TH1D *h_deri,double* sat_x);
 double first_P_lower_thres(TH1D* hist,double thres,double lowerbond,double upperbond);
 vector<double> findzeros(TH1D* hist);
 void LGTOT_result(TProfile& tpr, TSpline3 &s, TH1D& h_deri,TLine& left, TLine& right, TF1& Fit);
 

 vector<bool>   fit_remove;
 vector<double> HG_vec;
 vector<double> LG_vec;
 vector<double> TOT_vec;
 setup_config *mysetup;
 MakePlots    *myplots;
 string fname;
 
};


#endif
