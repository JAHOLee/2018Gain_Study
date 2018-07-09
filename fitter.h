#ifndef fitter_h
#define fitter_h

#include "TGraph.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH1.h"
#include "TF1.h"
#include <vector>

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
 void root_logon();
 void ratio_plot(TProfile *tpr,TF1 *fit,TH1D *hratio);
 TPad *pad1,*pad2;
 vector<bool>   fit_remove;
 vector<double> HG_vec;
 vector<double> LG_vec;
 vector<double> TOT_vec;
 
};
#endif
