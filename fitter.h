#ifndef fitter_h
#define fitter_h

#include "TGraph.h"
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
 void fit_Draw();
 vector<bool>   fit_remove;
 vector<double> HG_vec;
 vector<double> LG_vec;
 vector<double> TOT_vec;
 
};
#endif
