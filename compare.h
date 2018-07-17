#ifndef compare_h
#define compare_h

#include "TFile.h"
#include <string>
#include "fitter.h"
#include "TCanvas.h"
#include "TMultiGraph.h"

using namespace std;
const int MAXLABEL = 12;

class compare : public fitter{
 public:
  compare();
  ~compare();
  void compare_Ene(int method = 1);
  void compare_method(int Ene = 100);

 private:
  
  void store_GR(string fname,int label);
  void Set_mgr(TMultiGraph& mgr,string xtitle,string ytitle);
  
  TGraph *gr_p0[MAXLABEL][MAXBD];
  TGraph *gr_p1[MAXLABEL][MAXBD];
  TGraph *gr_sat[MAXLABEL][MAXBD];

  char title[100];
  TCanvas *c1;
  
};

#endif
