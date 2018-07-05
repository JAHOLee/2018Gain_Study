#include "TH1.h"
#include "TH2.h"

void makeClass(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  //TApplication *a = new TApplication("a", 0, 0);
  bool data = 0;
  TFile f("HexaOutput_91_noCalib.root");
  TTree *tt = (TTree *) f.Get("pulseshapeplotter/tree");
  tt->MakeClass("analysis");
}
