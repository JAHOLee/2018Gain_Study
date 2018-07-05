#include "makePlots.h"
#include "fitter.h"
#include <fstream>
#include <iostream>

int main(){
  TApplication *app = new TApplication("app",0,0);
  // fitter fit;
  // fit.fit(10);
  // return 0;

  TChain *chain = new TChain("pulseshapeplotter/tree");
  
  string filename;
  ifstream infile("input.txt");
  while(true){
    
    infile >> filename;
    if(infile.eof()) break;
    if( filename.length() > 2){
      cout << "input file: " << filename << endl;
      chain->Add(filename.c_str());
    }
    else{
      cout << "file " << filename << " is not available, please check input.txt" << endl;}
  }
  infile.close();

  makePlots M(chain,filename);
  M.Loop();
  return(0);
}
