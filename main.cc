#include "makePlots.h"
#include "fitter.h"
#include <fstream>
#include <iostream>

int main(){
  TApplication *app = new TApplication("app",0,0);
  // fitter fit;
  // int Ene_arr[6] = {10,30,50,80,100,150};
  // for(int i = 0 ;i < (int)sizeof(Ene_arr)/sizeof(int) ; ++i){
  //   fit.fit(Ene_arr[i]);}
  // return 0;

  TChain *chain = new TChain("pulseshapeplotter/tree");
  
  string filename;
  //string dirpath = "/afs/cern.ch/user/c/chchien/HG_LG/";
  string dirpath = "./";

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
  M.dirpath = dirpath;
  M.Loop();
  return(0);
}
