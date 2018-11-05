#include "makePlots.h"
#include "fitter.h"
#include "compare.h"
#include "single_module.h"
#include <fstream>
#include <iostream>

int main(){
  TApplication *app = new TApplication("app",0,0);
  fitter fit;
  //fit.fit_LGTOT(-1);
  //fit.fit_output(-1);
  //fit.fit_output(100);
  fit.fit_spline(-1);
  fit.fit_LGTOT(-1);
  fit.fit_LGTOT(100);
  fit.fit_spline(100);

  // fit.DEBUG();
  // cout <<"2" << endl;
  // fit.DEBUG();
  // int Ene_arr[6] = {10,30,50,80,100,150};
  // for(int i = 0 ;i < (int)sizeof(Ene_arr)/sizeof(int) ; ++i){
  //   fit.fit(Ene_arr[i]);
  //   fit.fit_spline(Ene_arr[i]);
  // }

  //compare com;
  // com.compare_Ene();
  //com.compare_method();
  return 0;
  
  TChain *chain  = new TChain("rechitntupler/hits");
  TChain *chain2 = new TChain("trackimpactntupler/impactPoints");

  //TChain *chain_single = new TChain("treeproducer/sk2cms");

  int TB_member     = 0;
  int single_member = 0;
  
  string filename;
  //string dirpath = "/afs/cern.ch/user/c/chchien/HG_LG/";
  string dirpath = "./";

  ifstream infile("input.txt");
  
  while(true){
    
    infile >> filename;
    if(infile.eof()) break;
    if( filename.length() > 2){
      cout << "input file: " << filename << endl;

      TFile f( filename.c_str() );
      //check if directories exist...
      bool single_tree,pulseshape_tree,TB_ntuple;
      single_tree = f.GetListOfKeys()->Contains("treeproducer");
      pulseshape_tree = f.GetListOfKeys()->Contains("pulseshapeplotter");
      TB_ntuple   = f.GetListOfKeys()->Contains("rechitntupler");

      if(single_tree && pulseshape_tree){
	single_member++;
	TChain *chain_single  = new TChain("pulseshapeplotter/tree");
        chain_single->Add(filename.c_str());
	single_module S(chain_single,filename);
	S.Loop();
	delete chain_single;
      }
      else if(TB_ntuple){
	TB_member++;
	chain ->Add(filename.c_str());
	chain2->Add(filename.c_str());}
      else{
	cout << filename.c_str() << " contains unknown tree to me ..." << endl;
      }
      f.Close();
    }
    else{
      cout << "file " << filename << " is not available, please check input.txt" << endl;}
  }
  infile.close();
  
  if(TB_member != 0){
    makePlots M(chain,chain2,filename);
    M.dirpath = dirpath;
    M.Loop();}
  //single_module S(chain,filename);
  //S.Loop();
  return(0);
}
