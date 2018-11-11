#include "single_module.h"
#include "TBReader.h"
#include <fstream>
#include <iostream>
#include "TCanvas.h"

int main(){
  TApplication *app = new TApplication("app",0,0);
  //fitter fit;
  //fit.fit_LGTOT(-1);
  //fit.fit_output(-1);
  //fit.fit_output(100);

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
  //return 0;
  
  //TChain *chain_single = new TChain("treeproducer/sk2cms");

  string filename;
  //string dirpath = "/afs/cern.ch/user/c/chchien/HG_LG/";
  string dirpath = "./";
  
  // Initialize output dirrectory
  TBReader init_dir;
  init_dir.dirpath = dirpath;
  init_dir.Make_dir();
  
  
  int TB_member     = 0;
  int single_member = 0;
  
  string inputfile = "./data_input.txt";
  ifstream infile(inputfile.c_str());

  string TProfileoutname = "TPro_test.root";
  MakePlots *M = new MakePlots;
  M->Init_TFile(TProfileoutname);
  //TCanvas *c1 = new TCanvas();
  
  while(true){
    
    infile >> filename;
    if(infile.eof()) {
      M-> Write_TProfile();
      break;}
    if( filename.length() > 2){
      cout << "input file: " << filename << endl;

      TFile f( filename.c_str() );
      //check if root directories exist...
      bool single_tree,pulseshape_tree,TB_ntuple,trackimpactntupler;
      single_tree = f.GetListOfKeys()->Contains("treeproducer");
      pulseshape_tree = f.GetListOfKeys()->Contains("pulseshapeplotter");
      TB_ntuple   = f.GetListOfKeys()->Contains("rechitntupler");
      trackimpactntupler = f.GetListOfKeys()->Contains("trackimpactntupler");
      
      if(single_tree && pulseshape_tree){
	single_member++;
	TChain *chain_single  = new TChain("pulseshapeplotter/tree");
        chain_single->Add(filename.c_str());
	single_module S(chain_single,filename);
	S.Loop();
	delete chain_single;
      }
      else if(TB_ntuple && trackimpactntupler){
	TB_member++;
	TChain *chain  = new TChain("rechitntupler/hits");
	TChain *chain2 = new TChain("trackimpactntupler/impactPoints");
	chain ->Add(filename.c_str());
	chain2->Add(filename.c_str());
	TBReader TBReader(chain,chain2,filename);
	TBReader.dirpath = dirpath;
	//TBReader.Ntuple_Maker();
	//TBReader.TProfile_Maker();
	TBReader.TProfile_Maker(M);
	if(TB_member %5 == 0){
	  M-> Write_TProfile();}
	delete chain;
	delete chain2;
      }
      else{
	cout << filename.c_str() << " contains unknown tree to me ..." << endl;
      }
      f.Close();
    }
    else{
      cout << "file " << filename << " is not available, please check "
	   << inputfile << endl;}
  }
  infile.close();
  
  return(0);
}

