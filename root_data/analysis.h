//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  5 09:50:22 2018 by ROOT version 6.10/08
// from TTree tree/Pulse shape fitter results
// found on file: HexaOutput_91_noCalib.root
//////////////////////////////////////////////////////////

#ifndef analysis_h
#define analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventID;
   vector<int>     *skirocID;
   vector<int>     *boardID;
   vector<int>     *moduleID;
   vector<int>     *channelID;
   vector<float>   *HighGainADC;
   vector<float>   *HighGainTmax;
   vector<float>   *HighGainChi2;
   vector<float>   *HighGainErrorADC;
   vector<float>   *HighGainErrorTmax;
   vector<int>     *HighGainStatus;
   vector<int>     *HighGainNCalls;
   vector<float>   *LowGainADC;
   vector<float>   *LowGainTmax;
   vector<float>   *LowGainChi2;
   vector<float>   *LowGainErrorADC;
   vector<float>   *LowGainErrorTmax;
   vector<int>     *LowGainStatus;
   vector<int>     *LowGainNCalls;
   vector<int>     *TotSlow;
   vector<int>     *ToaRise;
   vector<int>     *ToaFall;

   // List of branches
   TBranch        *b_eventID;   //!
   TBranch        *b_skirocID;   //!
   TBranch        *b_boardID;   //!
   TBranch        *b_moduleID;   //!
   TBranch        *b_channelID;   //!
   TBranch        *b_HighGainADC;   //!
   TBranch        *b_HighGainTmax;   //!
   TBranch        *b_HighGainChi2;   //!
   TBranch        *b_HighGainErrorADC;   //!
   TBranch        *b_HighGainErrorTmax;   //!
   TBranch        *b_HighGainStatus;   //!
   TBranch        *b_HighGainNCalls;   //!
   TBranch        *b_LowGainADC;   //!
   TBranch        *b_LowGainTmax;   //!
   TBranch        *b_LowGainChi2;   //!
   TBranch        *b_LowGainErrorADC;   //!
   TBranch        *b_LowGainErrorTmax;   //!
   TBranch        *b_LowGainStatus;   //!
   TBranch        *b_LowGainNCalls;   //!
   TBranch        *b_TotSlow;   //!
   TBranch        *b_ToaRise;   //!
   TBranch        *b_ToaFall;   //!

   analysis(TTree *tree=0);
   virtual ~analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef analysis_cxx
analysis::analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("HexaOutput_91_noCalib.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("HexaOutput_91_noCalib.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("HexaOutput_91_noCalib.root:/pulseshapeplotter");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

analysis::~analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t analysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   skirocID = 0;
   boardID = 0;
   moduleID = 0;
   channelID = 0;
   HighGainADC = 0;
   HighGainTmax = 0;
   HighGainChi2 = 0;
   HighGainErrorADC = 0;
   HighGainErrorTmax = 0;
   HighGainStatus = 0;
   HighGainNCalls = 0;
   LowGainADC = 0;
   LowGainTmax = 0;
   LowGainChi2 = 0;
   LowGainErrorADC = 0;
   LowGainErrorTmax = 0;
   LowGainStatus = 0;
   LowGainNCalls = 0;
   TotSlow = 0;
   ToaRise = 0;
   ToaFall = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("skirocID", &skirocID, &b_skirocID);
   fChain->SetBranchAddress("boardID", &boardID, &b_boardID);
   fChain->SetBranchAddress("moduleID", &moduleID, &b_moduleID);
   fChain->SetBranchAddress("channelID", &channelID, &b_channelID);
   fChain->SetBranchAddress("HighGainADC", &HighGainADC, &b_HighGainADC);
   fChain->SetBranchAddress("HighGainTmax", &HighGainTmax, &b_HighGainTmax);
   fChain->SetBranchAddress("HighGainChi2", &HighGainChi2, &b_HighGainChi2);
   fChain->SetBranchAddress("HighGainErrorADC", &HighGainErrorADC, &b_HighGainErrorADC);
   fChain->SetBranchAddress("HighGainErrorTmax", &HighGainErrorTmax, &b_HighGainErrorTmax);
   fChain->SetBranchAddress("HighGainStatus", &HighGainStatus, &b_HighGainStatus);
   fChain->SetBranchAddress("HighGainNCalls", &HighGainNCalls, &b_HighGainNCalls);
   fChain->SetBranchAddress("LowGainADC", &LowGainADC, &b_LowGainADC);
   fChain->SetBranchAddress("LowGainTmax", &LowGainTmax, &b_LowGainTmax);
   fChain->SetBranchAddress("LowGainChi2", &LowGainChi2, &b_LowGainChi2);
   fChain->SetBranchAddress("LowGainErrorADC", &LowGainErrorADC, &b_LowGainErrorADC);
   fChain->SetBranchAddress("LowGainErrorTmax", &LowGainErrorTmax, &b_LowGainErrorTmax);
   fChain->SetBranchAddress("LowGainStatus", &LowGainStatus, &b_LowGainStatus);
   fChain->SetBranchAddress("LowGainNCalls", &LowGainNCalls, &b_LowGainNCalls);
   fChain->SetBranchAddress("TotSlow", &TotSlow, &b_TotSlow);
   fChain->SetBranchAddress("ToaRise", &ToaRise, &b_ToaRise);
   fChain->SetBranchAddress("ToaFall", &ToaFall, &b_ToaFall);
   Notify();
}

Bool_t analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analysis_cxx
