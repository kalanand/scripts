// $Id: setupChains.cc,v 1.32 2006/06/28 02:00:15 fwinkl Exp $
// setup the global files

#include "TTree.h"
#include <iostream>
#include "../BToDKTo3piK/globals/chains.hh"

// forward
Bool_t checkMerge(TTree* mergeTree, TTree* t1, TTree* t2);


// Close files
void cleanupChains()
{
   // Close all files
  for (unsigned int i=0; i<sizeof(fileList)/sizeof(TFile**); i++) {
    TFile* f = *fileList[i];
    if (f) f->Close();
  }
}


void setupChains(Bool_t bestMes = kTRUE)
{
  cout << "--- START setupChains() ---" << std::endl;

  cleanupChains();

  if (bestMes) cout << "Using best mES selection ntuples."<<endl;
  else cout << "Using random selection ntuples."<<endl;

  //  string path = "/nfs/farm/babar/AWG17/BCK/SmallN/";
  //string path = "/nfs/farm/babar/AWG17/BCK/Frank/ntuples/R16R18-MC/";
  //  string r16path = "/nfs/farm/babar/AWGbreco02/BDK/R16Small/";
  //  string r18path = "/nfs/farm/babar/AWGbreco02/BDK/R18Small/";

  string path, r16path, r18path, r18DataPath, r18OffPath, ktPath;
  string suffix;

  if (bestMes) {
    path = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/SPmerge_F/";
    r16path = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/SP56_F/";
    r18path = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/SP8_F/";
    r18DataPath = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/RUN15_F/";
    r18OffPath = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/Off/";
    ktPath = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/KT/";
    suffix = "f_wt2-";
  }
  else {
    path = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/SPmerge_R/";
    r16path = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/SP56_R/";
    r18path = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/SP8_R/";
    r18DataPath = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/RUN15_R/";
    r18OffPath = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/Off/";
    ktPath = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/KT/";
    suffix = "rad_wt2-";
  }

  // These are the files used for the analysis
  // They are merged into one big root file since chains caused problems


  // phase space signal MC
  sigFlatFile = new TFile((r18path+"4854_wt2.root").c_str());
  sigFlatTree = (TTree*)sigFlatFile->Get("h1");
  sigFlatTree->SetName("sigFlatTree");
  sigFlatTree->SetTitle("flat signal");

  // phase space signal MC with PID weights
  sigFlatWFile = new TFile((ktPath+"4854killn_wt2.root").c_str());
  sigFlatWTree = (TTree*)sigFlatWFile->Get("h1");
  sigFlatWTree->SetName("sigFlatWTree");
  sigFlatWTree->SetTitle("flat signal with PID weights");

  /*
  // phase space signal MC after PID tweaking
  sigTweakFile = new TFile((ktPath+"tweak4854_wt2.root").c_str());
  sigTweakTree = (TTree*)sigTweakFile->Get("h1");
  sigTweakTree->SetName("sigTweakTree");
  sigTweakTree->SetTitle("flat signal after PID tweaking");

  // phase space signal MC after PID killing
  sigKillFile = new TFile((ktPath+"kill4854_wt2.root").c_str());
  sigKillTree = (TTree*)sigKillFile->Get("h1");
  sigKillTree->SetName("sigKillTree");
  sigKillTree->SetTitle("flat signal after PID killing");
  */


  // signal MC
  sigFile = new TFile((r18path+"6795_wt2.root").c_str());
  sigTree = (TTree*)sigFile->Get("h1");
  sigTree->SetName("sigTree");
  sigTree->SetTitle("signal");

  // signal MC with PID weights
  sigWFile = new TFile((ktPath+"6795killn_wt2.root").c_str());
  sigWTree = (TTree*)sigWFile->Get("h1");
  sigWTree->SetName("sigWTree");
  sigWTree->SetTitle("signal with PID weights");

  // ccbar + uds MC
  qqFile = new TFile((path+"qqbar.root").c_str());
  qqTree = (TTree*)qqFile->Get("h1");
  qqTree->SetName("qqTree");
  qqTree->SetTitle("qqbar");

  // B+B- + B0B0bar MC
  bbFile = new TFile((path+"bbbar.root").c_str());
  bbTree = (TTree*)bbFile->Get("h1");
  bbTree->SetName("bbTree");
  bbTree->SetTitle("B+/B0");

  // DPi events
  dpiFile = new TFile((path+"1235dpi_wt2.root").c_str());
  dpiTree = (TTree*)dpiFile->Get("h1");
  dpiTree->SetName("dpiTree");
  dpiTree->SetTitle("DPi");

  // on peak data
  dataFile = new TFile((r18DataPath+"onrun.root").c_str());
  dataTree = (TTree*)dataFile->Get("h1");
  dataTree->SetName("dataTree");
  dataTree->SetTitle("onPeak");  

  // off peak data
  offFile = new TFile((r18OffPath+"offrun.root").c_str());
  offTree = (TTree*)offFile->Get("h1");
  offTree->SetName("offTree");
  offTree->SetTitle("offPeak");  

  // These are the individual root files in chains for simple studies
  
  // B+B- MC
  bpTree = new TChain("h1","B+B-");
  bpTree->Add((r16path+"1235other"+suffix+"*.root").c_str());
  bpTree->Add((r18path+"1235other"+suffix+"*.root").c_str());
  
  // B0B0bar MC
  b0Tree = new TChain("h1","B0B0bar");
  b0Tree->Add((r16path+"1237other"+suffix+"*.root").c_str());
  b0Tree->Add((r18path+"1237other"+suffix+"*.root").c_str());
    
  // uds MC
  udsTree = new TChain("h1","uds");
  udsTree->Add((r16path+"998other"+suffix+"*.root").c_str());
  udsTree->Add((r18path+"998other"+suffix+"*.root").c_str());
    
  // ccbar MC
  ccTree = new TChain("h1","ccbar");
  ccTree->Add((r16path+"1005other"+suffix+"*.root").c_str());
  ccTree->Add((r18path+"1005other"+suffix+"*.root").c_str());
  

  //
  // Release 16 ntuples
  //
  
  // flat signal MC
  sigFlatFile16 = new TFile((r16path+"4854_wt2.root").c_str());
  sigFlatTree16 = (TTree*)sigFlatFile16->Get("h1");
  sigFlatTree16->SetName("sigFlatTree16");
  sigFlatTree16->SetTitle("flat signal R16");
  
   // signal MC
  sigFile16 = new TFile((r16path+"6795_wt2.root").c_str());
  sigTree16 = (TTree*)sigFile16->Get("h1");
  sigTree16->SetName("sigTree16");
  sigTree16->SetTitle("signal R16");

  // ccbar + uds MC
  qqFile16 = new TFile((r16path+"qqbar.root").c_str());
  qqTree16 = (TTree*)qqFile16->Get("h1");
  qqTree16->SetName("qqTree16");
  qqTree16->SetTitle("qqbar R16");

  // B+B- + B0B0bar MC
  bbFile16 = new TFile((r16path+"bbbar.root").c_str());
  bbTree16 = (TTree*)bbFile16->Get("h1");
  bbTree16->SetName("bbTree16");
  bbTree16->SetTitle("B+/B0 R16");

  // DPi events
  dpiFile16 = new TFile((r16path+"1235dpi_wt2.root").c_str());
  dpiTree16 = (TTree*)dpiFile16->Get("h1");
  dpiTree16->SetName("dpiTree16");
  dpiTree16->SetTitle("DPi R16");

  // on peak data
  /*
  dataFile16 = new TFile((r16path+"onrun.root").c_str());
  dataTree16 = (TTree*)dataFile16->Get("h1");
  dataTree16->SetName("dataTree16");
  dataTree16->SetTitle("onPeak");
  */
  

  // B+B- MC R16
  bpTree16 = new TChain("h1","B+B- R16");
  bpTree16->Add((r16path+"1235other"+suffix+"*.root").c_str());
  
  // B0B0bar MC R16
  b0Tree16 = new TChain("h1","B0B0bar R16");
  b0Tree16->Add((r16path+"1237other"+suffix+"*.root").c_str());
    
  // uds MC R16
  udsTree16 = new TChain("h1","uds R16");
  udsTree16->Add((r16path+"998other"+suffix+"*.root").c_str());
    
  // ccbar MC R16
  ccTree16 = new TChain("h1","ccbar R16");
  ccTree16->Add((r16path+"1005other"+suffix+"*.root").c_str());

  
    
  //
  // Release 18 ntuples
  //

  // flat signal MC
  sigFlatFile18 = sigFlatFile;
  sigFlatTree18 = sigFlatTree;

  // signal MC
  sigFile18 = sigFile;
  sigTree18 = sigTree;

  // ccbar + uds MC
  qqFile18 = new TFile((r18path+"qqbar.root").c_str());
  qqTree18 = (TTree*)qqFile18->Get("h1");
  qqTree18->SetName("qqTree18");
  qqTree18->SetTitle("qqbar R18");

  // B+B- + B0B0bar MC
  bbFile18 = new TFile((r18path+"bbbar.root").c_str());
  bbTree18 = (TTree*)bbFile18->Get("h1");
  bbTree18->SetName("bbTree18");
  bbTree18->SetTitle("B+/B0 R18");

  // DPi events
  dpiFile18 = new TFile((r18path+"1235dpi_wt2.root").c_str());
  dpiTree18 = (TTree*)dpiFile18->Get("h1");
  dpiTree18->SetName("dpiTree18");
  dpiTree18->SetTitle("DPi R18");

  // on peak data
  dataFile18 = dataFile;
  dataTree18 = dataTree;

  // B+B- MC R18
  bpTree18 = new TChain("h1","B+B- R18");
  bpTree18->Add((r18path+"1235other"+suffix+"*.root").c_str());
  
  // B0B0bar MC R18
  b0Tree18 = new TChain("h1","B0B0bar R18");
  b0Tree18->Add((r18path+"1237other"+suffix+"*.root").c_str());
    
  // uds MC R18
  udsTree18 = new TChain("h1","uds R18");
  udsTree18->Add((r18path+"998other"+suffix+"*.root").c_str());
    
  // ccbar MC R18
  ccTree18 = new TChain("h1","ccbar R18");
  ccTree18->Add((r18path+"1005other"+suffix+"*.root").c_str());



  // Make a list of all the data trees
  // This is used to allow unblinded fits on MC

  dataTreeList = new TList();
  if (dataTree) dataTreeList->Add(dataTree);
  if (dataTree18) dataTreeList->Add(dataTree18);

  std::cout << "--- END setupChains() ---" << std::endl;

}


// Do some consistency checks on the trees
Bool_t checkTrees()
{
  Bool_t OK = kTRUE;
  
  OK &= checkMerge(bbTree,bpTree,b0Tree);
  OK &= checkMerge(qqTree,ccTree,udsTree);
  //  OK &= checkMerge(sigFlatTree,sigFlatTree18,sigFlatTree16);
  //  OK &= checkMerge(sigTree,sigTree18,sigTree16);
  OK &= checkMerge(bbTree,bbTree18,bbTree16);
  OK &= checkMerge(qqTree,qqTree18,qqTree16);
  OK &= checkMerge(dpiTree,dpiTree18,dpiTree16);
  OK &= checkMerge(qqTree16,ccTree16,udsTree16);
  OK &= checkMerge(qqTree18,ccTree18,udsTree18);
  OK &= checkMerge(bbTree16,bpTree16,b0Tree16);
  OK &= checkMerge(bbTree18,bpTree18,b0Tree18);

  return OK;
}

// Check the merge of two trees
Bool_t checkMerge(TTree* mergeTree, TTree* t1, TTree* t2)
{
  if (mergeTree==0 || t1==0 || t2==0) {
    std::cout << "chechMerge: One of the trees is missing."<<std::endl;
    return kFALSE;
  }

  Long64_t N1 = t1->GetEntries();
  Long64_t N2 = t2->GetEntries();
  Long64_t N = mergeTree->GetEntries();

  if (N1+N2 != N) {
    std::cout << "checkMerge: ERROR: Tree "<<mergeTree->GetName()<<" ("<<N<<" entries)"
         << " is not a merge of:"<<std::endl
         << "Tree "<<t1->GetName()<<" ("<<N1<<" entries)"<<std::endl
         << "Tree "<<t2->GetName()<<" ("<<N2<<" entries)"<<std::endl;
    return kFALSE;
  }
  return kTRUE;
}
