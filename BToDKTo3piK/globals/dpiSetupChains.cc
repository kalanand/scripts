// $Id: dpiSetupChains.cc,v 1.1 2006/05/29 22:56:57 fwinkl Exp $
// setup the global DPI files

#include "TTree.h"
#include <iostream>
#include "../BToDKTo3piK/globals/chains.hh"


void dpiSetupChains()
{
  cout << "--- START dpiSetupChains() ---" << std::endl;

  string dpiPath = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/dpi_F/";

   //
  // B -> D pi ntuples
  //

  dpiDataFile = new TFile((dpiPath+"data/onrun.root").c_str());
  dpiDataTree = (TTree*)dpiDataFile->Get("h1");
  dpiDataTree->SetName("dpiDataTree");
  dpiDataTree->SetTitle("DPI data");

  dpiBbFile = new TFile((dpiPath+"SPmerge/bbbar.root").c_str());
  dpiBbTree = (TTree*)dpiBbFile->Get("h1");
  dpiBbTree->SetName("dpiBbTree");
  dpiBbTree->SetTitle("DPI B+/B0");

  dpiQqFile = new TFile((dpiPath+"SPmerge/qqbar.root").c_str());
  dpiQqTree = (TTree*)dpiQqFile->Get("h1");
  dpiQqTree->SetName("dpiQqTree");
  dpiQqTree->SetTitle("DPI qqbar");

  dpiSigFile = new TFile((dpiPath+"SPmerge/1235dpi_wt2.root").c_str());
  dpiSigTree = (TTree*)dpiSigFile->Get("h1");
  dpiSigTree->SetName("dpiSigTree");
  dpiSigTree->SetTitle("DPI signal");

  dpiB0Tree = new TChain("h1","DPI B0B0bar");
  dpiB0Tree->Add((dpiPath+"/R16/1237dpiother*.root").c_str());
  dpiB0Tree->Add((dpiPath+"/R18/1237dpiother*.root").c_str());

  dpiBpTree = new TChain("h1","DPI B+B-");
  dpiBpTree->Add((dpiPath+"/R16/1235dpiother*.root").c_str());
  dpiBpTree->Add((dpiPath+"/R18/1235dpiother*.root").c_str());

  dpiUdsTree = new TChain("h1","DPI uds");
  dpiUdsTree->Add((dpiPath+"/R16/998dpiother*.root").c_str());
  dpiUdsTree->Add((dpiPath+"/R18/998dpiother*.root").c_str());

  dpiCcTree = new TChain("h1","DPI ccbar");
  dpiCcTree->Add((dpiPath+"/R16/1005dpiother*.root").c_str());
  dpiCcTree->Add((dpiPath+"/R18/1005dpiother*.root").c_str());
  

  // Release 16
  dpiBbFile16 = new TFile((dpiPath+"R16/bbbar.root").c_str());
  dpiBbTree16 = (TTree*)dpiBbFile16->Get("h1");
  dpiBbTree16->SetName("dpiBbTree16");
  dpiBbTree16->SetTitle("DPI B+/B0 R16");

  dpiQqFile16 = new TFile((dpiPath+"R16/qqbar.root").c_str());
  dpiQqTree16 = (TTree*)dpiQqFile16->Get("h1");
  dpiQqTree16->SetName("dpiQqTree16");
  dpiQqTree16->SetTitle("DPI qqbar R16");

  dpiSigFile16 = new TFile((dpiPath+"R16/1235dpi_wt2.root").c_str());
  dpiSigTree16 = (TTree*)dpiSigFile16->Get("h1");
  dpiSigTree16->SetName("dpiSigTree16");
  dpiSigTree16->SetTitle("DPI DPi R16");

  dpiB0Tree16 = new TChain("h1","DPI B0B0bar R16");
  dpiB0Tree16->Add((dpiPath+"/R16/1237dpiother*.root").c_str());

  dpiBpTree16 = new TChain("h1","DPI B+B- R16");
  dpiBpTree16->Add((dpiPath+"/R16/1235dpiother*.root").c_str());

  dpiUdsTree16 = new TChain("h1","DPI uds R16");
  dpiUdsTree16->Add((dpiPath+"/R16/998dpiother*.root").c_str());

  dpiCcTree16 = new TChain("h1","DPI ccbar R16");
  dpiCcTree16->Add((dpiPath+"/R16/1005dpiother*.root").c_str());


  // Release 18 
  dpiBbFile18 = new TFile((dpiPath+"R18/bbbar.root").c_str());
  dpiBbTree18 = (TTree*)dpiBbFile18->Get("h1");
  dpiBbTree18->SetName("dpiBbTree18");
  dpiBbTree18->SetTitle("DPI B+/B0 R18");

  dpiQqFile18 = new TFile((dpiPath+"R18/qqbar.root").c_str());
  dpiQqTree18 = (TTree*)dpiQqFile18->Get("h1");
  dpiQqTree18->SetName("dpiQqTree18");
  dpiQqTree18->SetTitle("DPI qqbar R18");
  
  dpiSigFile18 = new TFile((dpiPath+"R18/1235dpi_wt2.root").c_str());
  dpiSigTree18 = (TTree*)dpiSigFile18->Get("h1");
  dpiSigTree18->SetName("dpiSigTree18");
  dpiSigTree18->SetTitle("DPI DPi R18");

  dpiB0Tree18 = new TChain("h1","DPI B0B0bar R18");
  dpiB0Tree18->Add((dpiPath+"/R18/1237dpiother*.root").c_str());

  dpiBpTree18 = new TChain("h1","DPI B+B- R18");
  dpiBpTree18->Add((dpiPath+"/R18/1235dpiother*.root").c_str());

  dpiUdsTree18 = new TChain("h1","DPI uds R18");
  dpiUdsTree18->Add((dpiPath+"/R18/998dpiother*.root").c_str());

  dpiCcTree18 = new TChain("h1","DPI ccbar R18");
  dpiCcTree18->Add((dpiPath+"/R18/1005dpiother*.root").c_str());

  cout << "--- END dpiSetupChains() ---" << std::endl;
}
