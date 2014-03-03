//$Id: chains.hh,v 1.24 2006/06/05 21:06:19 fwinkl Exp $


// Define chains/trees and related quantities

#ifndef CHAINS_HH
#define CHAINS_HH

#include "TTree.h"
#include "TChain.h"
#include "TList.h"

// Luminosity and B-counting number from
// http://www.slac.stanford.edu/BFROOT/www/Computing/DataQuality/datasets/R18-Skims/R18Skims-v03-all.html

// On-peak data luminosity used to calcualte the weights
const Double_t ON_PEAK_LUMI = 288.48;

// Exact B-counting number (used to set pdfOnResDK.nBB() in setupPdfs)
const Double_t N_BB = 3.24041437e08;
const Double_t N_BB_ERR = 0.011;         // 1.1% systematic error

// MC luminosities
const Double_t BP_LUMI16  = 1059.94;
const Double_t B0_LUMI16  =  981.38;
const Double_t CC_LUMI16  =  308.73;
const Double_t UDS_LUMI16 =  324.03;
const Double_t DPI_LUMI16 = BP_LUMI16;

const Double_t BP_LUMI18  = 725.70;
const Double_t B0_LUMI18  =  644.18;
const Double_t CC_LUMI18  =  179.03;
const Double_t UDS_LUMI18 =  158.75;
const Double_t DPI_LUMI18 = BP_LUMI18;

const Double_t BP_LUMI  = BP_LUMI16 + BP_LUMI18;
const Double_t B0_LUMI  = B0_LUMI16 + B0_LUMI18;
const Double_t CC_LUMI  = CC_LUMI16 + CC_LUMI18;
const Double_t UDS_LUMI = UDS_LUMI16 + UDS_LUMI18;
const Double_t DPI_LUMI = BP_LUMI;

// MC sample weights
const Double_t BP_WEIGHT16 = ON_PEAK_LUMI/BP_LUMI16;
const Double_t B0_WEIGHT16 = ON_PEAK_LUMI/B0_LUMI16;
const Double_t CC_WEIGHT16 = ON_PEAK_LUMI/CC_LUMI16;
const Double_t UDS_WEIGHT16 = ON_PEAK_LUMI/UDS_LUMI16;
const Double_t DPI_WEIGHT16 = ON_PEAK_LUMI/DPI_LUMI16;

const Double_t BP_WEIGHT18 = ON_PEAK_LUMI/BP_LUMI18;
const Double_t B0_WEIGHT18 = ON_PEAK_LUMI/B0_LUMI18;
const Double_t CC_WEIGHT18 = ON_PEAK_LUMI/CC_LUMI18;
const Double_t UDS_WEIGHT18 = ON_PEAK_LUMI/UDS_LUMI18;
const Double_t DPI_WEIGHT18 = ON_PEAK_LUMI/DPI_LUMI18;

const Double_t BP_WEIGHT = ON_PEAK_LUMI/BP_LUMI;
const Double_t B0_WEIGHT = ON_PEAK_LUMI/B0_LUMI;
const Double_t CC_WEIGHT = ON_PEAK_LUMI/CC_LUMI;
const Double_t UDS_WEIGHT = ON_PEAK_LUMI/UDS_LUMI;
const Double_t DPI_WEIGHT = ON_PEAK_LUMI/DPI_LUMI;


// Signal MC sample weights
const Double_t SIG_BR = 5.5e-6;

const Int_t SIGFLAT_EVENTS16 = 347000;
const Int_t SIGFLAT_EVENTS18 = 542000;
const Int_t SIGFLAT_EVENTS = SIGFLAT_EVENTS18;   // we only use SP8 for signal

const Int_t SIG_EVENTS16 = 351000;
const Int_t SIG_EVENTS18 = 542000;
const Int_t SIG_EVENTS = SIG_EVENTS18;

const Double_t SIGFLAT_WEIGHT16 = N_BB/(SIGFLAT_EVENTS16/SIG_BR);
const Double_t SIGFLAT_WEIGHT18 = N_BB/(SIGFLAT_EVENTS18/SIG_BR);
const Double_t SIGFLAT_WEIGHT = N_BB/(SIGFLAT_EVENTS/SIG_BR);

const Double_t SIG_WEIGHT16 = N_BB/(SIG_EVENTS16/SIG_BR);
const Double_t SIG_WEIGHT18 = N_BB/(SIG_EVENTS18/SIG_BR);
const Double_t SIG_WEIGHT = N_BB/(SIG_EVENTS/SIG_BR);


TFile *sigFlatFile = 0;
TTree *sigFlatTree = 0;
TFile *sigTweakFile = 0;
TTree *sigTweakTree = 0;
TFile *sigKillFile = 0;
TTree *sigKillTree = 0;
TFile *sigFlatWFile = 0;
TTree *sigFlatWTree = 0;

TFile *sigFile = 0;
TTree *sigTree = 0;
TFile *sigWFile = 0;
TTree *sigWTree = 0;

TFile *qqFile = 0;
TTree *qqTree = 0;

TFile *bbFile = 0;
TTree *bbTree = 0;

TFile *dpiFile = 0;
TTree *dpiTree = 0;

TFile *dataFile = 0;
TTree *dataTree = 0;

TFile *offFile = 0;
TTree *offTree = 0;

TChain* bpTree = 0;
TChain* b0Tree = 0;
TChain* ccTree = 0;
TChain* udsTree = 0;

TChain* bpTree16 = 0;
TChain* b0Tree16 = 0;
TChain* ccTree16 = 0;
TChain* udsTree16 = 0;

TChain* bpTree18 = 0;
TChain* b0Tree18 = 0;
TChain* ccTree18 = 0;
TChain* udsTree18 = 0;

TFile* sigFlatFile16 = 0;
TFile* sigFile16 = 0;
TFile* qqFile16 = 0;
TFile* bbFile16 = 0;
TFile* dpiFile16 = 0;
//TFile* dataFile16 = 0;

TTree* sigFlatTree16 = 0;
TTree* sigTree16 = 0;
TTree* qqTree16 = 0;
TTree* bbTree16 = 0;
TTree* dpiTree16 = 0;
//TTree* dataTree16 = 0;

TFile* sigFlatFile18 = 0;
TFile* sigFile18 = 0;
TFile* qqFile18 = 0;
TFile* bbFile18 = 0;
TFile* dpiFile18 = 0;
TFile* dataFile18 = 0;

TTree* sigFlatTree18 = 0;
TTree* sigTree18 = 0;
TTree* qqTree18 = 0;
TTree* bbTree18 = 0;
TTree* dpiTree18 = 0;
TTree* dataTree18 = 0;

TList *dataTreeList = 0;

// List of all files, trees and chains

TFile** fileList[] = {&sigFlatFile, &sigFile, &qqFile, &bbFile, &dpiFile, &dataFile, &offFile,
                      &sigFlatFile16, &sigFile16, &qqFile16, &bbFile16, &dpiFile16,
                      &sigFlatFile18, &sigFile18, &qqFile18, &bbFile18, &dpiFile18, &dataFile18};

TChain** chainList[] = {&bpTree, &b0Tree, &ccTree, &udsTree, 
                        &bpTree16, &b0Tree16, &ccTree16, &udsTree16, 
                        &bpTree18, &b0Tree18, &ccTree18, &udsTree18};

TTree** treeList[] = {&sigFlatTree, &sigTree, &qqTree, &bbTree, &dpiTree, &dataTree, &offTree,
                      &sigFlatTree16, &sigTree16, &qqTree16, &bbTree16, &dpiTree16,
                      &sigFlatTree18, &sigTree18, &qqTree18, &bbTree18, &dpiTree18, &dataTree18};



/////////////////////////////////////////////////////////////////////////
//                          Dpi tress and files
/////////////////////////////////////////////////////////////////////////

// Files
TFile* dpiDataFile = 0;
TFile* dpiBbFile = 0;
TFile* dpiQqFile = 0;
TFile* dpiSigFile = 0;

TFile* dpiBbFile16 = 0;
TFile* dpiQqFile16 = 0;
TFile* dpiSigFile16 = 0;

TFile* dpiBbFile18 = 0;
TFile* dpiQqFile18 = 0;
TFile* dpiSigFile18 = 0;

// Trees
TTree* dpiDataTree = 0;
TTree* dpiBbTree = 0;
TTree* dpiQqTree = 0;
TTree* dpiSigTree = 0;

TTree* dpiBbTree16 = 0;
TTree* dpiQqTree16 = 0;
TTree* dpiSigTree16 = 0;

TTree* dpiBbTree18 = 0;
TTree* dpiQqTree18 = 0;
TTree* dpiSigTree18 = 0;


// Chains
TChain* dpiB0Tree = 0;
TChain* dpiBpTree = 0;
TChain* dpiUdsTree = 0;
TChain* dpiCcTree = 0;

TChain* dpiB0Tree16 = 0;
TChain* dpiBpTree16 = 0;
TChain* dpiUdsTree16 = 0;
TChain* dpiCcTree16 = 0;

TChain* dpiB0Tree18 = 0;
TChain* dpiBpTree18 = 0;
TChain* dpiUdsTree18 = 0;
TChain* dpiCcTree18 = 0;





Double_t getTreeWeight(TTree* tree)
{
  if (tree==bpTree) return BP_WEIGHT;
  else if (tree==b0Tree) return B0_WEIGHT;
  else if (tree==ccTree) return CC_WEIGHT;
  else if (tree==udsTree) return UDS_WEIGHT;
  else if (tree==dpiTree) return DPI_WEIGHT;
  else if (tree==sigTree) return SIG_WEIGHT;
  else if (tree==sigFlatTree) return SIGFLAT_WEIGHT;

  else if (tree==bpTree16) return BP_WEIGHT16;
  else if (tree==b0Tree16) return B0_WEIGHT16;
  else if (tree==ccTree16) return CC_WEIGHT16;
  else if (tree==udsTree16) return UDS_WEIGHT16;
  else if (tree==dpiTree16) return DPI_WEIGHT16;
  else if (tree==sigTree16) return SIG_WEIGHT16;
  else if (tree==sigFlatTree16) return SIGFLAT_WEIGHT16;

  else if (tree==bpTree18) return BP_WEIGHT18;
  else if (tree==b0Tree18) return B0_WEIGHT18;
  else if (tree==ccTree18) return CC_WEIGHT18;
  else if (tree==udsTree18) return UDS_WEIGHT18;
  else if (tree==dpiTree18) return DPI_WEIGHT18;
  else if (tree==sigTree18) return SIG_WEIGHT18;
  else if (tree==sigFlatTree18) return SIGFLAT_WEIGHT18;

  else return -1;
}


Double_t getTreeLumi(TTree* tree)
{
  if (tree==bpTree) return BP_LUMI;
  else if (tree==b0Tree) return B0_LUMI;
  else if (tree==ccTree) return CC_LUMI;
  else if (tree==udsTree) return UDS_LUMI;
  else if (tree==dpiTree) return DPI_LUMI;

  else if (tree==bpTree16) return BP_LUMI16;
  else if (tree==b0Tree16) return B0_LUMI16;
  else if (tree==ccTree16) return CC_LUMI16;
  else if (tree==udsTree16) return UDS_LUMI16;
  else if (tree==dpiTree16) return DPI_LUMI16;

  else if (tree==bpTree18) return BP_LUMI18;
  else if (tree==b0Tree18) return B0_LUMI18;
  else if (tree==ccTree18) return CC_LUMI18;
  else if (tree==udsTree18) return UDS_LUMI18;
  else if (tree==dpiTree18) return DPI_LUMI18;

  else return -1;
}
      
#endif
