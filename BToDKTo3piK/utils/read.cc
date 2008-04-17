// $Id: read.cc,v 1.19 2006/05/08 23:13:19 fwinkl Exp $
// Reads the data into a RooDataSet and returns the dataset. The caller
// owns the returned pointer.

#ifndef READ_CC
#define READ_CC

#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCut.h"
#include "TRandom.h"

#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooRealVar.hh"

#include "../BToDKTo3piK/globals/cuts.hh"
#include "../BToDKTo3piK/globals/vars.hh"

//----------------------------------------------------------
// Creates a RooDataSet from a tree in an efficient way
// by first applying the cut on the tree and then reading
// in the variables.
RooDataSet *dataSetFromTree(TTree *tree, TCut cut, 
                            const RooArgSet *vars = 0) {
  gROOT->cd();
  TTree *cutTree = tree;

  if (strlen(cut.GetTitle())>0) {
    cout << "read(): Reading "<<tree->GetName();
    TDirectory* dir = tree->GetDirectory();   // gives NULL for TChains
    if (dir) cout << " from "<<dir->GetName();
    cout << endl << "read(): applying cut " << cut.GetTitle() << endl;
    cutTree = tree->CopyTree(cut);
  }

  const RooArgSet *varSet = allVars;
  if (vars) varSet = vars;

  RooDataSet* dataCut = new RooDataSet("data","data",cutTree,*varSet);
  cout << "read(): Read " << dataCut->numEntries() << " events" << endl;

  if (strlen(cut.GetTitle())>0) delete cutTree;

  dataCut->SetName(tree->GetName());
  dataCut->SetTitle(tree->GetTitle());

  return dataCut;
}


//------------------------------------------------------
// Change dataRead if possible, otherwise make a copy
RooDataSet * refineRead(RooDataSet* dataRead) {

  // Add a random number for making random chunks:
  RooDataSet * outData = 0;

  // Add randRep if it doesn't already exist:
  if (0 == dataRead->get()->find(randAdd->GetName())) {

    outData = (RooDataSet*)dataRead->emptyClone();
    outData->addColumn(*randAdd);
    
    TRandom ran;
    for (int i = 0; i < dataRead->numEntries(); ++i){
      // Get the original event:
      RooArgSet outEvt(*dataRead->get(i));
      outEvt.add(*randAdd);
      RooRealVar * tmpRan = (RooRealVar *)outEvt.find(randAdd->GetName());
      tmpRan->setVal( ran.Rndm() );
      outData->add(outEvt);
    }
  }
  else {
    outData = dataRead;
  }

  // Add columns with square Dalitz variables
  outData->addColumn(*s12);
  outData->addColumn(*s13);
  outData->addColumn(*s23);
  outData->addColumn(*s12mc);
  outData->addColumn(*s13mc);
  // Add transformed NN net variable
  outData->addColumn(*qprimeF);
  outData->addColumn(*dprimeF);

  return outData;
}


// Apply cut and read 'vars'
RooDataSet * read(TTree * tree, TTree * tree2,
                  TCut cut, RooArgSet* vars) {
  gROOT->cd();
  if (!vars) vars = allVars;
  RooDataSet* dataRead = dataSetFromTree(tree, cut, vars);

  if (0 != tree2) {
    gROOT->cd();    
    RooDataSet* dataRead2 = dataSetFromTree(tree2, cut, vars);
    dataRead->append(*dataRead2);
  }

  RooDataSet* refined = refineRead(dataRead);

  return refined;
}

// Apply global 'readCut' and read 'allVars' by default
RooDataSet * read(TTree * tree, TTree * tree2 = 0) {

  return read(tree, tree2, readCut, allVars);
}

#endif
