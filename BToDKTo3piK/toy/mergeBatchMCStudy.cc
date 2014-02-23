// $Id: mergeBatchMCStudy.cc,v 1.9 2007/05/11 14:27:34 fwinkl Exp $
// Merge the individual fit result files from BdkBatchMCStudy
// Only parameters floating in the fit are considered

#include <fstream>

#include "TFile.h"
#include "TString.h"
#include "TIterator.h"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooAbsArg.hh"

RooArgSet* getFloatingVars(const RooArgSet& set);


RooDataSet* mergeBatchMCStudy(const char* fileList, const char* mergeFile = 0,
                              Bool_t mergeAllVars = kFALSE)
{
  RooDataSet *data = 0;
  ifstream files;
  Int_t nFiles = 0;
  files.open(fileList);
  while (!files.eof()) {
    string fitFile;
    files >> fitFile;   
    if (fitFile=="") continue;
    TFile f(fitFile.c_str());   
    if (f.IsZombie()) {
      cout << "Cannot open file " << fitFile << endl;
      continue;
    }

    RooDataSet* d = (RooDataSet*)f.Get("fitParData");

    if (!data) {
      RooArgSet* set;
      if (mergeAllVars) set = (RooArgSet*)d->get();
      else set = getFloatingVars(*d->get());

      cout << "Merging the following variables:"<<endl;
      set->Print();
      data = new RooDataSet(d->GetName(),d->GetTitle(),*set);
      if (!mergeAllVars) delete set;
    }
    

    cout << fitFile <<" contains "<<d->numEntries()<<" fit results."<<endl;
    gROOT->cd();
    data->append(*d);

    nFiles++;
    f.Close();
  }
  cout << "Total: " << data->numEntries() <<" fit results from " << nFiles << " files."<<endl;
  if (mergeFile) {
    cout << "Writing all fit results to "<<mergeFile<<endl;
    TFile f(mergeFile,"recreate");
    data->Write();
    f.Close();
  }
  return data;
}

// Find the floating variables in set including pull and error
RooArgSet* getFloatingVars(const RooArgSet& set)
{
  RooArgSet* subset = new RooArgSet();
  TIterator* iter = set.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (!arg->isConstant()) {
      subset->add(*arg);
      TString name = arg->GetName();
      RooAbsArg* pull = set.find(name+"pull");
      if (!pull) cout << "Cannot find pull variable for "<<name<<endl;
      else subset->add(*pull);
      
      RooAbsArg* err = set.find(name+"err");
      if (!err) cout << "Cannot find error variable for "<<name<<endl;
      else subset->add(*err);
    }
  }
  return subset;
}
