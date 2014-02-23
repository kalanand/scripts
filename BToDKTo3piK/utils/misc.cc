// $Id: misc.cc,v 1.7 2007/04/02 14:14:31 fwinkl Exp $
// Misc routines

#ifndef MISC_CC
#define MISC_CC

#include "TTree.h"
#include "TIterator.h"
#include "TString.h"

#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAbsCollection.hh"


Double_t sqr(Double_t d)
{
  return d*d;
}

// Some statistics functions to be used on TTree or RooDataSet

// Median
Double_t getMedian(TTree& tree, const char* var)
{
  // We use draw to access the array of double's via TTree::GetV1()
  tree.Draw(var,"","goff");   // g(raphics)off
  return TMath::Median(tree.GetSelectedRows(),tree.GetV1());
}

Double_t getMedian(const RooDataSet& data, const RooRealVar& var)
{
  return getMedian((TTree&)data.tree(),var.GetName());
}



// Mean
Double_t getMean(TTree& tree, const char* var)
{
  tree.Draw(var,"","goff");
  return TMath::Mean(tree.GetSelectedRows(),tree.GetV1());
}

Double_t getMean(const RooDataSet& data, const RooRealVar& var)
{
  return getMean((TTree&)data.tree(),var.GetName());
}


// RMS
Double_t getRMS(TTree& tree, const char* var)
{
  tree.Draw(var,"","goff");
  return TMath::RMS(tree.GetSelectedRows(),tree.GetV1());
}

Double_t getRMS(const RooDataSet& data, const RooRealVar& var)
{
  return getRMS((TTree&)data.tree(),var.GetName());
}


// Set RooRealVar (and errors) to zero
void zero(RooRealVar& r)
{
  r.setVal(0);
  if (r.hasError()) r.setError(0);
  if (r.hasAsymError()) r.setAsymError(0,0);
}


// Set all RooRealVars in col (and errors) to zero
void zero(RooAbsCollection& coll)
{ 
  TIterator* iter = coll.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())!="RooRealVar") continue;
    zero(*(RooRealVar*)arg);
  }
  delete iter;
}

// Fix all RooRealVars in coll
void fixAll(RooAbsCollection& coll, Bool_t fix = kTRUE)
{
  TIterator* iter = coll.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())!="RooRealVar") continue;
    ((RooRealVar*)arg)->setConstant(fix);
  }
  delete iter;
}

// Remove all ranges from RooRealVars in coll
void removeAllRanges(RooAbsCollection& coll)
{
  TIterator* iter = coll.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())!="RooRealVar") continue;
    ((RooRealVar*)arg)->removeRange();
  }
  delete iter;
}


// Remove all errors from RooRealVars in coll
void removeAllErrors(RooAbsCollection& coll)
{
  TIterator* iter = coll.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())!="RooRealVar") continue;
    ((RooRealVar*)arg)->removeError();
  }
  delete iter;
}


// Copy error to value and remove error
void copyErrorToVal(RooAbsCollection& coll)
{
  TIterator* iter = coll.createIterator();
  RooRealVar* r;
  while ((r = (RooRealVar*)iter->Next())) {
    if (TString(r->ClassName())!="RooRealVar") continue;
    r->setVal(r->getError());
    r->removeError();
  }
  delete iter;
}


// Increment RooRealVar

void inc(RooRealVar& r, Double_t n = 1)
{
  r.setVal(r.getVal() + n);
}

// Add two RooRealVar
void add(RooRealVar& result, const RooRealVar& r1, const RooRealVar& r2,
         Bool_t addError = kFALSE)
{
  result.setVal(r1.getVal()+r2.getVal());
  if (addError) result.setError(r1.getError()+r2.getError());
  else result.setError(sqrt(sqr(r1.getError())+sqr(r2.getError())));
}

// Add three RooRealVar
void add(RooRealVar& result,
	 const RooRealVar& r1, const RooRealVar& r2, const RooRealVar& r3,
         Bool_t addError = kFALSE)
{
  result.setVal(r1.getVal()+r2.getVal()+r3.getVal());
  if (addError) result.setError(r1.getError()+r2.getError()+r3.getError());
  else result.setError(sqrt(sqr(r1.getError())
                            +sqr(r2.getError())
                            +sqr(r3.getError())));
}

// Subtract two RooRealVar: r1 - r2
void sub(RooRealVar& result, const RooRealVar& r1, const RooRealVar& r2,
         Bool_t subError = kFALSE)
{
  result.setVal(r1.getVal() - r2.getVal());
  if (subError) result.setError(r1.getError()-r2.getError());
  else result.setError(sqrt(sqr(r1.getError())+sqr(r2.getError())));
}


// multiply two RooRealVar
void multiply(RooRealVar& result, const RooRealVar& r1, const RooRealVar& r2)
{
  result.setVal(r1.getVal()*r2.getVal());

  if (r1.getVal()==0 || r2.getVal()==0) result.setError(0);
  else result.setError(result.getVal()* sqrt(sqr(r1.getError()/r1.getVal()) +
                                             sqr(r2.getError()/r2.getVal())));
}

// Divide two RooRealVar
Bool_t divide(RooRealVar& result, const RooRealVar& dividend, const RooRealVar& divisor)
{
  if (divisor.getVal()==0) return kFALSE;
  result.setVal(dividend.getVal()/divisor.getVal());
  if (result.getVal()==0) result.setError(0);
  else result.setError(result.getVal()*
                       sqrt(sqr(dividend.getError())/sqr(dividend.getVal()) +
                            sqr(divisor.getError())/sqr(divisor.getVal())));

  return kTRUE;
}


// Add RooRealVars in coll to sum.
// sum must have the same members as coll.
Bool_t add(RooAbsCollection& sum, const RooAbsCollection& coll,
           Bool_t addError = kFALSE)
{
  if (sum.getSize()!=coll.getSize()) return false;
  if (!sum.equals(coll)) return false;

  TIterator* iter = coll.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())!="RooRealVar") continue;
    
    RooRealVar* r = (RooRealVar*)arg;
    RooRealVar* rsum = (RooRealVar*)sum.find(r->GetName());
    if (rsum==0) continue;
    
    add(*rsum,*rsum,*r,addError);    
  }
  delete iter;
  return true;
}


// Subtract RooRealVars in coll from diff.
// diff must have the same members as coll.
Bool_t sub(RooAbsCollection& diff, const RooAbsCollection& coll,
           Bool_t subError = kFALSE)
{
  if (diff.getSize()!=coll.getSize()) return false;
  if (!diff.equals(coll)) return false;

  TIterator* iter = coll.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())!="RooRealVar") continue;
    
    RooRealVar* r = (RooRealVar*)arg;
    RooRealVar* rdiff = (RooRealVar*)diff.find(r->GetName());
    if (rdiff==0) continue;
    
    sub(*rdiff,*rdiff,*r,subError);    
  }
  delete iter;
  return true;
}


// Divide RooRealVars in divident by RooRealVars in divisor
// Members must be the same
Bool_t divide(RooAbsCollection& dividend, const RooAbsCollection& divisor)
{
  if (divisor.getSize()!=dividend.getSize()) return false;
  if (!divisor.equals(dividend)) return false;

  TIterator* iter = dividend.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())!="RooRealVar") continue;
    
    RooRealVar* rDividend = (RooRealVar*)arg;
    RooRealVar* rDivisor = (RooRealVar*)divisor.find(arg->GetName());
    if (rDivisor==0) continue;
    
    divide(*rDividend,*rDividend,*rDivisor);
  }
  delete iter;
  return true;
}



// Calculate mean of RooRealVars in list that have the same name.
// User owns return value.
RooArgSet* mean(const RooArgList& list)
{
  RooArgSet* sum = new RooArgSet();
  RooArgSet* N = new RooArgSet();

  for (int i=0; i<list.getSize(); i++) {
    RooRealVar* r;
    if (!(r = dynamic_cast<RooRealVar*>(list.at(i)))) return 0;
    
    if (!N->contains(*r)) {
      RooRealVar* rTemp;
      rTemp = (RooRealVar*)N->addClone(*r);
      rTemp->setVal(1);
      rTemp->removeRange();
      rTemp = (RooRealVar*)sum->addClone(*r);
      rTemp->setError(r->getError());
      rTemp->removeRange();
    }
    else {
      inc(*(RooRealVar*)N->find(r->GetName()));
      RooRealVar* rsum = (RooRealVar*)sum->find(r->GetName());
      add(*rsum,*rsum,*r,true);
    }
  }

  TIterator* iter = sum->createIterator();
  RooRealVar* r;
  while ((r = (RooRealVar*)iter->Next())) {
    Double_t n = ((RooRealVar*)N->find(r->GetName()))->getVal();

    r->setVal(r->getVal()/n);
    r->setError(r->getError()/n);
  }
  delete iter;
  return sum;
}

// Calculate mean of RooAbsCollections in list.
// User owns return value.
RooAbsCollection* mean(const TList& list)
{
  if (list.GetSize()==0) return 0;
  
  TIterator* iter = list.MakeIterator();
  TObject* obj;
  RooAbsCollection* sum = 0;
  int N = 0;
  while ((obj = iter->Next())) {
    RooAbsCollection* coll;
    if (!(coll = dynamic_cast<RooAbsCollection*>(obj))) return 0;
    if (sum==0) {
      sum = (RooAbsCollection*)coll->snapshot(kFALSE);
      zero(*sum);
    }
    if (add(*sum,*coll,true)) ++N;
    else return 0;
  }      
  delete iter;

  iter = sum->createIterator();
  while ((obj = iter->Next())) {
    RooRealVar* r;
    if (!(r = dynamic_cast<RooRealVar*>(obj))) return 0;
    r->setVal(r->getVal()/N);
    r->setError(r->getError()/N);
  }
  delete iter;

  return sum;
}


// Square all variables in coll
void sqr(RooAbsCollection& coll)
{
  TIterator* iter = coll.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())!="RooRealVar") continue;
    
    RooRealVar* r = (RooRealVar*)arg;
    r->setVal(r->getVal()*r->getVal());
  }
  delete iter;
}

// Take sqrt of all variables in coll
void sqrt(RooAbsCollection& coll)
{
  TIterator* iter = coll.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())!="RooRealVar") continue;
    
    RooRealVar* r = (RooRealVar*)arg;
    if (r->getVal()<0) {
      cout << "Cannot take sqrt of ";
      r->printToStream(cout);
    }
    r->setVal(sqrt(r->getVal()));
  }
  delete iter;
}

// Multiply all variables in coll with r
void multiply(RooAbsCollection& coll, const RooRealVar& r)
{
  TIterator* iter = coll.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())!="RooRealVar") continue;
    
    RooRealVar* var = (RooRealVar*)arg;
    multiply(*var, *var, r);
  }
  delete iter;
}

// Add suffix to all RooAbsArgs in coll
void addSuffix(RooAbsCollection& coll, const char* suffix)
{
  if (suffix==0) return;
  TIterator* iter = coll.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)iter->Next())) {
    arg->SetName(TString(arg->GetName())+TString(suffix));
  }
  delete iter;
}

#endif
