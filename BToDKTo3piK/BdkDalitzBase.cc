/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkDalitzBase.cc,v 1.18 2006/04/20 16:47:08 fwinkl Exp $
 * Description:
 *   Base class for Dalitz PDFs
 * History:
 *   26 Oct 2005, created, Frank Winklmeier
 *
 * Copyright (C) 2005 Colorado State University and SLAC
 *****************************************************************************/

#include "BToDKTo3piK/BdkDalitzBase.hh"
#include "BToDKTo3piK/BdkDalitzEff.hh"

ClassImp(BdkDalitzBase)


// constructors:

BdkDalitzBase::BdkDalitzBase(const char * name, const char * title,
                             BdkDalitzBase::Flavor flavor, 
			     BdkDalitzBase::Mode DdecMode) :
  RooAbsPdf(name, title),
  BdkDalitz(flavor, DdecMode)
{
  setEfficiencyFunc(0);
}

BdkDalitzBase::BdkDalitzBase(const BdkDalitzBase& other, const char * name) :
  RooAbsPdf(other, name),
  BdkDalitz(other),
  _efficiencyFunc(other._efficiencyFunc),
  _effFirstEval(other._effFirstEval)
{
}

// destructor:  
BdkDalitzBase::~BdkDalitzBase() 
{
}


void BdkDalitzBase::setEfficiencyFunc(const BdkDalitzEff * f) 
{
  _efficiencyFunc = f;
  _effFirstEval = kTRUE;
}

/// Evaluate the efficiency function
/// Return 1 if none is assigned.
Double_t BdkDalitzBase::efficiency(Double_t m12, Double_t m13) {

  // Print message about global efficiency function on first evaluation
  // Only if RooAbsPdf::verboseEval()
  if (_verboseEval) {
    if (_effFirstEval) {
      if (_efficiencyFunc) {
        cout << GetName()<<" uses this efficiency function:"<<endl;
        _efficiencyFunc->getParameters(RooArgSet())->Print("v");
      }
      else cout << GetName()<<" does not use an efficiency function."<<endl;
      _effFirstEval = kFALSE;
    }
  }
  
  if (0 != _efficiencyFunc) {    
    return _efficiencyFunc->evaluateAt(m12, m13);
  }
  return 1.0;
}
