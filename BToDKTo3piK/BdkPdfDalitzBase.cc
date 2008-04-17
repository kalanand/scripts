/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfDalitzBase.cc,v 1.9 2007/04/18 11:57:06 fwinkl Exp $
 * Description:
 *   Base class for Dalitz PDF wrappers
 * History:
 *   18 Oct 2005, created, Abi soffer
 *
 * Copyright (C) 2005 Colorado State University and SLAC
 *****************************************************************************/
// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
// 
// Abstract base class wrapper for Dalitz pdfs. 
// 

#include <iostream>
using namespace std;

#include "TString.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooArgSet.hh"
#include "BToDKTo3piK/BdkPdfDalitzBase.hh"

ClassImp(BdkPdfDalitzBase)

// constructors:
BdkPdfDalitzBase::BdkPdfDalitzBase() {
  setIsValid(kFALSE);
}


BdkPdfDalitzBase::BdkPdfDalitzBase(const char * theName, const char * theDesc,
                                   RooAbsReal & m12, RooAbsReal & m13,
                                   BdkDalitzBase::Flavor flavor,
                                   BdkDalitz::Mode DdecMode) {
  // initialize:
  init(theName, theDesc, m12, m13, flavor, DdecMode);
}

// destructor:  
BdkPdfDalitzBase::~BdkPdfDalitzBase() {}


// Initializer:  
void BdkPdfDalitzBase::init(const char * theName, const char * theDesc,
                            RooAbsReal & m12, RooAbsReal & m13,
                            BdkDalitzBase::Flavor flavor,
                            BdkDalitz::Mode DdecMode) {
  baseInit(theName, theDesc);
  _m12 = &m12;
  _m13 = &m13;
  _flavor = flavor;
  _DdecMode = DdecMode;
  _effFunc = 0;
  setIsValid(kFALSE);
}


RooArgSet BdkPdfDalitzBase::dependents() {
  if (0 != _m12 && 0 != _m13) {
    return RooArgSet(*_m12, *_m13);
  }
  return RooArgSet();
}

/// Get efficiency function.
const BdkDalitzEff* BdkPdfDalitzBase::efficiencyFunc() const
{
  return _effFunc;
}

/// Set efficiency function.
void BdkPdfDalitzBase::setEfficiencyFunc(const BdkDalitzEff * f)
{
  _effFunc = (BdkDalitzEff*)f;
  setIsValid(kFALSE);
}

/// Link efficiency functions of two objects.
void BdkPdfDalitzBase::linkEfficiency(const BdkPdfDalitzBase& other)
{
  setEfficiencyFunc(efficiencyFunc());
  setIsValid(kFALSE);
}
