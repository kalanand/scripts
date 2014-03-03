/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfGeneric.cc,v 1.1 2005/10/09 22:50:37 abi Exp $
 * Authors:
 *   WW, Wolfgang Walkowiak, UC Santa Cruz, walkowia@scipp.ucsc.edu
 * Description:
 *   Generic pdf wrapper for any RooAbsPdf.
 *   To be used in conjunction with RooFitCore/Models.
 * History:
 *   20-Mar-2002 WW Created initial version
 *   20-Sep-2002 WW Adjust code for compilation on SunOS 5.8
 *   04-Nov-2002 Abi Added dependents() capability
 *
 * Copyright (C) 2002 University of California and SLAC
 *****************************************************************************/

// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
//
// This class provides a generic wrapper for any RooAbsPdf.
//

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include "TString.h"

#include "RooFitCore/RooAbsPdf.hh"
#include "RooFitCore/RooDataSet.hh"

#include "BToDKTo3piK/BdkPdfGeneric.hh"

ClassImp(BdkPdfGeneric)

BdkPdfGeneric::BdkPdfGeneric() : _setDependentsCalled(kFALSE){
//
// Default constructor
//
}

BdkPdfGeneric::BdkPdfGeneric(const char * theName,
			     const char * theDesc,
			     RooAbsPdf & rooPdf, 
			     Bool_t owned,
			     RooArgSet * dependents)
{
// 
// Constructor
//
    init(theName, theDesc, rooPdf, owned, dependents);
}

BdkPdfGeneric::~BdkPdfGeneric() {
//
// Destructor
//
}

void BdkPdfGeneric::init(const char * theName,
			 const char * theDesc,
			 RooAbsPdf & rooPdf, 
			 Bool_t owned,
			 RooArgSet * dependents) {
// 
// Initialization
//
    baseInit(theName, theDesc);
    setPdf(rooPdf, owned);
    setDependents(dependents);
    initParameters();
}

void BdkPdfGeneric::setDependents(RooArgSet * deps) {
//
// Set the dependents
//
  if (0 == deps) {
    _setDependentsCalled = kFALSE;
    _dependents.removeAll();
  }
  else {
    _dependents.removeAll();
    _dependents.add(*deps);
  }
}

void BdkPdfGeneric::setPdf(RooAbsPdf & rooPdf, Bool_t owned) {
//
// Set the RooAbsPdf
//
    setIsValid(kFALSE);
    BdkPdfAbsBase::setPdf(rooPdf, owned);
}
  
void BdkPdfGeneric::initParameters() {
//
// Initialize parameters -- none for this one
//
    setIsValid(kFALSE);
}
 
void BdkPdfGeneric::createPdf() {
//
// Create the pdf -- nothing to be done here
//
    setIsValid(kTRUE);
}

RooArgSet BdkPdfGeneric::dependents() {
//
// Get dependents.
//
  if (kFALSE == _setDependentsCalled) {
    cerr 
      << GetName() << "::dependents(): BdkPdfGeneric knows nothing of" 
      << endl 
      << "   underlying PDF dependents unless setDependes() is called."
      << endl 
      << "   Returning empty RooArgSet." << endl;
    return RooArgSet();
  }

  return _dependents;
}
    


