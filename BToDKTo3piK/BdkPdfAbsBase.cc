/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfAbsBase.cc,v 1.4 2006/04/18 01:45:41 kalanand Exp $
 * Authors:
 *   WW, Wolfgang Walkowiak, UC Santa Cruz, walkowia@scipp.ucsc.edu
 * Description:
 *   Abstract base class for pdfs which define their own variables.
 *   To be used in conjunction with RooFitCore/Models.
 * History:
 *   18-Mar-2002 WW Created initial version
 *   17-Sep-2002 WW Bug fix for SunOS 5.8 compiler.
 *
 * Copyright (C) 2001 University of California and SLAC
 *****************************************************************************/

// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
//
// This class is the base class for RooFit pdf wrappers.
// It provides holds the pdf pointer and 
// provides the basic accessor methods getPdf(), getArgSet() and generate().
// In addition, utility methods used for the generation of 
// events for simultaneous pdfs are implemented in this base class.
//
#include <ctype.h>
#include <iostream>
using std::cout;
using std::endl;
#include <assert.h>

#include "TString.h"
#include "RooFitCore/RooAbsPdf.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooAbsArg.hh"
#include "RooFitCore/RooAbsCollection.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooRandom.hh"
#include "RooFitCore/RooRealVar.hh"

#include "BToDKTo3piK/BdkPdfAbsBase.hh"

ClassImp(BdkPdfAbsBase)

static const char rcsid[] = "$Id: BdkPdfAbsBase.cc,v 1.4 2006/04/18 01:45:41 kalanand Exp $";


const RooArgSet * BdkPdfAbsBase::_allAnalysisVars = 0;


BdkPdfAbsBase::BdkPdfAbsBase() : 
  _thePdf(0), 
  _isValid(kFALSE) , 
  _pdfOwned(kFALSE) // _thePdf is not yet owned since there is no _thePDf yet. 
{
    //
    // Constructor
    //

    // cout << "BdkPdfAbsBase::BdkPdfAbsBase - constructor" << endl;
}

BdkPdfAbsBase::~BdkPdfAbsBase() {
    //
    // Destructor
    //

  if (verbose().Contains("d")) {
    cout << GetName() << "::~BdkPdfAbsBase called." << endl;
  }

  deletePdf();
}

BdkPdfAbsBase::BdkPdfAbsBase(BdkPdfAbsBase & pdf) {
    //
    // Dummy copy constructor
    //
    cout << "BdkPdfAbsBase::BdkPdfAbsBase: dummy copy constructor!" << endl;
    assert(0);
}

BdkPdfAbsBase & BdkPdfAbsBase::operator=(BdkPdfAbsBase & pdf) {
    //
    // Dummy assignment operator
    //
    cout << "forbidden BdkPdfAbsBase::operator= !" << endl;
    assert(0);
    return pdf;
}

void BdkPdfAbsBase::baseInit(const char *theName,
			     const char *theDesc) {
//
//  base initialization
// 
    setIsValid(kFALSE);
    SetNameTitle(theName,theDesc);
    deletePdf();
}         

void BdkPdfAbsBase::initParameters() {
}


RooAbsPdf* BdkPdfAbsBase::getPdf() {
// 
// public accessor to resulting pdf
//
  if (verbose().Contains("g")) {
    cout << GetName() << ": getPdf() called." << endl;
  }
  if ( !getIsValid() ) {
    if (verbose().Contains("c")) {
      cout << GetName() << ": Calling createPdf()." << endl;
    }
    deletePdf();
    createPdf();
  }

  return _thePdf;
}

RooArgSet *BdkPdfAbsBase::getArgSet() {
// 
// public accessor to all parameters + dependents in RooArgSet
//
  return getPdf()->getParameters(RooArgSet());
}

RooAbsArg *BdkPdfAbsBase::find(const char *name) {
// 
// public accessor to one parameter of pdf by name
//
  return getArgSet()->find(name);
}

RooArgSet BdkPdfAbsBase::parameters() {
//
// Returns all the parameters (non-dependents)
//

  RooArgSet deps = dependents();
  // Add the full list of all variables:
  if (0 != allAnalysisVars()){
    deps.add(*allAnalysisVars());
  }

  return *(getPdf()->getParameters(&deps));
}

RooArgSet BdkPdfAbsBase::parametersFree(Bool_t freeOnly) {
//
// Returns only the free (floating) parameters if freeOnly = kTRUE,
// and only the fixed parameters if freeOnly = kFALSE.
//
  RooArgSet result;

  RooArgSet allPars = parameters();
  TIterator * iter = allPars.createIterator();
  RooRealVar * var;
  while(0!= (var= (RooRealVar*)iter->Next())) {
    if ((freeOnly && !var->isConstant()) ||
	(!freeOnly && var->isConstant())) {
      result.add(*var);
    }
  }
  delete iter;

  return result;
}


RooDataSet *BdkPdfAbsBase::generate(Int_t nEvents, RooDataSet* protoData) {
// 
// Generic generate method.  
// Calls RooAbsPdf::generate().
//
  RooDataSet * data = 0;
  RooArgSet deps = dependents();
  if (0 == protoData) {
    data = getPdf()->generate(deps, nEvents);
  }
  else {
    data = getPdf()->generate(deps, *protoData, nEvents);
  }
  
  return data;
}



RooAbsArg* BdkPdfAbsBase::findReplace(RooAbsArg& arg, 
				      RooAbsCollection& set) {
//
// Replace arg by argument with the same name ending "xxx.ending"
// from set.
//
    RooAbsArg *resArg = 0;
    TString name = arg.GetName();
    TString suf  = ( name.Contains('.') ? 
		     "*"+ name(name.Last('.'),name.Length()-1) : name );
    RooAbsCollection* set2 = set.selectByName(suf.Data());

    if ( set2->getSize() > 0 ) {
	Bool_t success = kFALSE;
	if ( set2->getSize() > 1 ) {
	    cout << "BdkPdfAbsBase::findReplace: too many matching args!" 
		 << endl;
	    cout << "BdkPdfAbsBase::findReplace: taking first of " 
		 << set2->getSize() << " args." << endl;
	} 
	TIterator *iter = set2->createIterator();
	resArg = (RooAbsArg*)iter->Next();
	success = kTRUE;
	delete iter;
	cout << GetName() << ": Link " << arg.GetName() << " --> " ; 
	if ( success ) 
	    cout << resArg->GetName() << " ..."  << " success !" << endl;
	else
	    cout << " ... failed !!!" << endl;
    } else {
	resArg = &arg;
	/// cout << "BdkPdfAbsBase::findReplace: no args in set2!" << endl;
    }
    delete set2;
    
    return resArg;
}

Bool_t BdkPdfAbsBase::getIsValid() const {
// 
// get validity
//
    return _isValid;
}

void BdkPdfAbsBase::setIsValid(Bool_t isValid) {
//
// set validity
//
    _isValid = isValid;

    if (verbose().Contains("v")) {
      cout << GetName() << "::setIsValid(" << (int)isValid << ") called."
	   << endl;
    }
}


void BdkPdfAbsBase::fixAll(Bool_t val) {
//
// Fix or float all the parameters
//
  RooArgSet pars = parameters();
  TIterator * iter = pars.createIterator();
  RooRealVar * var;
  while(0 != (var= (RooRealVar*)iter->Next())) {
    // Fix/free the parameter:
    var->setConstant(val);
    
    if (kFALSE == val) {
      // Always fix the dummy parameters which are numerical only:
      TString theName = var->GetName();
      Bool_t numeric = kTRUE;
      int numDots = 0;
      for (int d = 0; d < theName.Length(); ++d){
	if (!isdigit(theName(d))) {
	  // the character is not a digit. See if it's a dot:
	  if (theName(d) == '.') {
	    ++numDots;
	  }
	  else {
	    // otherwise, found a non-digit and non-'.', so not numerical
	    numeric = kFALSE;
	    break; 
	  }
	} // end if character is not a digit
      } // end loop on characters in theName
      if (numDots <= 1 && numeric == kTRUE){
	// Only digits or one '.' found. So numerical. Fix it:
	cout << "fixAll() fixing numerical parameter " << theName << endl;
	var->setConstant();
      }
    } // end if (kFALSE == val)
  } // end while iterating on parameters
  delete iter;
}

void BdkPdfAbsBase::setVerbose(const char * val) {
  //
  // set verbose flag:
  // c: Message every time createPdf is called
  // i: Informational messages
  // d: Message every time object is deleted
  // g: Message when getPdf() is called
  //
  _verbose = val;
  _verbose.ToLower();
}

void BdkPdfAbsBase::setPdf(RooAbsPdf & thePdf, Bool_t owned) {
  // Set _thePdf.
  // First, clean up the old PDF:
  deletePdf();

  // Then assign it to the argument:
  _thePdf = &thePdf;
  _pdfOwned = owned;
}


void BdkPdfAbsBase::deletePdf() {
  // Delete _thePdf if it's owned and non-0:
  if (kTRUE == _pdfOwned && 0 != _thePdf) {
    if (verbose().Contains("d")) {
      cout << GetName() << "::deletePdf(): deleting _thePdf." << endl;
    }

    delete _thePdf;
    
    // Since this function can get called not in the destructor, make
    // sure that _thePdf is null:
    _thePdf = 0;
  }
}

