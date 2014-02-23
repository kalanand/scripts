/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfArgus.cc,v 1.2 2006/05/20 02:50:33 fwinkl Exp $
 * Authors:
 *   Abi Soffer, Colorado State University, abi@slac.stanford.edu
 * Description:
 *   Class for an ARGUS pdf with it's own variable definition
 *   To be used in conjunction with RooFitCore/Models.
 * History:
 *   21-Mar-2002 abi Created initial version
 *   28-Mar-2002 ww  add dependents(), delete generate()
 *   17-Sep-2002 WW Bug fixes for SunOS 5.8 compilation.
 *
 * Copyright (C) 2002 Colorado State University and SLAC
 *****************************************************************************/

// -- CLASS DESCRIPTION [IHFPDFWRAPPER] --
// 
// Wrapper for argus pdf.
// 

#include <iostream>
using std::cout;
using std::endl;
#include "TString.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitModels/RooArgusBG.hh"

#include "BToDKTo3piK/BdkPdfArgus.hh"

ClassImp(BdkPdfArgus)

BdkPdfArgus::BdkPdfArgus() {
}

BdkPdfArgus::BdkPdfArgus(const char * theName,
			 const char * theDesc,
			 RooRealVar & dependent)
{
  init(theName, theDesc, dependent);
}

BdkPdfArgus::~BdkPdfArgus() {
}

BdkPdfArgus::BdkPdfArgus(BdkPdfArgus & pdf) {
//
// Dummy copy constructor.
//
    cout << "BdkPdfArgus::BdkPdfArgus -- dummy copy constructor" << endl;
}

void BdkPdfArgus::init(const char * theName,
		       const char * theDesc,
		       RooRealVar & dependent) {
  baseInit(theName, theDesc);
  setDependent(dependent);
  initParameters();
}

void BdkPdfArgus::setDependent(RooRealVar & dependent) {
  setIsValid(kFALSE);
  _dependent = &dependent;
}

void BdkPdfArgus::initParameters() {
  setIsValid(kFALSE);
  
  // now define the RooVars. The default parameters are good 
  // for missing mass:
  _endPoint = new RooRealVar(TString(GetName())+".endPoint" ," endPoint" , 5.289);
  
  _exp = new RooRealVar(TString(GetName())+".exp" ," exp", -30);

  ((RooRealVar*)_endPoint)->setConstant(kFALSE);
  ((RooRealVar*)_exp)->setConstant(kFALSE);

  // give resonable errors as first steps
  ((RooRealVar*)_endPoint)->setError(0.05);
  ((RooRealVar*)_exp)->setError(0.05);
}


void BdkPdfArgus::createPdf() {  
  _thePdf =
    new RooArgusBG(TString(GetName())+".pdf", TString(GetTitle()) + " Argus",
		   *_dependent, *_endPoint, *_exp);

  setIsValid(kTRUE);
}

RooArgSet BdkPdfArgus::dependents() {
  
  return RooArgSet(*_dependent);
}

void BdkPdfArgus::linkParameters(const BdkPdfArgus & pdf2) {
  setIsValid(kFALSE);
  linkParameters(pdf2._endPoint , pdf2._exp);
}

void BdkPdfArgus::linkParameters(RooAbsReal *endPoint,
				 RooAbsReal *exp) {
  setIsValid(kFALSE);
  if ( endPoint != 0 ) {_endPoint = endPoint;}
  if ( exp != 0 ) {_exp = exp;} 
}

