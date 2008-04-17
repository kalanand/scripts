/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfPolyn.cc,v 1.2 2006/03/02 01:06:03 fwinkl Exp $
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

// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
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
#include "RooFitCore/RooArgList.hh"
#include "RooFitModels/RooPolynomial.hh"
#include "RooFitCore/RooFormulaVar.hh" 

#include "BToDKTo3piK/BdkPdfPolyn.hh"

ClassImp(BdkPdfPolyn)

BdkPdfPolyn::BdkPdfPolyn() {
}

BdkPdfPolyn::BdkPdfPolyn(const char * theName,
			 const char * theDesc,
			 RooRealVar & dependent,
                         Int_t lowestOrder,
			 Int_t highestOrder) {

  init(theName, theDesc, dependent, lowestOrder, highestOrder);
}

BdkPdfPolyn::~BdkPdfPolyn() {
}

BdkPdfPolyn::BdkPdfPolyn(BdkPdfPolyn & pdf) {
//
// Dummy copy constructor.
//
    cout << "BdkPdfPolyn::BdkPdfPolyn -- dummy copy constructor" << endl;
}

void BdkPdfPolyn::init(const char * theName,
		       const char * theDesc,
		       RooRealVar& dependent,
                       Int_t lowestOrder,
                       Int_t highestOrder) {
  baseInit(theName, theDesc);
  useOffset(kTRUE);
  setDependent(dependent);
  setlowestOrder(lowestOrder);
  sethighestOrder(highestOrder);
  initParameters();
}

void BdkPdfPolyn::initParameters() {
  setIsValid(kFALSE);
  _offset = new RooRealVar(TString(GetName())+".offset",
			  TString(GetTitle())+" Polynomial offset", 
                          0, -100., 100.);  

  _diff = new RooFormulaVar(TString(GetName())+ ".diff", 
			    TString(GetTitle()) + " offset", "@0-@1",  
			    RooArgList(*_dependent, *_offset)); 

  for(Int_t i=_lowestOrder; i<=_highestOrder;i++) {
	TString pdfname = TString(GetName()) + ".coef";
        pdfname += i ;
	RooRealVar *temp = new RooRealVar(pdfname,"coef", 1);
        temp->setConstant(kFALSE);
        temp->setError(0.01);
        _coefList.addOwned(*temp);
  } 
}

void BdkPdfPolyn::setDependent(RooRealVar & dependent) {
  setIsValid(kFALSE);
  _dependent = &dependent;
}

void BdkPdfPolyn::setlowestOrder(Int_t lowestOrder) {
   setIsValid(kFALSE);
   _lowestOrder = lowestOrder;
} 

void BdkPdfPolyn::sethighestOrder(Int_t highestOrder) {
   setIsValid(kFALSE);
   _highestOrder = highestOrder;
}
	

RooRealVar* BdkPdfPolyn::getCoef(Int_t n) {
  return (RooRealVar*)_coefList.at(n - _lowestOrder);
}
 
void BdkPdfPolyn::setcoefConst(Int_t n, Double_t val,Bool_t trfk ) {
  setIsValid(kFALSE);
  RooRealVar* tmp = getCoef(n);
   tmp->setConstant(trfk);
   tmp->setVal(val);
}

void BdkPdfPolyn::useOffset(Bool_t use) {
  _useOffset = use;
  setIsValid(kFALSE);
}

void BdkPdfPolyn::createPdf() { 
  RooAbsReal * var = _diff;
  if (kFALSE == _useOffset) {
    var = _dependent;
  }

  _thePdf =
    new RooPolynomial(TString(GetName())+".pdf", 
		      TString(GetTitle()) + ".Polyn",
		      *var, _coefList , _lowestOrder);
  
  setIsValid(kTRUE);
}

RooArgSet BdkPdfPolyn::dependents() {
  
  return RooArgSet(*_dependent);
}

void BdkPdfPolyn::linkParameters(const BdkPdfPolyn & pdf2) {
  setIsValid(kFALSE);
  setParameters(pdf2._coefList , pdf2._lowestOrder);
}

void BdkPdfPolyn::linkParameters(Int_t lowestOrder, Int_t highestOrder) {
  setIsValid(kFALSE);
  setlowestOrder(lowestOrder);
  sethighestOrder(highestOrder);
}
	
void BdkPdfPolyn::setParameters(const RooArgList & coefList,Int_t lowestOrder) {
  setIsValid(kFALSE);
  if ( coefList.getSize() > 0 ) {_coefList = coefList;}
  if ( lowestOrder >= 0 ) {_lowestOrder = lowestOrder;} 
  sethighestOrder(_coefList.getSize()+_lowestOrder+1);
}

Double_t BdkPdfPolyn::analyticalIntegral(Int_t code) 
{
  return _thePdf->analyticalIntegral(code);
}

