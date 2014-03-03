/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 *
 * Package: 
 *
 *    File: $Id: BdkPdfCBShape.cc,v 1.2 2006/05/20 02:50:34 fwinkl Exp $
 *
 * Authors:
 *
 *   R. de Sangro, Laboratori Nazionali di Frascati INFN - riccardo.desangro@slac.stanford.edu
 *
 * Description:
 *   Class for a Crystall Ball shaped pdf with it's own variable definition
 *   To be used in conjunction with RooFitCore/Models.
 * History:
 *   04-Mar-2003 rid Created initial version from BdkPdfArgus from Abi Soffer
 *
 * Copyright (C) 2002 INFN - LNF 
 *****************************************************************************/
// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
// 
// Wrapper for Crystal Ball Shaped pdf.
// 

#include <iostream>
using std::cout;
using std::endl;

#include "TString.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitModels/RooCBShape.hh"

#include "BToDKTo3piK/BdkPdfCBShape.hh"

ClassImp(BdkPdfCBShape)

BdkPdfCBShape::BdkPdfCBShape() {
}

BdkPdfCBShape::BdkPdfCBShape(const char * theName,
			 const char * theDesc,
			 RooRealVar & dependent)
{
  init(theName, theDesc, dependent);
}

BdkPdfCBShape::~BdkPdfCBShape() {
}

BdkPdfCBShape::BdkPdfCBShape(BdkPdfCBShape & pdf) {
//
// Dummy copy constructor.
//
    cout << "BdkPdfCBShape::BdkPdfCBShape -- dummy copy constructor" << endl;
}

void BdkPdfCBShape::init(const char * theName,
		       const char * theDesc,
		       RooRealVar & dependent) {
  baseInit(theName, theDesc);
  setDependent(dependent);
  initParameters();
}

void BdkPdfCBShape::setDependent(RooRealVar & dependent) {
  setIsValid(kFALSE);
  _dependent = &dependent;
}

void BdkPdfCBShape::initParameters() {
  setIsValid(kFALSE);
  
  // now define the RooVars. The default parameters are good 
  // for missing mass:
  _m0 = new RooRealVar(TString(GetName())+".m0" ,"m0" , 1);
  
  _sigma = new RooRealVar(TString(GetName())+".sigma" ,"sigma", 1);

  _alpha = new RooRealVar(TString(GetName())+".alpha" ,"alpha", 1);

  _enne = new RooRealVar(TString(GetName())+".enne" ,"enne", 1);

  ((RooRealVar*)_m0)->setConstant(kFALSE);  
  ((RooRealVar*)_sigma)->setConstant(kFALSE);
  ((RooRealVar*)_alpha)->setConstant(kFALSE);
  ((RooRealVar*)_enne)->setConstant(kFALSE);

  // give resonable errors as first steps
  ((RooRealVar*)_m0)->setError(0.05);
  ((RooRealVar*)_sigma)->setError(0.05);
  ((RooRealVar*)_alpha)->setError(0.05);
  ((RooRealVar*)_enne)->setError(0.05);
}


void BdkPdfCBShape::createPdf() {  
  _thePdf =
    new RooCBShape(TString(GetName())+".pdf", TString(GetTitle()) + " CBShape",
		   *_dependent, *_m0, *_sigma, *_alpha, *_enne);
  setIsValid(kTRUE);
}

RooArgSet BdkPdfCBShape::dependents() {
  return RooArgSet(*_dependent);
}

void BdkPdfCBShape::linkParameters(const BdkPdfCBShape & pdf2) {
  setIsValid(kFALSE);
  linkParameters(pdf2._m0, pdf2._sigma, pdf2._alpha, pdf2._enne);
}

void BdkPdfCBShape::linkParameters(RooAbsReal *m0,
  				   RooAbsReal *sigma, 
				   RooAbsReal *alpha, 
				   RooAbsReal *enne){
  setIsValid(kFALSE);
  if ( m0 != 0 ) { _m0= m0;}
  if ( sigma != 0 ) { _sigma = sigma;} 
  if ( alpha != 0 ) { _alpha = alpha;} 
  if ( enne != 0 ) { _enne = enne;} 
}

