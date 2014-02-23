/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 *
 * Package: 
 *
 *    File: $Id: BdkPdfCBArgus.cc,v 1.1 2005/10/09 22:50:36 abi Exp $
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
// Wrapper for Crystall Ball shape plus argus pdf.
// 

#include <iostream>
using std::cout;
using std::endl;

#include "TString.h"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooArgSet.hh"

#include "BToDKTo3piK/BdkPdfCBArgus.hh"

ClassImp(BdkPdfCBArgus)

//--------------------------------------------------------------
BdkPdfCBArgus::BdkPdfCBArgus() : _fracCBShape(0) {
  // Default Constructor:
  addComponents();
}

//--------------------------------------------------------------
BdkPdfCBArgus::BdkPdfCBArgus(const char * theName,
				   const char * theDesc,
				   RooRealVar & dependent)  : _fracCBShape(0) {
  // Initializing Constructor:
  addComponents();
  init(theName, theDesc, dependent);
}

//--------------------------------------------------------------
BdkPdfCBArgus::~BdkPdfCBArgus() {
  // Destructor:
  deleteFrac();
}

//--------------------------------------------------------------
BdkPdfCBArgus::BdkPdfCBArgus(const BdkPdfCBArgus & pdf) {
// Dummy copy constructor:
    cout << "BdkPdfCBArgus::BdkPdfCBArgus -- dummy copy constructor"
	 << endl;
}


//--------------------------------------------------------------
void BdkPdfCBArgus::init(const char * theName,
			    const char * theDesc,
			    RooRealVar & dependent) {
  // Initializer:
  _fracCBShape = 0;

  baseInit(theName, theDesc);
  setDependent(dependent);
  initParameters();
}

//--------------------------------------------------------------
void BdkPdfCBArgus::setDependent(RooRealVar & dependent) {
  // Sets the dependent variable:
  _dependent = &dependent;
  _cbShape.setDependent(dependent);
  _argus.setDependent(dependent);

  setIsValid(kFALSE);
}

//--------------------------------------------------------------
void BdkPdfCBArgus::initParameters() {
  // Initializes the contained RRV's and PDF wrappers:
  _fracCBShape = new RooRealVar(GetName() + TString(".fracCBShape"),
			      GetTitle () + TString(" fracCBShape"),
			      0.1, 0, 1);

  ((RooRealVar*)_fracCBShape)->setError(0.01);

  _cbShape.init(GetName() + TString(".cbShape"), 
	      GetTitle() + TString(" Crystall Ball"), 
	      *_dependent);

  _argus.init(GetName() + TString(".argus"), 
	      GetTitle() + TString(" Argus"), 
	      *_dependent);
}

//--------------------------------------------------------------
void BdkPdfCBArgus::createPdf() {
  // Builds the PDF. The PDFs of the contained wrappers are automatically built
  _thePdf = new RooAddPdf(GetName() + TString(".thePdf"), 
			  GetTitle() + TString(" thePdf"), 
			  *_cbShape.getPdf(), 
			  *_argus.getPdf(),
			  *_fracCBShape);

  setIsValid(kTRUE);
}

//--------------------------------------------------------------
void BdkPdfCBArgus::deleteFrac() {
  // Deletes pointers:
  setIsValid(kFALSE);
  if ( _fracCBShape != 0 )   delete _fracCBShape;
}

//--------------------------------------------------------------
void BdkPdfCBArgus::addComponents() {
  // Adds the wrappers to the composite's list:
  addComponent(_cbShape);
  addComponent(_argus);
}

//--------------------------------------------------------------
void BdkPdfCBArgus::linkParameters( BdkPdfCBArgus & pdf2)
{
  // Copies RRV pointers:
  _cbShape.linkParameters(pdf2.cbShape());
  _argus.linkParameters(pdf2.argus());
  _fracCBShape = pdf2.fracCBShape();
}

void BdkPdfCBArgus::linkParameters(RooAbsReal * m0,
                                   RooAbsReal * sigma,
                                   RooAbsReal * alpha,
                                   RooAbsReal * enne,
                                   RooAbsReal * endPoint,
                                   RooAbsReal * exp,
                                   RooAbsReal * fracB)
{
  // Copies RRV pointers:
  _cbShape.linkParameters(m0, sigma, alpha, enne);
  _argus.linkParameters(endPoint, exp);
                                                                                                                                          
  if (0 != fracB) {
    _fracCBShape = fracB;
    setIsValid(kFALSE);
  }
}

