/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 *
 * Package: 
 *
 *    File: $Id: BdkPdfCB2A.cc,v 1.1 2005/10/09 22:50:36 abi Exp $
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

#include "BToDKTo3piK/BdkPdfCB2A.hh"

ClassImp(BdkPdfCB2A)

//--------------------------------------------------------------
BdkPdfCB2A::BdkPdfCB2A() {
  // Default Constructor:
  addComponents();
}

//--------------------------------------------------------------
BdkPdfCB2A::BdkPdfCB2A(const char * theName,
				   const char * theDesc,
				   RooRealVar & dependent) {
  // Initializing Constructor:
  addComponents();
  init(theName, theDesc, dependent);
}

//--------------------------------------------------------------
BdkPdfCB2A::~BdkPdfCB2A() {
  // Destructor:
}

//--------------------------------------------------------------
BdkPdfCB2A::BdkPdfCB2A(const BdkPdfCB2A & pdf) {
// Dummy copy constructor:
    cout << "BdkPdfCB2A::BdkPdfCB2A -- dummy copy constructor"
	 << endl;
}


//--------------------------------------------------------------
void BdkPdfCB2A::init(const char * theName,
			    const char * theDesc,
			    RooRealVar & dependent) {
  // Initializer:

  baseInit(theName, theDesc);
  setDependent(dependent);
  initParameters();
}

//--------------------------------------------------------------
void BdkPdfCB2A::setDependent(RooRealVar & dependent) {
  // Sets the dependent variable:
  _dependent = &dependent;
  _cbargus.setDependent(dependent);
  _argus.setDependent(dependent);

  setIsValid(kFALSE);
}

//--------------------------------------------------------------
void BdkPdfCB2A::initParameters() {
  // Initializes the contained RRV's and PDF wrappers:

  _frac = new RooRealVar(TString(GetName()) + ".CBAfrac",
                          TString(GetTitle()) + ".CBAFrac",
                          0.5, 0, 1);

  _cbargus.init(GetName() + TString(".cbargusPdf"), 
	      GetTitle() + TString(" Argus Crystall Ball"), 
	      *_dependent);

  _argus.init(GetName() + TString(".argus"), 
	      GetTitle() + TString(" Argus"), 
	      *_dependent);
}

//--------------------------------------------------------------
void BdkPdfCB2A::createPdf() {

  // Builds the PDF. The PDFs of the contained wrappers are automatically built
  RooAddPdf * thePdf = new RooAddPdf(GetName() + TString(".thePdf"), 
			  GetTitle() + TString(" thePdf"), 
                          *_cbargus.getPdf(), 
			  *_argus.getPdf(),
			  *_frac);

  setPdf(*thePdf, kTRUE);
  setIsValid(kTRUE);
}

//--------------------------------------------------------------
void BdkPdfCB2A::addComponents() {
  // Adds the wrappers to the composite's list:
  addComponent(_cbargus);
  addComponent(_argus);
}

//--------------------------------------------------------------
void BdkPdfCB2A::linkParameters(BdkPdfCB2A & pdf2)
{
  // Copies RRV pointers:
  _cbargus.linkParameters(pdf2.cbargus());
  _argus.linkParameters(pdf2.argus());
  _frac = pdf2.fracCBA();
  setIsValid(kFALSE);
}

//--------------------------------------------------------------
void BdkPdfCB2A::linkEndPoints() {
  _argus.linkParameters(_cbargus.argus().endPoint());
  setIsValid(kFALSE);
}
