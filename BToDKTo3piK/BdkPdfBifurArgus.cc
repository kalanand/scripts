/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfBifurArgus.cc,v 1.1 2005/10/09 22:50:36 abi Exp $
 * Authors:
 *   Abi Soffer, Colorado State University, abi@slac.stanford.edu
 * Description:
 *   Class for a bifur + ARGUS pdf with it's own variable definition
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
// Wrapper for bifurcated gaussian plus argus pdf.
// 

#include <iostream>
using std::cout;
using std::endl;

#include "TString.h"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooArgSet.hh"

#include "BToDKTo3piK/BdkPdfBifurArgus.hh"

ClassImp(BdkPdfBifurArgus)

//--------------------------------------------------------------
BdkPdfBifurArgus::BdkPdfBifurArgus() : _fracBifur(0) {
  // Default Constructor:
  addComponents();
}

//--------------------------------------------------------------
BdkPdfBifurArgus::BdkPdfBifurArgus(const char * theName,
				   const char * theDesc,
				   RooRealVar & dependent)  : _fracBifur(0) {
  // Initializing Constructor:
  addComponents();
  init(theName, theDesc, dependent);
}

//--------------------------------------------------------------
BdkPdfBifurArgus::~BdkPdfBifurArgus() {
  // Destructor:
  deleteFrac();
}

//--------------------------------------------------------------
BdkPdfBifurArgus::BdkPdfBifurArgus(const BdkPdfBifurArgus & pdf) {
// Dummy copy constructor:
    cout << "BdkPdfBifurArgus::BdkPdfBifurArgus -- dummy copy constructor"
	 << endl;
}


//--------------------------------------------------------------
void BdkPdfBifurArgus::init(const char * theName,
			    const char * theDesc,
			    RooRealVar & dependent) {
  // Initializer:
  _fracBifur = 0;

  SetName(theName);
  SetTitle(theDesc);

  baseInit(theName, theDesc);
  setDependent(dependent);
  initParameters();
}

//--------------------------------------------------------------
void BdkPdfBifurArgus::setDependent(RooRealVar & dependent) {
  // Sets the dependent variable:
  _dependent = &dependent;
  _bifur.setDependent(dependent);
  _argus.setDependent(dependent);

  setIsValid(kFALSE);
}

//--------------------------------------------------------------
void BdkPdfBifurArgus::initParameters() {
  // Initializes the contained RRV's and PDF wrappers:
  _fracBifur = new RooRealVar(GetName() + TString(".fracBifur"),
			      GetTitle () + TString(" fracBifur"),
			      0.8, 0, 1);

  _fracBifur->setError(0.02);

  _bifur.init(GetName() + TString(".bifur"), 
	      GetTitle() + TString(" BifurGauss"), 
	      *_dependent);

  _argus.init(GetName() + TString(".argus"), 
	      GetTitle() + TString(" Argus"), 
	      *_dependent);
}

//--------------------------------------------------------------
void BdkPdfBifurArgus::createPdf() {
  // Builds the PDF. The PDFs of the contained wrappers are automatically built
  _thePdf = new RooAddPdf(GetName() + TString(".thePdf"), 
			  GetTitle() + TString(" thePdf"), 
			  *_bifur.getPdf(), 
			  *_argus.getPdf(),
			  *_fracBifur);

  setIsValid(kTRUE);
}

//--------------------------------------------------------------
void BdkPdfBifurArgus::deleteFrac() {
  // Deletes pointers:
  setIsValid(kFALSE);
  if ( _fracBifur != 0 )   delete _fracBifur;
}

//--------------------------------------------------------------
void BdkPdfBifurArgus::addComponents() {
  // Adds the wrappers to the composite's list:
  addComponent(_bifur);
  addComponent(_argus);
}

//--------------------------------------------------------------
void BdkPdfBifurArgus::linkParameters(RooRealVar * mean,
				      RooRealVar * sigl, 
				      RooRealVar * sigr,
				      RooRealVar * endPoint, 
				      RooRealVar * exp,
				      RooRealVar * fracB) {
  // Copies RRV pointers:
  _bifur.linkParameters(mean, sigl, sigr);
  _argus.linkParameters(endPoint, exp);

  if (0 != fracB) {
    _fracBifur = fracB;
    setIsValid(kFALSE);
  }
}


