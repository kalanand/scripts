/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfGaussArgus.cc,v 1.1 2005/10/09 22:50:37 abi Exp $
 * Authors:
 *   Abi Soffer, Colorado State University, abi@slac.stanford.edu
 * Description:
 *   Class for a Gaussian + ARGUS pdf with it's own variable definition
 *   To be used in conjunction with RooFitCore/Models.
 * History:
 *   6-Mar-2004 abi Created initial version
 *
 * Copyright (C) 2004 Colorado State University and SLAC
 *****************************************************************************/

// -- CLASS DESCRIPTION [IHFPDFWRAPPER] --
// 
// Wrapper for gaussian plus argus pdf.
// 

#include <iostream>
using std::cout;
using std::endl;

#include "TString.h"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooArgSet.hh"

#include "BToDKTo3piK/BdkPdfGaussArgus.hh"

ClassImp(BdkPdfGaussArgus)

//--------------------------------------------------------------
BdkPdfGaussArgus::BdkPdfGaussArgus() : _fracGauss(0) {
  // Default Constructor:
  addComponents();
}

//--------------------------------------------------------------
BdkPdfGaussArgus::BdkPdfGaussArgus(const char * theName,
				   const char * theDesc,
				   RooRealVar & dependent)  : _fracGauss(0) {
  // Initializing Constructor:
  addComponents();
  init(theName, theDesc, dependent);
}

//--------------------------------------------------------------
BdkPdfGaussArgus::~BdkPdfGaussArgus() {
  // Destructor:
  deleteFrac();
}

//--------------------------------------------------------------
BdkPdfGaussArgus::BdkPdfGaussArgus(const BdkPdfGaussArgus & pdf) {
// Dummy copy constructor:
    cout << "BdkPdfGaussArgus::BdkPdfGaussArgus -- dummy copy constructor"
	 << endl;
}


//--------------------------------------------------------------
void BdkPdfGaussArgus::init(const char * theName,
			    const char * theDesc,
			    RooRealVar & dependent) {
  // Initializer:
  _fracGauss = 0;

  SetName(theName);
  SetTitle(theDesc);

  baseInit(theName, theDesc);
  setDependent(dependent);
  initParameters();
}

//--------------------------------------------------------------
void BdkPdfGaussArgus::setDependent(RooRealVar & dependent) {
  // Sets the dependent variable:
  _dependent = &dependent;
  _gauss.setDependent(dependent);
  _argus.setDependent(dependent);

  setIsValid(kFALSE);
}

//--------------------------------------------------------------
void BdkPdfGaussArgus::initParameters() {
  // Initializes the contained RRV's and PDF wrappers:
  _fracGauss = new RooRealVar(GetName() + TString(".fracGauss"),
			      GetTitle () + TString(" fracGauss"),
			      0.8, 0, 1);

  _fracGauss->setError(0.02);

  _gauss.init(GetName() + TString(".gauss"), 
	      GetTitle() + TString(" Gauss"), 
	      *_dependent);

  _argus.init(GetName() + TString(".argus"), 
	      GetTitle() + TString(" Argus"), 
	      *_dependent);
}

//--------------------------------------------------------------
void BdkPdfGaussArgus::createPdf() {
  // Builds the PDF. The PDFs of the contained wrappers are automatically built
  _thePdf = new RooAddPdf(GetName() + TString(".thePdf"), 
			  GetTitle() + TString(" thePdf"), 
			  *_gauss.getPdf(), 
			  *_argus.getPdf(),
			  *_fracGauss);

  setIsValid(kTRUE);
}

//--------------------------------------------------------------
void BdkPdfGaussArgus::deleteFrac() {
  // Deletes pointers:
  setIsValid(kFALSE);
  if ( _fracGauss != 0 )   delete _fracGauss;
}

//--------------------------------------------------------------
void BdkPdfGaussArgus::addComponents() {
  // Adds the wrappers to the composite's list:
  addComponent(_gauss);
  addComponent(_argus);
}

//--------------------------------------------------------------
void BdkPdfGaussArgus::linkParameters(RooRealVar * mean,
				      RooRealVar * sigma,
				      RooRealVar * endPoint, 
				      RooRealVar * exp,
				      RooRealVar * fracB) {
  // Copies RRV pointers:
  _gauss.linkResolution(mean, sigma);
  _argus.linkParameters(endPoint, exp);

  if (0 != fracB) {
    _fracGauss = fracB;
    setIsValid(kFALSE);
  }
}


