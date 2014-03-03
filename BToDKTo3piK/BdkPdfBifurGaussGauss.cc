/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfBifurGaussGauss.cc,v 1.1 2006/03/01 21:39:55 fwinkl Exp $
 * Authors:
 *   Abi Soffer, Colorado State University, abi@slac.stanford.edu
 * Description:
 *   Class for a Gaussian + Bifurcated Gaussian pdf with it's own 
 *   variable definition to be used in conjunction with RooFitCore/Models.
 * History:
 *   1-Mar-2006 fwinkl Created initial version
 *
 * Copyright (C) 2006 Colorado State University and SLAC
 *****************************************************************************/
// -- CLASS DESCRIPTION [IHFPDFWRAPPER] --
// 
// Wrapper for gaussian plus bifurcated gaussian pdf
// 

#include <iostream>
using std::cout;
using std::endl;

#include "TString.h"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooArgSet.hh"

#include "BToDKTo3piK/BdkPdfBifurGaussGauss.hh"

ClassImp(BdkPdfBifurGaussGauss)

//--------------------------------------------------------------
BdkPdfBifurGaussGauss::BdkPdfBifurGaussGauss() : _fracGauss(0) {
  // Default Constructor:
  addComponents();
}

//--------------------------------------------------------------
BdkPdfBifurGaussGauss::BdkPdfBifurGaussGauss(const char * theName,
				   const char * theDesc,
				   RooRealVar & dependent)  : _fracGauss(0) {
  // Initializing Constructor:
  addComponents();
  init(theName, theDesc, dependent);
}

//--------------------------------------------------------------
BdkPdfBifurGaussGauss::~BdkPdfBifurGaussGauss() {
  // Destructor:
  deleteFrac();
}

//--------------------------------------------------------------
BdkPdfBifurGaussGauss::BdkPdfBifurGaussGauss(const BdkPdfBifurGaussGauss & pdf) {
// Dummy copy constructor:
  cout << "BdkPdfBifurGaussGauss::BdkPdfBifurGaussGauss -- dummy copy constructor"
       << endl;
}


//--------------------------------------------------------------
void BdkPdfBifurGaussGauss::init(const char * theName,
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
void BdkPdfBifurGaussGauss::setDependent(RooRealVar & dependent) {
  // Sets the dependent variable:
  _dependent = &dependent;
  _gauss.setDependent(dependent);
  _bifur.setDependent(dependent);

  setIsValid(kFALSE);
}

//--------------------------------------------------------------
void BdkPdfBifurGaussGauss::initParameters() {
  // Initializes the contained RRV's and PDF wrappers:
  _fracGauss = new RooRealVar(GetName() + TString(".fracGauss"),
			      GetTitle () + TString(" fracGauss"),
			      0.8, 0, 1);

  _fracGauss->setError(0.02);

  _gauss.init(GetName() + TString(".gauss"), 
	      GetTitle() + TString(" Gauss"), 
	      *_dependent);

  _bifur.init(GetName() + TString(".bifur"), 
	      GetTitle() + TString(" Bifurcated Gauss"), 
	      *_dependent);
}

//--------------------------------------------------------------
void BdkPdfBifurGaussGauss::createPdf() {
  // Builds the PDF. The PDFs of the contained wrappers are automatically built
  _thePdf = new RooAddPdf(GetName() + TString(".thePdf"), 
			  GetTitle() + TString(" thePdf"), 
			  *_gauss.getPdf(), 
			  *_bifur.getPdf(),
			  *_fracGauss);

  setIsValid(kTRUE);
}

//--------------------------------------------------------------
void BdkPdfBifurGaussGauss::deleteFrac() {
  // Deletes pointers:
  setIsValid(kFALSE);
  if ( _fracGauss != 0 )   delete _fracGauss;
}

//--------------------------------------------------------------
void BdkPdfBifurGaussGauss::addComponents() {
  // Adds the wrappers to the composite's list:
  addComponent(_gauss);
  addComponent(_bifur);
}

//--------------------------------------------------------------
void BdkPdfBifurGaussGauss::linkParameters(RooRealVar * b1,
                                           RooRealVar * s1,
                                           RooRealVar * b2, 
                                           RooRealVar * s2l,
                                           RooRealVar * s2r,
                                           RooRealVar * fracGauss) {
  // Copies RRV pointers:
  _gauss.linkResolution(b1, s1);
  _bifur.linkParameters(b2, s2l, s2r);

  if (0 != fracGauss) {
    _fracGauss = fracGauss;
    setIsValid(kFALSE);
  }
}


