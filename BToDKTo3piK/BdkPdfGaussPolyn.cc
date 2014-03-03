/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfGaussPolyn.cc,v 1.1 2005/10/09 22:50:37 abi Exp $
 * Authors:
 *   Abi Soffer, Colorado State University, abi@slac.stanford.edu
 * Description:
 *   Class for a Gaussian + POLYN pdf with it's own variable definition
 *   To be used in conjunction with RooFitCore/Models.
 * History:
 *   6-Mar-2004 abi Created initial version
 *
 * Copyright (C) 2004 Colorado State University and SLAC
 *****************************************************************************/

// -- CLASS DESCRIPTION [IHFPDFWRAPPER] --
// 
// Wrapper for gaussian plus polyn pdf.
// 

#include <iostream>
using std::cout;
using std::endl;

#include "TString.h"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooArgSet.hh"

#include "BToDKTo3piK/BdkPdfGaussPolyn.hh"

ClassImp(BdkPdfGaussPolyn)

//--------------------------------------------------------------
BdkPdfGaussPolyn::BdkPdfGaussPolyn() : _fracGauss(0) {
  // Default Constructor:
  addComponents();
}

//--------------------------------------------------------------
BdkPdfGaussPolyn::BdkPdfGaussPolyn(const char * theName,
				   const char * theDesc,
				   RooRealVar & dependent,
				   Int_t lowestOrder,
				   Int_t highestOrder)  : _fracGauss(0) {
  // Initializing Constructor:
  addComponents();
  init(theName, theDesc, dependent,
       lowestOrder, highestOrder);
}

//--------------------------------------------------------------
BdkPdfGaussPolyn::~BdkPdfGaussPolyn() {
  // Destructor:
  deleteFrac();
}

//--------------------------------------------------------------
BdkPdfGaussPolyn::BdkPdfGaussPolyn(const BdkPdfGaussPolyn & pdf) {
// Dummy copy constructor:
    cout << "BdkPdfGaussPolyn::BdkPdfGaussPolyn -- dummy copy constructor"
	 << endl;
}


//--------------------------------------------------------------
void BdkPdfGaussPolyn::init(const char * theName,
			    const char * theDesc,
			    RooRealVar & dependent,
			    Int_t lowestOrder,
			    Int_t highestOrder) {
  // Initializer:
  _fracGauss = 0;

  SetName(theName);
  SetTitle(theDesc);

  baseInit(theName, theDesc);
  setDependent(dependent);

  _lowestOrder = lowestOrder;
  _highestOrder = highestOrder;
  initParameters();
}

//--------------------------------------------------------------
void BdkPdfGaussPolyn::setDependent(RooRealVar & dependent) {
  // Sets the dependent variable:
  _dependent = &dependent;
  _gauss.setDependent(dependent);
  _polyn.setDependent(dependent);

  setIsValid(kFALSE);
}

//--------------------------------------------------------------
void BdkPdfGaussPolyn::initParameters() {
  // Initializes the contained RRV's and PDF wrappers:
  _fracGauss = new RooRealVar(GetName() + TString(".fracGauss"),
			      GetTitle () + TString(" fracGauss"),
			      0.5, 0, 1);

  _fracGauss->setError(0.01);

  _gauss.init(GetName() + TString(".gauss"), 
	      GetTitle() + TString(" Gauss"), 
	      *_dependent);

  _polyn.init(GetName() + TString(".polyn"), 
	      GetTitle() + TString(" Polyn"), 
	      *_dependent, _lowestOrder, 
              _highestOrder);
}

//--------------------------------------------------------------
void BdkPdfGaussPolyn::createPdf() {
  // Builds the PDF. The PDFs of the contained wrappers are automatically built
  _thePdf = new RooAddPdf(GetName() + TString(".thePdf"), 
			  GetTitle() + TString(" thePdf"), 
			  *_gauss.getPdf(), 
			  *_polyn.getPdf(),
			  *_fracGauss);

  setIsValid(kTRUE);
}

//--------------------------------------------------------------
void BdkPdfGaussPolyn::deleteFrac() {
  // Deletes pointers:
  setIsValid(kFALSE);
  if ( _fracGauss != 0 )   delete _fracGauss;
}

//--------------------------------------------------------------
void BdkPdfGaussPolyn::addComponents() {
  // Adds the wrappers to the composite's list:
  addComponent(_gauss);
  addComponent(_polyn);
}

//--------------------------------------------------------------
void BdkPdfGaussPolyn::linkParameters(RooRealVar * mean,
				      RooRealVar * sigma,
				      Int_t lowestOrder,
				      Int_t highestOrder,	 
				      RooRealVar * fracB) {
  // Copies RRV pointers:
  _gauss.linkResolution(mean, sigma);
  _polyn.linkParameters(lowestOrder, highestOrder);

  if (0 != fracB) {
    _fracGauss = fracB;
    setIsValid(kFALSE);
  }
}

//---------------------------------------------------
void BdkPdfGaussPolyn::linkParameters( BdkPdfGaussPolyn & pdf2) {
	_gauss.linkResolution(pdf2.gauss());
        _polyn.linkParameters(pdf2.polyn());
	_fracGauss = pdf2.fracGauss();
} 
