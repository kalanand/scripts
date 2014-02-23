/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfCBPolyn.cc,v 1.1 2005/10/09 22:50:36 abi Exp $
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

#include "BToDKTo3piK/BdkPdfCBPolyn.hh"
#include "TString.h"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooArgSet.hh"

ClassImp(BdkPdfCBPolyn)

//--------------------------------------------------------------
BdkPdfCBPolyn::BdkPdfCBPolyn() : _fracCBShape(0) {
  // Default Constructor:
  addComponents();
}

//--------------------------------------------------------------
BdkPdfCBPolyn::BdkPdfCBPolyn(const char * theName,
				   const char * theDesc,
				   RooRealVar & dependent,
				   Int_t lowestOrder,
				   Int_t highestOrder
				   )  : _fracCBShape(0) {
  // Initializing Constructor:
  addComponents();
  init(theName, theDesc, dependent, 
	lowestOrder, highestOrder);
}

//--------------------------------------------------------------
BdkPdfCBPolyn::~BdkPdfCBPolyn() {
  // Destructor:
  deleteFrac();
}

//--------------------------------------------------------------
BdkPdfCBPolyn::BdkPdfCBPolyn(const BdkPdfCBPolyn & pdf) {
// Dummy copy constructor:
    cout << "BdkPdfCBPolyn::BdkPdfCBPolyn -- dummy copy constructor"
	 << endl;
}


//--------------------------------------------------------------
void BdkPdfCBPolyn::init(const char * theName,
			    const char * theDesc,
			    RooRealVar & dependent,
			    Int_t lowestOrder,
			    Int_t highestOrder) {
  // Initializer:
  _fracCBShape = 0;

  SetName(theName);
  SetTitle(theDesc);

  baseInit(theName, theDesc);
  setDependent(dependent);

  _lowestOrder = lowestOrder;
  _highestOrder = highestOrder;
  initParameters();
}

//--------------------------------------------------------------
void BdkPdfCBPolyn::setDependent(RooRealVar & dependent) {
  // Sets the dependent variable:
  _dependent = &dependent;
  _cbShape.setDependent(dependent);
  _polyn.setDependent(dependent);

  setIsValid(kFALSE);
}

//--------------------------------------------------------------
void BdkPdfCBPolyn::initParameters() {
  // Initializes the contained RRV's and PDF wrappers:
  _fracCBShape = new RooRealVar(GetName() + TString(".fracCBShape"),
			      GetTitle () + TString(" fracCBShape"),
			      0.5, 0, 1);

  _fracCBShape->setError(0.01);

  _cbShape.init(GetName() + TString(".cbShape"), 
	      GetTitle() + TString(" cbShape"), 
	      *_dependent);

  _polyn.init(GetName() + TString(".polyn"), 
	      GetTitle() + TString(" Polyn"), 
	      *_dependent, _lowestOrder, 
              _highestOrder);
}

//--------------------------------------------------------------
void BdkPdfCBPolyn::createPdf() {
  // Builds the PDF. The PDFs of the contained wrappers are automatically built
  _thePdf = new RooAddPdf(GetName() + TString(".cbpolyn"), 
			  GetTitle() + TString(" cbpolyn"), 
			  *_cbShape.getPdf(), 
			  *_polyn.getPdf(),
			  *_fracCBShape);

  setIsValid(kTRUE);
}

//--------------------------------------------------------------
void BdkPdfCBPolyn::deleteFrac() {
  // Deletes pointers:
  setIsValid(kFALSE);
  if ( _fracCBShape != 0 )   delete _fracCBShape;
}

//--------------------------------------------------------------
void BdkPdfCBPolyn::addComponents() {
  // Adds the wrappers to the composite's list:
  addComponent(_cbShape);
  addComponent(_polyn);
}

//--------------------------------------------------------------
void BdkPdfCBPolyn::linkParameters(RooRealVar * m0,
				   RooRealVar * sigma,
 				   RooRealVar * alpha,
				   RooRealVar * enne,
				   Int_t lowestOrder,
				   Int_t highestOrder,	 
				   RooRealVar * fracCBShape) {
  // Copies RRV pointers:
  _cbShape.linkParameters(m0, sigma, alpha, enne);
  _polyn.linkParameters(lowestOrder, highestOrder);

  if (0 != fracCBShape) {
    _fracCBShape = fracCBShape;
    setIsValid(kFALSE);
  }
}


