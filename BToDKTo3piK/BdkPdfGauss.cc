/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfGauss.cc,v 1.4 2006/03/02 01:06:03 fwinkl Exp $
 * Authors:
 *   WW, Wolfgang Walkowiak, UC Santa Cruz, walkowia@scipp.ucsc.edu
 * Description:
 *   Class for a gaussian pdf with it's own variable definitions.
 *   To be used in conjunction with RooFitTools.
 * History:
 *   20-Mar-2002 WW Created initial version
 *
 * Copyright (C) 2002 University of California and SLAC
 *****************************************************************************/
// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
// 
// Wrapper for gaussian pdf.
// 

#include "TString.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitModels/RooGaussian.hh"

#include "BToDKTo3piK/BdkPdfGauss.hh"

ClassImp(BdkPdfGauss)
    
BdkPdfGauss::BdkPdfGauss() {
}

BdkPdfGauss::BdkPdfGauss(const char * theName,
			 const char * theDesc,
			 RooRealVar & dt) {
//
// Constructor
//
    init(theName, theDesc, dt);
}

BdkPdfGauss::~BdkPdfGauss() {
//
// Destructor
//
}

void BdkPdfGauss::init(const char * theName,
		       const char * theDesc,
		       RooRealVar & dt) {
//
// Initialization
//
    baseInit(theName, theDesc);
    setDependent(dt);
    initParameters();
}

void BdkPdfGauss::setDependent(RooRealVar & dt) {
//
// set dependent
//
    setIsValid(kFALSE);
    _dt = &dt;
}

void BdkPdfGauss::initParameters() {
//
// Initialize parameters 
// The parameter variables are actually created here.
//
    setIsValid(kFALSE);
    
    // now define the RooVars
    _b = new RooRealVar(TString(GetName())+".b",
			 TString(GetTitle())+" Bias",
			 0);
    _s = new RooRealVar(TString(GetName())+".s",
			TString(GetTitle())+" Scale Factor",
			 1);

  // set reasonable errors as first steps
    _b->setConstant(kFALSE);
    _s->setConstant(kFALSE);
    _b->setError(0.001);
    _s->setError(0.001);
}

void BdkPdfGauss::createPdf() {
//
// Build the pdf.
//    
  RooAbsPdf * thePdf = new RooGaussian(TString(GetName())+".pdf",
				       TString(GetTitle())+" Pdf", 
				       *_dt, *_b, *_s);
  setPdf(*thePdf);
  setIsValid(kTRUE);
}



void BdkPdfGauss::linkResolution(const BdkPdfGauss & pdf2) {
//
// Link all resolution parameters of the two pdfs.
//
    setIsValid(kFALSE);
    linkResolution(pdf2._b , pdf2._s);
}

void BdkPdfGauss::linkResolution(RooRealVar *b,  RooRealVar *s) {
//
// Link the resolution parameters which are not set to 0.
//
    setIsValid(kFALSE);
    if ( b != 0 ) {_b = b;}
    if ( s != 0 ) {_s = s;}
}

RooArgSet BdkPdfGauss::dependents() {
//
// Get dependent.
//    
    return RooArgSet(*_dt);
}

