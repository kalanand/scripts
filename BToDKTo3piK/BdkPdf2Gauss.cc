/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdf2Gauss.cc,v 1.2 2006/03/04 02:03:10 fwinkl Exp $
 * Authors:
 *   WW, Wolfgang Walkowiak, UC Santa Cruz, walkowia@scipp.ucsc.edu
 * Description:
 *   Class for a double gaussian pdf with it's own variable definitions.
 *   To be used in conjunction with RooFitTools.
 * History:
 *   20-Mar-2002 WW Created initial version
 *
 * Copyright (C) 2002 University of California and SLAC
 *****************************************************************************/

// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
// 
// Wrapper for double gaussian pdf.
// 

#include "TString.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooArgList.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooRealConstant.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooAddModel.hh"
#include "RooFitModels/RooGaussian.hh"

#include "BToDKTo3piK/BdkPdf2Gauss.hh"

ClassImp(BdkPdf2Gauss)
    
BdkPdf2Gauss::BdkPdf2Gauss() {
}

BdkPdf2Gauss::BdkPdf2Gauss(const char * theName,
			   const char * theDesc,
			   RooRealVar & dt) {
//
// Constructor
//
    init(theName, theDesc, dt);
}

BdkPdf2Gauss::~BdkPdf2Gauss() {
//
// Destructor
//
}

void BdkPdf2Gauss::init(const char * theName,
			const char * theDesc,
			RooRealVar & dt) {
//
// Initialization
//
    baseInit(theName, theDesc);
    setDependent(dt);
    initParameters();
}

void BdkPdf2Gauss::setDependent(RooRealVar & dt) {
//
// set dependent
//
    setIsValid(kFALSE);
    _dt = &dt;
}

void BdkPdf2Gauss::initParameters() {
//
// Initialize parameters 
// The parameter variables are actually created here.
//
    setIsValid(kFALSE);
    
    // now define the RooVars
    _f1 = new RooRealVar(TString(GetName())+".f1",
			 "First Gaussian Fraction",
			 0.5,0.,1.);
    _b1 = new RooRealVar(TString(GetName())+".b1",
			 "First Gaussian Bias",
			 0.);
    _s1 = new RooRealVar(TString(GetName())+".s1",
			 "First Gaussian Scale Factor",
			 1.);
    _b2 = new RooRealVar(TString(GetName())+".b2",
			 "Second Gaussian Bias",
			 0.);
    _s2 = new RooRealVar(TString(GetName())+".s2",
			 "Second Gaussian Scale Factor",
			 2.);
  
    // set reasonable errors as first steps
    _f1->setConstant(kFALSE);
    _b1->setConstant(kFALSE);
    _s1->setConstant(kFALSE);
    _b2->setConstant(kFALSE);
    _s2->setConstant(kFALSE);

    _f1->setError(0.05);
    _b1->setError(0.05);
    _s1->setError(0.05);
    _b2->setError(0.05);
    _s2->setError(0.05);
}

void BdkPdf2Gauss::createPdf() {
//
// Build the pdf.
//    
    
    _g1 = new RooGaussian(TString(GetName())+".g1","First Gaussian", 
			 *_dt, *_b1, *_s1);
    _g2 = new RooGaussian(TString(GetName())+".g2","Second Gaussian", 
			 *_dt, *_b2, *_s2);

    _thePdf = new RooAddModel(TString(GetName())+".Pdf",
			      TString(GetTitle())+" Pdf",
			      RooArgList(*_g1,*_g2),
			      RooArgList(*_f1));
    
    setIsValid(kTRUE);
}



void BdkPdf2Gauss::linkResolution(const BdkPdf2Gauss & pdf2) {
//
// Link all resolution parameters of the two pdfs.
//
    setIsValid(kFALSE);
    linkResolution(pdf2._f1 ,
		   pdf2._b1 , pdf2._s1,
		   pdf2._b2 , pdf2._s2);
}

void BdkPdf2Gauss::linkResolution(RooRealVar *f1,
				  RooRealVar *b1,  RooRealVar *s1,
				  RooRealVar *b2,  RooRealVar *s2) {
//
// Link the resolution parameters which are not set to 0.
//
    setIsValid(kFALSE);
    if ( f1 != 0 ) {_f1 = f1;}
    if ( b1 != 0 ) {_b1 = b1;}
    if ( s1 != 0 ) {_s1 = s1;}
    if ( b2 != 0 ) {_b2 = b2;}
    if ( s2 != 0 ) {_s2 = s2;}
}

RooArgSet BdkPdf2Gauss::dependents() {
//
// Get dependent.
//    
    return RooArgSet(*_dt);
}
