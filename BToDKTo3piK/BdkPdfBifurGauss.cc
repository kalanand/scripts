/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfBifurGauss.cc,v 1.4 2006/03/01 22:23:38 fwinkl Exp $
 * Authors:
 *   WW, Wolfgang Walkowiak, UC Santa Cruz, walkowia@scipp.ucsc.edu
 * Description:
 *   Class for a BifurGauss pdf with it's own variable definition
 *   To be used in conjunction with RooFitCore/Models.
 * History:
 *   18-Mar-2002 WW Created initial version
 *   28-Mar-2002 WW add dependents(), delete generate()
 *   16-Aug-2002 WW Force numerical integration for RooBifurGauss until 
 *                  normalization problem in RooBifurGauss is solved.
 *
 * Copyright (C) 2002 University of California and SLAC
 *****************************************************************************/

// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
// 
// Wrapper for bifurcated gaussian pdf.
// 

#include <iostream>
using std::cout;
using std::endl;

#include "TString.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitModels/RooBifurGauss.hh"

#include "BToDKTo3piK/BdkPdfBifurGauss.hh"

ClassImp(BdkPdfBifurGauss)

BdkPdfBifurGauss::BdkPdfBifurGauss() {
//
// default constructor
//
}

BdkPdfBifurGauss::BdkPdfBifurGauss(const char * theName,
				   const char * theDesc,
				   RooRealVar & dependent)
{
    init(theName, theDesc, dependent);
}

BdkPdfBifurGauss::~BdkPdfBifurGauss() {

}

BdkPdfBifurGauss::BdkPdfBifurGauss(BdkPdfBifurGauss & pdf) {
//
// Dummy copy constructor
//
    cout << "BdkPdfBifurGauss::BdkPdfBifurGauss -- dummy copy constructor"
	 << endl;
}

void BdkPdfBifurGauss::init(const char * theName,
			    const char * theDesc,
			    RooRealVar & dependent) {
    baseInit(theName, theDesc);
    setDependent(dependent);
    initParameters();
}

void BdkPdfBifurGauss::setDependent(RooRealVar & dependent) {
    setIsValid(kFALSE);
    _dependent = &dependent;
}
  
void BdkPdfBifurGauss::initParameters() {
    setIsValid(kFALSE);
    
    // Now define the RooVars. The default values are those appropriate
    // for missing mass: 
    _mean = new RooRealVar(TString(GetName())+".mean" ,"mean", 1);
    _sigl = new RooRealVar(TString(GetName())+".sigl" ,"sig left" , 1);
    _sigr = new RooRealVar(TString(GetName())+".sigr" ,"sig right", 1);
    
    // Initial errors should be much less than an MeV to be useful for mRec:
    _mean->setConstant(kFALSE);
    _sigl->setConstant(kFALSE);
    _sigr->setConstant(kFALSE);
    _mean->setError(0.1);
    _sigl->setError(0.1);
    _sigr->setError(0.1);    
}
 
void BdkPdfBifurGauss::createPdf() {
  
    RooBifurGauss *f = 
	new RooBifurGauss(TString(GetName())+".pdf","Bifur Gaussian",
			  *_dependent,*_mean,*_sigl,*_sigr);

    // 
    // ww 08/16/2002
    // Force numerical integration for RooBifurGauss until 
    // normalization problem in RooBifurGauss is solved.
    //
    // fw 03/01/2006
    // Seems to be fixed now
    //    f->forceNumInt(kTRUE);
    //    cout << "BdkPdfBifurGauss::createPdf: forcing numerical integration "
    //	       << "for RooBifurGauss! (temp. fix)" << endl;
    //
    _thePdf = f;
    setIsValid(kTRUE);
}

RooArgSet BdkPdfBifurGauss::dependents() {
//
// Get dependents.
//
    
    return RooArgSet(*_dependent);
}

void BdkPdfBifurGauss::linkParameters(const BdkPdfBifurGauss & pdf2) {
    setIsValid(kFALSE);
    linkParameters(pdf2._mean , pdf2._sigl, pdf2._sigr);
}

void BdkPdfBifurGauss::linkParameters(RooRealVar *mean,
				      RooRealVar *sigl, RooRealVar *sigr) {
    setIsValid(kFALSE);
    if ( mean != 0 ) {_mean = mean;}
    if ( sigl != 0 ) {_sigl = sigl;} 
    if ( sigr != 0 ) {_sigr = sigr;}
}




