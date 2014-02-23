/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfDDalitzInc.cc,v 1.16 2006/07/11 21:24:31 fwinkl Exp $
 * Description:
 *   Signal Dalitz PDF wrapper
 * History:
 *   18 Oct 2005, created, Abi soffer
 *
 * Copyright (C) 2005 Colorado State University and SLAC
 *****************************************************************************/
// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
// 
// Wrapper for an incoherent sum of BW's
// 

#include "TString.h"

#include "RooFitCore/RooAbsPdf.hh"  
#include "RooFitCore/RooAddPdf.hh"  
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooRealProxy.hh"

#include "BToDKTo3piK/BdkPdfDDalitzInc.hh"
#include "BToDKTo3piK/BdkDDalitzAmp.hh"
#include "BToDKTo3piK/BdkDDalitz.hh"


ClassImp(BdkPdfDDalitzInc)


// Constructors:
BdkPdfDDalitzInc::BdkPdfDDalitzInc() : 
  _externalResPdfs(kFALSE)
{
  setIsValid(kFALSE);
}


BdkPdfDDalitzInc::BdkPdfDDalitzInc(const char * theName, const char * theDesc,
				   RooAbsReal & m12, RooAbsReal & m13, 
				   BdkDalitzBase::Flavor flavor, 
                                   Int_t spinResComp,
                                   BdkPdfDDalitzInc* externalAmps) : 
  _externalResPdfs(kFALSE)
{
  init(theName, theDesc, m12, m13, flavor, spinResComp);
}

// destructor:
BdkPdfDDalitzInc::~BdkPdfDDalitzInc() {
}

// initializer:  
void BdkPdfDDalitzInc::init(const char * theName, const char * theDesc,
			    RooAbsReal & m12, RooAbsReal & m13, 
			    BdkDalitzBase::Flavor flav,
                            Int_t spinResComp,
                            BdkPdfDDalitzInc* externalAmps) {
  
  // base class initialization:
  BdkPdfDalitzBase::init(theName, theDesc, m12, m13, flav);
  
  _spinResComp = spinResComp;  
  
  // initialize the parameters and the nonres PDF:
  initParameters();
  
  // initialize the owned res PDFs:
  _ownedRhoP.init(TString(GetName())+".rho+",
		  TString(GetTitle())+" rho+ Pdf", 
		  *_m12, *_m13, 
		  flavor(), 
		  externalAmps ? externalAmps->rhoP()->pdfType()->dalitzAmp() : 0,
		  BdkDDalitzAmp::RHOP,
		  _spinResComp);
  
  _ownedRhoM.init(TString(GetName())+".rho-",
		  TString(GetTitle())+" rho- Pdf", 
		  *_m12, *_m13, 
		  flavor(), 
		  externalAmps ? externalAmps->rhoM()->pdfType()->dalitzAmp() : 0,
		  BdkDDalitzAmp::RHOM,
		  _spinResComp);
  
  _ownedRho0.init(TString(GetName())+".rho0",
		  TString(GetTitle())+" rho0 Pdf", 
		  *_m12, *_m13, 
		  flavor(), 
		  externalAmps ? externalAmps->rho0()->pdfType()->dalitzAmp() : 0,
		  BdkDDalitzAmp::RHO0,
		  _spinResComp);
  
  // Point to the owned PDFs as the ones we're using:
  _externalResPdfs = kFALSE;
  _rhoP = &_ownedRhoP;
  _rhoM = &_ownedRhoM;
  _rho0 = &_ownedRho0;
}


// initialize with external resonant PDFs:
void BdkPdfDDalitzInc::init(const char * theName, const char * theDesc,
			    BdkPdfDDalitzInc & resonanceSource) {
  // base class init using info from resonanceSource:
  BdkPdfDalitzBase::init(theName, theDesc, 
			 *(resonanceSource.m12()),
			 *(resonanceSource.m13()),
			 resonanceSource.flavor());
  
  // Point to the resonanceSource's resonance PDFs:
  _externalResPdfs = kTRUE;
  _rhoP = resonanceSource.rhoP();
  _rhoM = resonanceSource.rhoM();
  _rho0 = resonanceSource.rho0();

  // Initialize the parameters and the nonres PDF:
  initParameters();
}


void BdkPdfDDalitzInc::initParameters() {
  // initialize the fractions:
  _fracRhoP = new RooRealVar(TString(GetName()) + ".fracRho+",
			     TString(GetTitle()) + " fracRho+",
			     0.3, 0, 1);
			     
  _fracRho0 = new RooRealVar(TString(GetName()) + ".fracRho0",
			     TString(GetTitle()) + " fracRho0",
			     0.3, 0, 1);
			     
  _fracNonres = new RooRealVar(TString(GetName()) + ".fracNonres",
			     TString(GetTitle()) + " fracNonres",
			     0.3, 0, 1);
			     
  _nonres.init(TString(GetName())+".nonres_pdf",
	       TString(GetTitle())+" nonres Pdf", 
	       *_m12, *_m13, 
	       flavor());

  // Don't need the constant term in the polynomial
  _nonres.c(0)->setVal(1);
  _nonres.c(0)->setConstant();

  setIsValid(kFALSE);
}


/// normalize all rho components
void BdkPdfDDalitzInc::calNorm (int events, Bool_t enforceNormalizeAll)
{
  _rhoP->pdfType()->dalitzAmp()->calNorm(events, enforceNormalizeAll);
  _rhoM->pdfType()->dalitzAmp()->calNorm(events, enforceNormalizeAll);
  _rho0->pdfType()->dalitzAmp()->calNorm(events, enforceNormalizeAll);
}


// Build the PDF:
void BdkPdfDDalitzInc::createPdf() {
  RooAbsPdf* sum1 = new RooAddPdf(TString(GetName()) + "rho+-_pdf",
                                  TString(GetTitle()) + " rho+- Pdf",
                                  *(_rhoP->getPdf()), 
				  *(_rhoM->getPdf()), 
				  *_fracRhoP);

  RooAbsPdf* sum2 = new RooAddPdf(TString(GetName()) + "rho+-0_pdf",
                                  TString(GetTitle()) + "rho+-0 Pdf",
                                  *(_rho0->getPdf()), 
				  *sum1, 
				  *_fracRho0);

  RooAbsPdf* thePdf = new RooAddPdf(TString(GetName()) + ".pdf",
                                    TString(GetTitle()) + " pdf",
                                    *(_nonres.getPdf()), 
				    *sum2, 
				    *_fracNonres);
 
  setPdf(*thePdf);
  _rhoP->setEfficiencyFunc(_effFunc);
  _rhoM->setEfficiencyFunc(_effFunc);
  _rho0->setEfficiencyFunc(_effFunc);

  setIsValid(kTRUE);
}










