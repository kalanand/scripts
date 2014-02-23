/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfDDalitz.cc,v 1.14 2007/08/07 18:42:06 kalanand Exp $
 * Description:
 *   Signal Dalitz PDF wrapper
 * History:
 *   18 Oct 2005, created, Abi soffer
 *
 * Copyright (C) 2005 Colorado State University and SLAC
 *****************************************************************************/
// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
// 
// Wrapper for B->DK Dalitz plot PDF
// 

#include "TString.h"

#include "RooFitCore/RooRealVar.hh"
#include "BToDKTo3piK/BdkDDalitz.hh"
#include "BToDKTo3piK/BdkPdfDDalitz.hh"


ClassImp(BdkPdfDDalitz)


// Constructors:
BdkPdfDDalitz::BdkPdfDDalitz() {
  setIsValid(kFALSE);
}


BdkPdfDDalitz::BdkPdfDDalitz(const char * theName, const char * theDesc,
                             RooAbsReal & m12, RooAbsReal & m13,
                             BdkDalitzBase::Flavor flavor, 
                             BdkAbsDDalitzAmp * externalAmp,
                             Int_t comps,
                             Int_t spinResComp,
                             BdkDalitz::Mode DdecMode) {
  init(theName, theDesc, m12, m13, flavor, externalAmp, comps, spinResComp, DdecMode);
}


// destructor:
BdkPdfDDalitz::~BdkPdfDDalitz() {}

// initializer:  
void BdkPdfDDalitz::init(const char * theName, const char * theDesc,
                         RooAbsReal & m12, RooAbsReal & m13, 
                         BdkDalitzBase::Flavor flavor,
                         BdkAbsDDalitzAmp * externalAmp,
                         Int_t comps,
                         Int_t spinResComp,
                         BdkDalitz::Mode DdecMode) {

  // base class initialization:
  BdkPdfDalitzBase::init(theName, theDesc, m12, m13, flavor, DdecMode);
  setComponentsBit(comps);
  _externalAmp = externalAmp;
  _spinResComp = spinResComp;
  setIsValid(kFALSE);
}

void BdkPdfDDalitz::setComponentsBit(Int_t comps) {
  _componentsBit = comps;
}


const BdkDalitzEff* BdkPdfDDalitz::efficiencyFunc() const
{
  return ((BdkDDalitz*)((BdkPdfDDalitz*)this)->getPdf())->efficiencyFunc();
}

void BdkPdfDDalitz::setEfficiencyFunc(const BdkDalitzEff *f)
{
  ((BdkDDalitz*)((BdkPdfDDalitz*)this)->getPdf())->setEfficiencyFunc(f);
}


// Build the PDF:
void BdkPdfDDalitz::createPdf() {
  RooAbsPdf * thePdf = new BdkDDalitz(TString(GetName())+".pdf",
                                      TString(GetTitle())+" Pdf", 
                                      *_m12, *_m13, 
                                      flavor(),
                                      _externalAmp,
                                      _componentsBit,
                                      getDdecMode(),
                                      _spinResComp);
  
  // Save pointer to dalitzAmp for reuse on subsequent createPdf() calls
  _externalAmp = ((BdkDDalitz*)thePdf)->dalitzAmp();

  setPdf(*thePdf);
  ((BdkDDalitz*)thePdf)->setEfficiencyFunc(_effFunc);

  setIsValid(kTRUE);
}

RooArgSet BdkPdfDDalitz::fitFractions() const {
  BdkDDalitz * pdf = ((BdkPdfDDalitz*)this)->pdfType(); // cast out const
  return pdf->fitFractions();
}


RooArgSet BdkPdfDDalitz::BreitWignerNormalizationCoefficients() const {
  BdkDDalitz * pdf = ((BdkPdfDDalitz*)this)->pdfType(); // cast out const
  return pdf->BreitWignerNormalizationCoefficients();
}

RooArgSet BdkPdfDDalitz::IntegralOverIsospin() const {
  BdkDDalitz * pdf = ((BdkPdfDDalitz*)this)->pdfType(); 
  return pdf->IntegralOverIsospin();
}

// Send verbosity flag to dalitzAmp:
void BdkPdfDDalitz::setVerbose(const char * val) {
  BdkPdfAbsBase::setVerbose(val); // do base class action
  if (verbose().Contains("+")) {  // set to amp
    pdfType()->dalitzAmp()->setVerbose(val);
  }
}








