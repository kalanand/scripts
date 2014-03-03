/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Description:
 *   Signal Dalitz PDF wrapper
 * History:
 *   18 Oct 2005, created, Abi soffer
 *   31 May 2006, adapted for KKpi0 DP, Kalanand Mishra
 *
 *****************************************************************************/


#include "TString.h"

#include "RooFitCore/RooRealVar.hh"
#include "BToDKTo3piK/KKPDalitz.hh"
#include "BToDKTo3piK/PdfKKPDalitz.hh"


ClassImp(PdfKKPDalitz)


// Constructors:
PdfKKPDalitz::PdfKKPDalitz() {
  setIsValid(kFALSE);
}


PdfKKPDalitz::PdfKKPDalitz(const char * theName, const char * theDesc,
			     RooAbsReal & m12, RooAbsReal & m13,
			     BdkDalitzBase::Flavor flavor, 
			     KKPDalitzAmp * externalAmp,
			     Int_t comps,
			     Int_t spinResComp) {
  init(theName, theDesc, m12, m13, flavor, externalAmp, comps, spinResComp);
}
  
// destructor:
PdfKKPDalitz::~PdfKKPDalitz() {}

// initializer:  
void PdfKKPDalitz::init(const char * theName, const char * theDesc,
			 RooAbsReal & m12, RooAbsReal & m13, 
			 BdkDalitzBase::Flavor flavor,
			 KKPDalitzAmp * externalAmp,
			 Int_t comps,
			 Int_t spinResComp) {

  // base class initialization:
  BdkPdfDalitzBase::init(theName, theDesc, m12, m13, flavor);
  setComponentsBit(comps);
  _externalAmp = externalAmp;
  _spinResComp = spinResComp;
  setIsValid(kFALSE);
}

void PdfKKPDalitz::setComponentsBit(Int_t comps) {
  _componentsBit = comps;
}


const BdkDalitzEff* PdfKKPDalitz::efficiencyFunc() const
{
  return ((KKPDalitz*)((PdfKKPDalitz*)this)->getPdf())->efficiencyFunc();
}

void PdfKKPDalitz::setEfficiencyFunc(const BdkDalitzEff *f)
{
  ((KKPDalitz*)((PdfKKPDalitz*)this)->getPdf())->setEfficiencyFunc(f);
}


// Build the PDF:
void PdfKKPDalitz::createPdf() {
  RooAbsPdf * thePdf = new KKPDalitz(TString(GetName())+".pdf",
				      TString(GetTitle())+" Pdf", 
				      *_m12, *_m13, 
				      flavor(),
				      _externalAmp,
				      _componentsBit,
				      BdkDalitzBase::KKP0,
				      _spinResComp);
  
  setPdf(*thePdf);
  setIsValid(kTRUE);
}

RooArgSet PdfKKPDalitz::fitFractions() const {
  KKPDalitz * pdf = ((PdfKKPDalitz*)this)->pdfType(); // cast out const
  return pdf->fitFractions();
}

// Send verbosity flag to dalitzAmp:
void PdfKKPDalitz::setVerbose(const char * val) {
  BdkPdfAbsBase::setVerbose(val); // do base class action
  if (verbose().Contains("+")) {  // set to amp
    pdfType()->dalitzAmp()->setVerbose(val);
  }
}








