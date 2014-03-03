#include "TString.h"
#include "RooFitCore/RooRealVar.hh"
#include "BToDKTo3piK/Bdk2DpolyDalitz.hh"
#include "BToDKTo3piK/PdfKKP2DpolyDalitz.hh"


ClassImp(PdfKKP2DpolyDalitz)
  

/// Constructors:
PdfKKP2DpolyDalitz::PdfKKP2DpolyDalitz() {
  setIsValid(kFALSE);
}


PdfKKP2DpolyDalitz::PdfKKP2DpolyDalitz(const char * theName, 
				       const char * theDesc, 
				       RooAbsReal & m12, 
				       RooAbsReal & m13,
				       BdkDalitzBase::Flavor flavor,
				       RooRealVar* c[10])
{
  init(theName, theDesc, m12, m13, flavor, c);
}

/// destructor:
PdfKKP2DpolyDalitz::~PdfKKP2DpolyDalitz() {}

/// initializer:  
void PdfKKP2DpolyDalitz::init(const char * theName, const char * theDesc,
			      RooAbsReal & m12, RooAbsReal & m13,
			      BdkDalitzBase::Flavor flavor,
			      RooRealVar* c[10])
{
  
  // base class initialization:
  BdkPdfDalitzBase::init(theName, theDesc, m12, m13, flavor);
  
  // Use external coefficients or create them locally if not supplied
  for (int i=0;i<10;i++) {
    if (c && c[i]) _c[i] = c[i];
    else {
      TString num, name;

      // Build the names of the parameters: c0, s1...s5, a1...a4
      if (i==0) name = TString("c0");
      else if (i<=5) name = TString("s") + (num+=i);
      else name = TString("a") + (num+=(i-5));

      _c[i] = new RooRealVar(TString(GetName()) + TString(".") + name,
			     TString(GetTitle()) + name +" coefficient", 
                             0.1);
      _c[i]->setConstant(kFALSE);
      _c[i]->setError(0.05);
    }
  }				
  
  setIsValid(kFALSE);
}

/// Build the PDF using the symmetric/asymmetric parametrization.
void PdfKKP2DpolyDalitz::createPdf() {
  Bdk2DpolyDalitz * thePdf = new Bdk2DpolyDalitz(TString(GetName())+".pdf",
						 TString(GetTitle())+" Pdf",
						 flavor(), 
						 BdkDalitzBase::KKP0,
						 *_m12, *_m13, *_c[0],
						 *_c[1], *_c[2], *_c[3],
						 *_c[4], *_c[5], *_c[6],
						 *_c[7], *_c[8], *_c[9],
                                                 true);
  
  setPdf(*thePdf);
  setIsValid(kTRUE);
}

