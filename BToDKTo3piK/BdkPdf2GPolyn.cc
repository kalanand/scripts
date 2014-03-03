/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 *
 * Package: 
 * Gauss Poly + Gauss 
 *
 * Description:
 *   To be used in conjunction with RooFitCore/Models.
 *
 *****************************************************************************/

// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
// 

#include <iostream>
using std::cout;
using std::endl;

#include "TString.h"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooArgSet.hh"

#include "BToDKTo3piK/BdkPdf2GPolyn.hh"

ClassImp(BdkPdf2GPolyn)

//--------------------------------------------------------------
BdkPdf2GPolyn::BdkPdf2GPolyn() {
  // Default Constructor:
  addComponents();
}

//--------------------------------------------------------------
BdkPdf2GPolyn::BdkPdf2GPolyn(const char * theName,
				   const char * theDesc,
				   RooRealVar & dependent) {
  // Initializing Constructor:
  addComponents();
  init(theName, theDesc, dependent);
}

//--------------------------------------------------------------
BdkPdf2GPolyn::~BdkPdf2GPolyn() {
  // Destructor:
}

//--------------------------------------------------------------
BdkPdf2GPolyn::BdkPdf2GPolyn(const BdkPdf2GPolyn & pdf) {
// Dummy copy constructor:
    cout << "BdkPdf2GPolyn::BdkPdf2GPolyn -- dummy copy constructor"
	 << endl;
}


//--------------------------------------------------------------
void BdkPdf2GPolyn::init(const char * theName,
			    const char * theDesc,
			    RooRealVar & dependent) {
  // Initializer:

  baseInit(theName, theDesc);
  setDependent(dependent);
  initParameters();
}

//--------------------------------------------------------------
void BdkPdf2GPolyn::setDependent(RooRealVar & dependent) {
  // Sets the dependent variable:
  _dependent = &dependent;
  _gaussPolyn.setDependent(dependent);
  _gauss.setDependent(dependent);

  setIsValid(kFALSE);
}

//--------------------------------------------------------------
void BdkPdf2GPolyn::initParameters() {
  // Initializes the contained RRV's and PDF wrappers:

  _fracGP = new RooRealVar(TString(GetName()) + ".gpFrac",
                          TString(GetTitle()) + " gpFrac",
                          0.5, 0, 1);

  _gaussPolyn.init(GetName() + TString(".gaussPolynPdf"), 
	      GetTitle() + TString(" Gausss Polyn Pdf "), 
	      *_dependent, 0 , 2);
  _gaussPolyn.gauss().b()->setVal(0);
  _gaussPolyn.gauss().s()->setVal(0.011);

  _gauss.init(GetName() + TString(".gauss"), 
	      GetTitle() + TString(" Gauss"), 
	      *_dependent);
  _gauss.b()->setVal(-0.1);
  _gauss.s()->setVal(0.018);
}

//--------------------------------------------------------------
void BdkPdf2GPolyn::createPdf() {

  // Builds the PDF. The PDFs of the contained wrappers are automatically built
  RooAddPdf * thePdf = new RooAddPdf(GetName() + TString(".thePdf"), 
			  GetTitle() + TString(" thePdf"), 
                          *_gaussPolyn.getPdf(), 
			  *_gauss.getPdf(),
			  *_fracGP);

  setPdf(*thePdf, kTRUE);
  setIsValid(kTRUE);
}

//--------------------------------------------------------------
void BdkPdf2GPolyn::addComponents() {
  // Adds the wrappers to the composite's list:
  addComponent(_gaussPolyn);
  addComponent(_gauss);
}

//--------------------------------------------------------------
void BdkPdf2GPolyn::linkParameters(BdkPdf2GPolyn & pdf2)
{
  // Copies RRV pointers:
  _gaussPolyn.linkParameters(pdf2.gaussPolyn());
  _gauss.linkResolution(pdf2.gauss());
  _fracGP = pdf2.fracGP();
  setIsValid(kFALSE);
}

//--------------------------------------------------------------
void BdkPdf2GPolyn::linkMeanVal(Bool_t lk ) {
  if( kTRUE == lk ) { 
    _gauss.linkResolution(_gaussPolyn.gauss().b());
    setIsValid(kFALSE);
  }
}
