/**************************************************************************
*                    BdkPdfHist.cc                                        *
***************************************************************************/

#include "RooFitCore/RooDataHist.hh"
#include "RooFitCore/RooHistPdf.hh"
#include "RooFitCore/RooArgSet.hh"

#include "BToDKTo3piK/BdkPdfHist.hh"

ClassImp(BdkPdfHist);

BdkPdfHist::BdkPdfHist() {
}

BdkPdfHist::BdkPdfHist(	const char * theName,
               		const char * theDesc,
               		const RooArgSet & vars,
               		const RooDataHist & dhist,
               		Int_t intOrder) {
//
// Constructor			
// 
   init(theName, theDesc, vars, dhist, intOrder);
}

BdkPdfHist::~BdkPdfHist() {} //Destructor

void BdkPdfHist::init(const char * theName,
		      const char * theDesc,
		      const RooArgSet & vars,
		      const RooDataHist & dhist,
		      Int_t intOrder) {
//
// Initialization
//
   baseInit(theName, theDesc);
   setDependents(vars);
   setOrder(intOrder);
   setDataHist(dhist);
   initParameters();
}

void BdkPdfHist::initParameters() {
//
//Initialize parameters
//no parameters are needed 
//   
}

void BdkPdfHist::setDependents(const RooArgSet & vars) {
	_vars.add(vars);
}


void BdkPdfHist::setOrder(Int_t intOrder) {
        _order = intOrder;
}

void BdkPdfHist::setDataHist(const RooDataHist& dhist) {
	_dataHist =(RooDataHist *) &dhist;
} 

void BdkPdfHist::createPdf() {
//
// Build the HistPdf
//
   _thePdf = new RooHistPdf(TString(GetName())+".pdf",
			    TString(GetTitle())+" Pdf",
			    _vars, *_dataHist, _order);
   setIsValid(kTRUE);
}

RooArgSet BdkPdfHist::dependents() {
//
// get dependents( list )
//
   return _vars;
}	  

