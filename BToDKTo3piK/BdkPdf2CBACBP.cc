/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdf2CBACBP.cc,v 1.1 2005/10/09 22:50:36 abi Exp $
 * Authors:
 *   R. de Sangro INFN - Frascati riccardo.desangro@slac.stanford.edu
 * Description:
 *   Class for a Crystall Bass + ARGUS pdf with it's own 
 *   variable definition to be used in conjunction with RooFitCore/Models.
 * History:
 *   17-Mar-2003 rid Created initial version from BdkPdfBifurArgus by Abi Soffer
 *
 * Copyright (C) 2003 INFN - LNF 
 *****************************************************************************/

#include "BToDKTo3piK/BdkPdf2CBACBP.hh"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooProdPdf.hh"
#include "RooFitCore/RooAddPdf.hh"

ClassImp(BdkPdf2CBACBP)

//------------------------------------------------------------
BdkPdf2CBACBP::BdkPdf2CBACBP() :
  _depCBA(0), 
  _depCBP(0), 
  _frac1(0), 
  _prod1(0), 
  _prod2(0)
{
}

//------------------------------------------------------------
BdkPdf2CBACBP::BdkPdf2CBACBP(const char * theName,
			   const char * theDesc,
			   RooRealVar & depCBA,
			   RooRealVar & depCBP,
			   int loOrder1, int hiOrder1,
			   int loOrder2, int hiOrder2) {
  
  init(theName,
       theDesc,
       depCBP,
       depCBA);
}
  
//------------------------------------------------------------
BdkPdf2CBACBP::~BdkPdf2CBACBP() {
  deleteAll();
}
  
//------------------------------------------------------------
void BdkPdf2CBACBP::init(const char * theName,
			const char * theDesc,
			RooRealVar & depCBA,
			RooRealVar & depCBP,
			int loOrder1, int hiOrder1,
			int loOrder2, int hiOrder2) {

  baseInit(theName, theDesc);
  setDependents(depCBA, depCBP);

  _loOrder1 = loOrder1;
  _hiOrder1 = hiOrder1;
  _loOrder2 = loOrder2;
  _hiOrder2 = hiOrder2;

  initParameters();
}
  
//------------------------------------------------------------
void BdkPdf2CBACBP::setDependents(RooRealVar & depCBA, RooRealVar & depCBP) {
  _depCBP = &depCBP;
  _depCBA = &depCBA;
  setIsValid(kFALSE);
}
  
//------------------------------------------------------------
void BdkPdf2CBACBP::initParameters() {
  _CBA1.init(TString(GetName()) + ".CBA1", 
	     TString(GetTitle()) + ".CBA1", 
	     *_depCBA);

  _CBA2.init(TString(GetName()) + ".CBA2", 
	     TString(GetTitle()) + ".CBA2", 
	     *_depCBA);

  _CBP1.init(TString(GetName()) + ".CBP1", 
	    TString(GetTitle()) + ".CBP1", 
	    *_depCBP, _loOrder1, _hiOrder1);

  _CBP2.init(TString(GetName()) + ".CBP2", 
	    TString(GetTitle()) + ".CBP2", 
	    *_depCBP, _loOrder2, _hiOrder2);

  _frac1 = new RooRealVar(TString(GetName()) + ".frac1", 
			  TString(GetTitle()) + ".CBP1", 
			  0.5, 0, 1);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::deleteAll() {
  delete _frac1;
  delete _prod1;
  delete _prod2;
}

//------------------------------------------------------------
void BdkPdf2CBACBP::createPdf() {
  _prod1 = new RooProdPdf(TString(GetName()) + ".prod1", 
			  TString(GetTitle()) + ".prod1", 
			  *(_CBA1.getPdf()), 
			  *(_CBP1.getPdf()));

  _prod2 = new RooProdPdf(TString(GetName()) + ".prod2", 
			  TString(GetTitle()) + ".prod2", 
			  *(_CBA2.getPdf()), 
			  *(_CBP2.getPdf()));

  RooAddPdf * thePdf = 
    new RooAddPdf(TString(GetName()) + ".pdf", 
		  TString(GetTitle()) + ".pdf", 
		  *_prod1,
		  *_prod2,
		  *_frac1);

  setPdf(*thePdf, kTRUE);
  setIsValid(kTRUE);
}


//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBAEndPoint() {
  _CBA2.linkParameters(0,0,0,0,_CBA1.argus().endPoint());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBAExp() {
  _CBA2.linkParameters(0,0,0,0,0,_CBA1.argus().exp());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBAAlpha() {
  _CBA2.linkParameters(0,0,_CBA1.cbShape().alpha());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBAEnne() {
  _CBA2.linkParameters(0,0,0,_CBA1.cbShape().enne());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBAM0() {
  _CBA2.linkParameters(_CBA1.cbShape().m0());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBASigma() {
  _CBA2.linkParameters(0,_CBA1.cbShape().sigma());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBAFracCBShape() {
  _CBA2.linkParameters(0,0,0,0,0,0,_CBA1.fracCBShape());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBPAlpha() {
  _CBP2.linkParameters(0,0,(RooRealVar*)_CBP1.cbShape().alpha());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBPEnne() {
  _CBP2.linkParameters(0,0,0,(RooRealVar*)_CBP1.cbShape().enne());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBPM0() {
  _CBP2.linkParameters((RooRealVar*)_CBP1.cbShape().m0());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBPSigma() {
  _CBP2.linkParameters(0,(RooRealVar*)_CBP1.cbShape().sigma());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::linkCBPFracCBShape() {
  _CBP2.linkParameters(0,0,0,0,0,0,(RooRealVar*)_CBP1.cbShape().m0());
  setIsValid(kFALSE);
}

//------------------------------------------------------------
void BdkPdf2CBACBP::defaultLinksCBA() {
  linkCBAEndPoint();
  linkCBAM0();
  linkCBAAlpha();
  linkCBAEnne();
  linkCBASigma();
}

//------------------------------------------------------------
void BdkPdf2CBACBP::defaultLinksCBP() {
  linkCBPAlpha();
  linkCBPEnne();
  linkCBPM0();
  linkCBPSigma();
}

//------------------------------------------------------------
void BdkPdf2CBACBP::defaultLinks() {
  defaultLinksCBA();
  defaultLinksCBP();
};
