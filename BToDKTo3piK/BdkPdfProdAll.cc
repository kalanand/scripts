/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfProdAll.cc,v 1.8 2006/05/30 04:47:25 fwinkl Exp $
 * Authors:
 *   Abi Soffer, Colorado State University, abi@slac.stanford.edu
 * Description:
 *   See BdkPdfProd.rdl
 * History:
 *   6-Mar-2004 abi Created initial version
 *
 * Copyright (C) 2004 Colorado State University and SLAC
 *****************************************************************************/

// -- CLASS DESCRIPTION [IHFPDFWRAPPER] --
// 
// Wrapper for product pdf 
// 

#include "BToDKTo3piK/BdkPdfProdAll.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooArgList.hh"
#include "RooFitCore/RooProdPdf.hh"

ClassImp(BdkPdfProdAll)


//-----------------------------------------------------------------
BdkPdfProdAll::BdkPdfProdAll() {
  // constructor
  // common initialization:
  commonInit();
}
 
//-----------------------------------------------------------------
BdkPdfProdAll::BdkPdfProdAll(const char * theName,
			     const char * theDesc,
			     RooRealVar& deltaE,
			     RooRealVar& mD,
			     RooRealVar& mes,
			     RooRealVar& nnCont,
			     RooRealVar& nnComb,
			     BdkPdfAbsBase* nnContPdf,
			     BdkPdfAbsBase* nnCombPdf,
			     BdkPdfAbsBase* dalitzPdf) {
  // Initializing constructor
  init(theName, 
       theDesc, 
       deltaE,
       mD,
       mes,
       nnCont,
       nnComb,
       nnContPdf,
       nnCombPdf,
       dalitzPdf);
}

//-----------------------------------------------------------------
BdkPdfProdAll::BdkPdfProdAll(const char * theName,
			     const char * theDesc,
			     BdkPdfAbsBase* deltaEPdf,
			     BdkPdfAbsBase* mDPdf,
			     BdkPdfAbsBase* mesPdf,
			     BdkPdfAbsBase* nnContPdf,
			     BdkPdfAbsBase* nnCombPdf,
			     BdkPdfAbsBase* dalitzPdf) {
  // Initializing constructor
  init(theName,
       theDesc,
       deltaEPdf,
       mDPdf,
       mesPdf,
       nnContPdf,
       nnCombPdf,
       dalitzPdf);
}
  
//-----------------------------------------------------------------
BdkPdfProdAll::~BdkPdfProdAll() {
  // Destructor. The pointer data members never need to be deleted.
}

//-----------------------------------------------------------------
void BdkPdfProdAll::init(const char * theName,
			 const char * theDesc,
			 RooRealVar& deltaE,
			 RooRealVar& mD,	
			 RooRealVar& mes,	
			 RooRealVar& nnCont,	
			 RooRealVar& nnComb,	
                         BdkPdfAbsBase* nnContPdf,
                         BdkPdfAbsBase* nnCombPdf,
			 BdkPdfAbsBase* dalitzPdf) {

  // Do base class initializer
  baseInit(theName, theDesc);

  // common initialization:
  commonInit();

  // set the dependents of the owned objects:
  _ownedDeltaEPdf.setDependent(deltaE);
  _ownedMDPdf.setDependent(mD);
  _ownedMesPdf.setDependent(mes);
  _ownedNnContPdf.setDependent(nnCont);
  _ownedNnCombPdf.setDependent(nnComb);

  // if external PDFs are supplied use them
  if (nnContPdf) _nnContPdf = nnContPdf;
  if (nnCombPdf) _nnCombPdf = nnCombPdf;

  // Set the external pointer PDFs;
  _dalitzPdf = dalitzPdf;

  initParameters();
}

//-----------------------------------------------------------------
void BdkPdfProdAll::init(const char * theName,
			 const char * theDesc,
			 BdkPdfAbsBase* deltaEPdf,
			 BdkPdfAbsBase* mDPdf,
			 BdkPdfAbsBase* mesPdf,
			 BdkPdfAbsBase* nnContPdf,
			 BdkPdfAbsBase* nnCombPdf,
			 BdkPdfAbsBase* dalitzPdf) {
  // Do base class initializer
  baseInit(theName, theDesc);

  // common initialization:
  commonInit();

  _deltaEPdf = deltaEPdf;
  _mDPdf = mDPdf;
  _mesPdf = mesPdf;
  _nnContPdf = nnContPdf;
  _nnCombPdf = nnCombPdf;
  _dalitzPdf = dalitzPdf;
}

//-----------------------------------------------------------------
void BdkPdfProdAll::initParameters() {
  // Initialize the owned PDFs and set the units of their parameters:  
  // deltae:
  _ownedDeltaEPdf.init(TString(GetName()) + ".deltaEPdf",
                       TString(GetTitle()) + ".deltaEPdf OWNED",
                       *_ownedDeltaEPdf.dependent(),0, 2);

  _ownedDeltaEPdf.gauss().b()->setVal(0.0);
  _ownedDeltaEPdf.gauss().s()->setVal(0.019);
  
  // mD:
  _ownedMDPdf.init(TString(GetName()) + ".mDPdf",
		   TString(GetTitle()) + ".mDPdf OWNED",
		   *_ownedMDPdf.dependent(),0, 2);   

  // mes:
  _ownedMesPdf.init(TString(GetName()) + ".mesPdf",
		   TString(GetTitle()) + ".mesPdf OWNED",
		   *_ownedMesPdf.dependent());   

  // nnCont:
  _ownedNnContPdf.init(TString(GetName()) + ".nnContPdf",
                       TString(GetTitle()) + ".nnContPdf OWNED",
                       *_ownedNnContPdf.dependent());  

  // nnComb:
  _ownedNnCombPdf.init(TString(GetName()) + ".nnCombPdf",
                       TString(GetTitle()) + ".nnCombPdf OWNED",
                       *_ownedNnCombPdf.dependent());  


}

//-----------------------------------------------------------------
void BdkPdfProdAll::commonInit() {
  // initialize external pointers:
  _dalitzPdf = 0;
  _nnContPdf = 0;
  _nnCombPdf = 0;
  _mDPdf = 0;
  _mesPdf = 0;

  // initialize these PDFs to the owned ones:
  useOwnedDeltaEPdf();
  useOwnedMDPdf();
  useOwnedMesPdf();
  useOwnedNnContPdf();
  useOwnedNnCombPdf();

  // by default, include all variablesin the product PDF:
  _includeDeltaE = kTRUE;
  _includeMD = kTRUE;
  _includeMes = kTRUE;
  _includeDalitz = kTRUE;
  _includeNnCont = kTRUE;
  _includeNnComb = kTRUE;

}

//-----------------------------------------------------------------
void BdkPdfProdAll::setDeltaEPdf(BdkPdfAbsBase& theDeltaEPdf) {
  _deltaEPdf = &theDeltaEPdf;
  setIsValid(kFALSE);
}

//-----------------------------------------------------------------
void BdkPdfProdAll::setMDPdf(BdkPdfAbsBase& theMDPdf) {
  _mDPdf = &theMDPdf;
  setIsValid(kFALSE);
}

//-----------------------------------------------------------------
void BdkPdfProdAll::setMesPdf(BdkPdfAbsBase& theMesPdf) {
  _mesPdf = &theMesPdf;
  setIsValid(kFALSE);
}

//-----------------------------------------------------------------
void BdkPdfProdAll::setNnContPdf(BdkPdfAbsBase& theNnContPdf) {
  _nnContPdf = &theNnContPdf;
  setIsValid(kFALSE);
}

//-----------------------------------------------------------------
void BdkPdfProdAll::setNnCombPdf(BdkPdfAbsBase& theNnCombPdf) {
  _nnCombPdf = &theNnCombPdf;
  setIsValid(kFALSE);
}

//-----------------------------------------------------------------
void BdkPdfProdAll::setDalitzPdf(BdkPdfAbsBase& theDalitzPdf) {
  _dalitzPdf = &theDalitzPdf;
  setIsValid(kFALSE);
}

//-----------------------------------------------------------------
void BdkPdfProdAll::createPdf() {
  // Add the components of this PDF:
  removeComponents();
  if (_includeDeltaE && 0 != getDeltaEPdf()) addComponent(getDeltaEPdf());
  if (_includeDalitz && 0 != getDalitzPdf()) addComponent(getDalitzPdf());
  if (_includeNnCont && 0 != getNnContPdf()) addComponent(getNnContPdf());
  if (_includeNnComb && 0 != getNnCombPdf()) addComponent(getNnCombPdf());
  if (_includeMD &&     0 != getMDPdf())     addComponent(getMDPdf());
  if (_includeMes &&    0 != getMesPdf())    addComponent(getMesPdf());

  // Makes both PDFs of this wrapper
  RooAbsPdf * thePdf = new RooProdPdf(GetName() + TString(".pdf"),
				      GetTitle() + TString("PDF"),
				      componentPdfs());

  setPdf(*thePdf, kTRUE);   // owned

  setIsValid(kTRUE);
}

//-----------------------------------------------------------------
BdkPdfAbsBase * BdkPdfProdAll::getVarPdf(BdkPdfProdAll::Var var) {
  switch(var) {
  case DELTAE: return getDeltaEPdf();
  case MD: return getMDPdf();
  case MES: return getMesPdf();
  case DALITZ: return getDalitzPdf();
  case NNCONT: return getNnContPdf();
  case NNCOMB: return getNnCombPdf();
  }
  return 0;
}

//-----------------------------------------------------------------
Bool_t BdkPdfProdAll::deltaEPdfIsOwned() const {
  Bool_t result = kFALSE;
  if((BdkPdfAbsBase *)(&_ownedDeltaEPdf) == _deltaEPdf) {
        result = kTRUE;
  }
  return result;
}

//-----------------------------------------------------------------
Bool_t BdkPdfProdAll::mDPdfIsOwned() const {
  Bool_t result = kFALSE;
  if((BdkPdfAbsBase *)(&_ownedMDPdf) == _mDPdf) {
        result = kTRUE;
  }
  return result;
}

//-----------------------------------------------------------------
Bool_t BdkPdfProdAll::mesPdfIsOwned() const {
  Bool_t result = kFALSE;
  if((BdkPdfAbsBase *)(&_ownedMesPdf) == _mesPdf) {
        result = kTRUE;
  }
  return result;
}

//-----------------------------------------------------------------
Bool_t BdkPdfProdAll::nnContPdfIsOwned() const {
  Bool_t result = kFALSE;
  if((BdkPdfAbsBase *)(&_ownedNnContPdf) == _nnContPdf) {
        result = kTRUE;
  }
  return result;
}

//-----------------------------------------------------------------
Bool_t BdkPdfProdAll::nnCombPdfIsOwned() const {
  Bool_t result = kFALSE;
  if((BdkPdfAbsBase *)(&_ownedNnCombPdf) == _nnCombPdf) {
        result = kTRUE;
  }
  return result;
}

//-----------------------------------------------------------------
void BdkPdfProdAll::includeDeltaE(Bool_t include) {
  _includeDeltaE = include;
  setIsValid(kFALSE);
}

//-----------------------------------------------------------------
void BdkPdfProdAll::includeMD(Bool_t include) {
  _includeMD = include;
  setIsValid(kFALSE);
}

//-----------------------------------------------------------------
void BdkPdfProdAll::includeMes(Bool_t include) {
  _includeMes = include;
  setIsValid(kFALSE);
}

//-----------------------------------------------------------------
void BdkPdfProdAll::includeDalitz(Bool_t include) {
  _includeDalitz = include;
  setIsValid(kFALSE);
}

//-----------------------------------------------------------------
void BdkPdfProdAll::includeNnCont(Bool_t include) {
  _includeNnCont = include;
  setIsValid(kFALSE);
}

//-----------------------------------------------------------------
void BdkPdfProdAll::includeNnComb(Bool_t include) {
  _includeNnComb = include;
  setIsValid(kFALSE);
}


//-----------------------------------------------------------------
void BdkPdfProdAll::printVariables() {
  cout << "(";
  if (_includeDeltaE) cout << "DeltaE";
  if (_includeMD) cout << ", MD";
  if (_includeMes) cout << ", Mes";
  if (_includeDalitz) cout << ", Dalitz";
  if (_includeNnCont) cout << ", NnCont";
  if (_includeNnComb) cout << ", NnComb";

  cout << ")" << endl;
}
