/***************
* class BdkPdfOnRes
***************/

#include <math.h>
#include <iostream>
using std::cout;
using std::endl;

#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooCategory.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooFormulaVar.hh"
#include "BToDKTo3piK/BdkPdfHolder.hh"
#include "BToDKTo3piK/BdkPdfOnRes.hh"
#include "RooFitCore/RooExtendPdf.hh"
#include "RooFitCore/RooSimultaneous.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooListProxy.hh"
#include "RooFitCore/RooFitResult.hh"
#include "RooFitCore/RooMinuit.hh"
#include "RooFitCore/RooNLLVar.hh"
#include "RooFitCore/RooGlobalFunc.hh"
#include "RooFitCore/RooRandom.hh"
#include "BToDKTo3piK/BdkOnResNLLYields.hh"
#include "BToDKTo3piK/BdkPdfDKDalitz.hh"

using namespace RooFit;

ClassImp(BdkPdfOnRes)


//------------------------------------------
BdkPdfOnRes::BdkPdfOnRes() : // constructor
  _charge(0),
  _yieldFitResult(0),
  _xyFitResult(0),
  _useMinos(kFALSE),
  _minuit(0),
  _minuitPrintLevel(1),
  _nllYieldsSystBit(-1)
{ 
} 

//------------------------------------------
BdkPdfOnRes::BdkPdfOnRes(const char * theName,
			 const char * theDesc,
                         int maxRange,
			 RooCategory& charge,
			 RooRealVar& deltaE,
			 RooRealVar& mD0,
			 RooRealVar& mes,
                         RooRealVar& nnCont,                           
                         RooRealVar& nnComb,  
                         BdkPdfHolder* nnContHolder,
                         BdkPdfHolder* nnCombHolder,
			 BdkPdfHolder* dalitzHolderN,
			 BdkPdfHolder* dalitzHolderP,
                         BdkPdfOnRes::Options opt) : // constructor
  _yieldFitResult(0),
  _xyFitResult(0),
  _useMinos(kFALSE),
  _minuit(0),
  _minuitPrintLevel(1),
  _nllYieldsSystBit(-1)
{ 
  init(theName, theDesc,
       maxRange, 
       charge, deltaE,
       mD0,
       mes,
       nnCont,
       nnComb,
       nnContHolder,
       nnCombHolder,
       dalitzHolderN,
       dalitzHolderP,
       opt);
}

//------------------------------------------
BdkPdfOnRes::~BdkPdfOnRes() {
} //deconstruct;

//-----------------------------------------------
void BdkPdfOnRes::init( const char * theName,
            		const char * theDesc,
	                int maxRange,
			RooCategory& charge,
            		RooRealVar& deltaE,
            		RooRealVar& mD0,
            		RooRealVar& mes,
                        RooRealVar& nnCont,  
                        RooRealVar& nnComb,
                        BdkPdfHolder* nnContHolder,
                        BdkPdfHolder* nnCombHolder,
            		BdkPdfHolder * dalitzHolderN,
            		BdkPdfHolder * dalitzHolderP,
                        BdkPdfOnRes::Options opt) {

  baseInit(theName, theDesc);  
  setCharge(charge);
     
  // Initialize the B+ PDFs:

  _sigGoodD0ProdN.init(TString(GetName())+".sigGoodD0N",
                       TString(GetTitle())+" sigGoodD0 for B-",
                       deltaE, mD0, mes, nnCont, nnComb,
                       nnContHolder ? nnContHolder->sigGoodD0() : 0,
                       nnCombHolder ? nnCombHolder->sigGoodD0() : 0,
                       dalitzHolderN ? dalitzHolderN->sigGoodD0(): 0);
  
  _sigBadD0ProdN.init(TString(GetName())+".sigBadD0N",
                      TString(GetTitle())+" sigBadD0 for B-",
                      deltaE, mD0, mes, nnCont, nnComb,
                      nnContHolder ? nnContHolder->sigBadD0() : 0,
                      nnCombHolder ? nnCombHolder->sigBadD0() : 0,
                      dalitzHolderN ? dalitzHolderN->sigBadD0(): 0);
  
  _DpiGoodD0ProdN.init(TString(GetName())+".DpiGoodD0N",
                       TString(GetTitle())+" DpiGoodD0 for B-",
                       deltaE, mD0, mes, nnCont, nnComb,
                       nnContHolder ? nnContHolder->DpiGoodD0() : 0,
                       nnCombHolder ? nnCombHolder->DpiGoodD0() : 0,
                       dalitzHolderN ? dalitzHolderN->DpiGoodD0(): 0);
  
  _DpiBadD0ProdN.init(TString(GetName())+".DpiBadD0N",
                      TString(GetTitle())+" DpiBadD0 for B-",
                      deltaE, mD0, mes, nnCont, nnComb,
                      nnContHolder ? nnContHolder->DpiBadD0() : 0,
                      nnCombHolder ? nnCombHolder->DpiBadD0() : 0,
                      dalitzHolderN ? dalitzHolderN->DpiBadD0(): 0);
  
  _DPiXProdN.init(TString(GetName())+".DPiXN",
		  TString(GetTitle())+" DPiX for B-",
		  deltaE, mD0, mes, nnCont, nnComb,
                  nnContHolder ? nnContHolder->DPiX() : 0,
                  nnCombHolder ? nnCombHolder->DPiX() : 0,
		  dalitzHolderN ? dalitzHolderN->DPiX(): 0);
    
  _DKXProdN.init(TString(GetName())+".DKXN",
                 TString(GetTitle())+" DKX for B-",
                 deltaE, mD0, mes, nnCont, nnComb,
                 nnContHolder ? nnContHolder->DKX() : 0,
                 nnCombHolder ? nnCombHolder->DKX() : 0,
                 dalitzHolderN ? dalitzHolderN->DKX(): 0);
  
  _BBGoodD0ProdN.init(TString(GetName())+".BBgoodD0N",
                      TString(GetTitle())+" BBgoodD0 for B-",
                      deltaE, mD0, mes, nnCont, nnComb,
                      nnContHolder ? nnContHolder->BBGoodD0() : 0,
                      nnCombHolder ? nnCombHolder->BBGoodD0() : 0,
                      dalitzHolderN ? dalitzHolderN->BBGoodD0(): 0);
  
  _BBBadD0ProdN.init(TString(GetName())+".BBbadD0N",  
                     TString(GetTitle())+" BBbadD0 for B-",
                     deltaE, mD0, mes, nnCont, nnComb,
                     nnContHolder ? nnContHolder->BBBadD0() : 0,
                     nnCombHolder ? nnCombHolder->BBBadD0() : 0,
                     dalitzHolderN ? dalitzHolderN->BBBadD0(): 0);
  
  _qqGoodD0ProdN.init(TString(GetName())+".ContGoodD0N",
                      TString(GetTitle())+" ContGooddD0 for B-",
                      deltaE, mD0, mes, nnCont, nnComb,
                      nnContHolder ? nnContHolder->qqGoodD0() : 0,
                      nnCombHolder ? nnCombHolder->qqGoodD0() : 0,
                      dalitzHolderN ? dalitzHolderN->qqGoodD0(): 0);

  _qqBadD0ProdN.init(TString(GetName())+".ContBadD0N",
                     TString(GetTitle())+" ContBadD0 for B-",
                     deltaE,mD0, mes, nnCont, nnComb,
                     nnContHolder ? nnContHolder->qqBadD0() : 0,
                     nnCombHolder ? nnCombHolder->qqBadD0() : 0,
                     dalitzHolderN ? dalitzHolderN->qqBadD0(): 0);

  // Initialize the B+ PDFs, using the same 1-D PDFs as for the B-,
  // except for the Dalitz part:

  _sigGoodD0ProdP.init(TString(GetName())+".sigGoodD0P",
                       TString(GetTitle())+" sigGoodD0 for B+",
		       _sigGoodD0ProdN.getDeltaEPdf(),
		       _sigGoodD0ProdN.getMDPdf(),
		       _sigGoodD0ProdN.getMesPdf(),
		       _sigGoodD0ProdN.getNnContPdf(),
		       _sigGoodD0ProdN.getNnCombPdf(),
                       dalitzHolderP ? dalitzHolderP->sigGoodD0(): 0);
  
  _sigBadD0ProdP.init(TString(GetName())+".sigBadD0P",
                      TString(GetTitle())+" sigBadD0 for B+",
		      _sigBadD0ProdN.getDeltaEPdf(),
		      _sigBadD0ProdN.getMDPdf(),
		      _sigBadD0ProdN.getMesPdf(),
		      _sigBadD0ProdN.getNnContPdf(),
		      _sigBadD0ProdN.getNnCombPdf(),
                      dalitzHolderP ? dalitzHolderP->sigBadD0(): 0);
  
  _DpiGoodD0ProdP.init(TString(GetName())+".DpiGoodD0P",
                       TString(GetTitle())+" DpiGoodD0 for B+",
		       _DpiGoodD0ProdN.getDeltaEPdf(),
		       _DpiGoodD0ProdN.getMDPdf(),
		       _DpiGoodD0ProdN.getMesPdf(),
		       _DpiGoodD0ProdN.getNnContPdf(),
		       _DpiGoodD0ProdN.getNnCombPdf(),
                       dalitzHolderP ? dalitzHolderP->DpiGoodD0(): 0);
  
  _DpiBadD0ProdP.init(TString(GetName())+".DpiBadD0P",
                      TString(GetTitle())+" DpiBadD0 for B+",
		      _DpiBadD0ProdN.getDeltaEPdf(),
		      _DpiBadD0ProdN.getMDPdf(),
		      _DpiBadD0ProdN.getMesPdf(),
		      _DpiBadD0ProdN.getNnContPdf(),
		      _DpiBadD0ProdN.getNnCombPdf(),
                      dalitzHolderP ? dalitzHolderP->DpiBadD0(): 0);
  
  _DPiXProdP.init(TString(GetName())+".DPiXP",
		  TString(GetTitle())+" DPiX for B+",
		  _DPiXProdN.getDeltaEPdf(),
		  _DPiXProdN.getMDPdf(),
		  _DPiXProdN.getMesPdf(),
		  _DPiXProdN.getNnContPdf(),
		  _DPiXProdN.getNnCombPdf(),
		  dalitzHolderP ? dalitzHolderP->DPiX(): 0);
    
  _DKXProdP.init(TString(GetName())+".DKXP",
                 TString(GetTitle())+" DKX for B+",
		 _DKXProdN.getDeltaEPdf(),
		 _DKXProdN.getMDPdf(),
		 _DKXProdN.getMesPdf(),
		 _DKXProdN.getNnContPdf(),
		 _DKXProdN.getNnCombPdf(),
                 dalitzHolderP ? dalitzHolderP->DKX(): 0);
  
  _BBGoodD0ProdP.init(TString(GetName())+".BBgoodD0P",
                      TString(GetTitle())+" BBgoodD0 for B+",
		      _BBGoodD0ProdN.getDeltaEPdf(),
		      _BBGoodD0ProdN.getMDPdf(),
		      _BBGoodD0ProdN.getMesPdf(),
		      _BBGoodD0ProdN.getNnContPdf(),
		      _BBGoodD0ProdN.getNnCombPdf(),
                      dalitzHolderP ? dalitzHolderP->BBGoodD0(): 0);
  
  _BBBadD0ProdP.init(TString(GetName())+".BBbadD0P",  
                     TString(GetTitle())+" BBbadD0 for B+",
		     _BBBadD0ProdN.getDeltaEPdf(),
		     _BBBadD0ProdN.getMDPdf(),
		     _BBBadD0ProdN.getMesPdf(),
		     _BBBadD0ProdN.getNnContPdf(),
		     _BBBadD0ProdN.getNnCombPdf(),
                     dalitzHolderP ? dalitzHolderP->BBBadD0(): 0);
  
  _qqGoodD0ProdP.init(TString(GetName())+".ContGoodD0P",
                      TString(GetTitle())+" ContGooddD0 for B+",
		      _qqGoodD0ProdN.getDeltaEPdf(),
		      _qqGoodD0ProdN.getMDPdf(),
		      _qqGoodD0ProdN.getMesPdf(),
		      _qqGoodD0ProdN.getNnContPdf(),
		      _qqGoodD0ProdN.getNnCombPdf(),
                      dalitzHolderP ? dalitzHolderP->qqGoodD0(): 0);
  
  _qqBadD0ProdP.init(TString(GetName())+".ContBadD0P",
                     TString(GetTitle())+" ContBadD0 for B+",
		     _qqBadD0ProdN.getDeltaEPdf(),
		     _qqBadD0ProdN.getMDPdf(),
		     _qqBadD0ProdN.getMesPdf(),
		     _qqBadD0ProdN.getNnContPdf(),
		     _qqBadD0ProdN.getNnCombPdf(),
                     dalitzHolderP ? dalitzHolderP->qqBadD0(): 0);


  setMaxRange(maxRange);	

  initParameters(opt);
} 

//------------------------------------------
void BdkPdfOnRes::initParameters(BdkPdfOnRes::Options opt) {
  _sigGoodD0NumEvts = new RooRealVar(TString(GetName())+".sigGoodD0NumEvts", 
				    TString(GetTitle())+" sigGoodD0 # events", 
                                     _maxRange, 0, _maxRange);
  _sigGoodD0NumEvts->setError(10);

  _sigBadD0Frac = new RooRealVar(TString(GetName())+".sigBadD0Frac", 
				   TString(GetTitle())+" sigBadD0/sigGoodD0", 
				   _maxRange, 0, 100);
  _sigBadD0Frac->setError(0.01);

  _DpiGoodD0NumEvts = new RooRealVar(TString(GetName())+".DpiGoodD0NumEvts", 
				    TString(GetTitle())+" DpiGoodD0 # events", 
				    _maxRange, 0, _maxRange);
  _DpiGoodD0NumEvts->setError(10);

  _DpiBadD0Frac = new RooRealVar(TString(GetName())+".DpiBadD0Frac", 
				   TString(GetTitle())+" DpiBadD0/DpiGoodD0", 
				   _maxRange, 0, 100);
  _DpiBadD0Frac->setError(0.01);  

//  _DPiXNumEvts = new RooRealVar(TString(GetName())+".DPiXNumEvts", 
//				TString(GetTitle())+" DPiX # events", 
//				_maxRange, 0, _maxRange);
//  _DPiXNumEvts->setError(10);

//  _DKXNumEvts = new RooRealVar(TString(GetName())+".DKXNumEvts", 
//				    TString(GetTitle())+" DKX # events", 
//				    _maxRange, 0, _maxRange);
//  _DKXNumEvts->setError(10); 

  _totBBNumEvts = new RooRealVar(TString(GetName())+".totBBNumEvts",
                               TString(GetTitle())+" All other BB # events",
                               _maxRange, 0, _maxRange);
  _totBBNumEvts->setError(10);

//  _DPiXFrac = new RooRealVar(TString(GetName())+".DPiXFrac",
//                            TString(GetTitle())+" DKOther/totBBDKPi",
//                            _maxRange, 0, 100);
//  _DPiXFrac->setError(0.01);

  _DKXFrac = new RooRealVar(TString(GetName())+".DKXFrac",
		            TString(GetTitle())+" DKOther/numDKP",
			    _maxRange, 0, _maxRange);
  _DKXFrac->setError(0.01);

//  _BBBadD0NumEvts = new RooRealVar(TString(GetName())+".BBBadD0NumEvts", 
//				   TString(GetTitle())+" BBBadD0 # events", 
//				   _maxRange, 0, _maxRange);
//  _BBBadD0NumEvts->setError(10); 

  _DPiXFrac = new RooRealVar(TString(GetName())+".DPiXFrac",
                            TString(GetTitle())+" D0 pi Other/totBBDKPi",
                            _maxRange, 0, _maxRange);
  _DPiXFrac->setError(0.01);
  
  _BBGoodD0Frac = new RooRealVar(TString(GetName())+".BBGoodD0Frac", 
				  TString(GetTitle())+" BBGoodD0/BBBadD0", 
				  _maxRange, 0, _maxRange);
  _BBGoodD0Frac->setError(0.01); 

  _qqBadD0NumEvts = new RooRealVar(TString(GetName())+".qqBadD0NumEvts", 
				   TString(GetTitle())+" qqBadD0 # events", 
				   _maxRange, 0, _maxRange);
  _qqBadD0NumEvts->setError(10);

  _qqGoodD0Frac = new RooRealVar(TString(GetName())+".qqGoodD0Frac", 
				  TString(GetTitle())+" qqGoodD0/qqBadD0", 
				  _maxRange, 0, _maxRange);
  _qqGoodD0Frac->setError(0.01);

  // Setup the formulas for the #'s of bad D's:
    _sigBadD0NumEvts = 
    new RooFormulaVar(TString(GetName())+".sigBadD0NumEvts",
		      TString(GetTitle())+" sigBadD0NumEvts",
		      TString(_sigGoodD0NumEvts->GetName()) + "*" 
		         + TString(_sigBadD0Frac->GetName()), 
		      RooArgList(*_sigGoodD0NumEvts, *_sigBadD0Frac));

  _DpiBadD0NumEvts = 
    new RooFormulaVar(TString(GetName())+".DpiBadD0NumEvts",
		      TString(GetTitle())+" DpiBadD0NumEvts",
		      TString(_DpiGoodD0NumEvts->GetName()) + "*" 
		         + TString(_DpiBadD0Frac->GetName()), 
		      RooArgList(*_DpiGoodD0NumEvts, *_DpiBadD0Frac));

  _BBGoodD0NumEvts = 
    new RooFormulaVar(TString(GetName())+".BBGoodD0NumEvts",
		      TString(GetTitle())+" BBGoodD0NumEvts",
		      TString(_totBBNumEvts->GetName()) + "*" 
		         + TString(_BBGoodD0Frac->GetName()), 
		      RooArgList(*_totBBNumEvts, *_BBGoodD0Frac));

  _qqGoodD0NumEvts = 
    new RooFormulaVar(TString(GetName())+".qqGoodD0NumEvts",
		      TString(GetTitle())+" qqGoodD0NumEvts",
		      TString(_qqBadD0NumEvts->GetName()) + "*" 
		         + TString(_qqGoodD0Frac->GetName()), 
		      RooArgList(*_qqBadD0NumEvts, *_qqGoodD0Frac));

  _DKXNumEvts =
    new RooFormulaVar(TString(GetName())+".DKXNumEvts",
                      TString(GetTitle())+" DKXNumEvts",
                      TString(_totBBNumEvts->GetName()) + "*"
 			 + TString(_DKXFrac->GetName()) + "*" 
                         + TString(_DPiXFrac->GetName()),
		      RooArgList(*_totBBNumEvts, *_DKXFrac, *_DPiXFrac));  

  _DPiXNumEvts =
    new RooFormulaVar(TString(GetName())+".DPiXNumEvts",
                      TString(GetTitle())+" DPiXNumEvts",
                      TString(_totBBNumEvts->GetName())+"*"
                        + TString(_DPiXFrac->GetName()),
                      RooArgList(*_totBBNumEvts, *_DPiXFrac));

 
  _BBBadD0NumEvts =
    new RooFormulaVar(TString(GetName())+".BBBadD0NumEvts",
                      TString(GetTitle())+" BBBadD0NumEvts",
                      TString(_totBBNumEvts->GetName())+"*(1.0 - ("
                         + TString(_DKXFrac->GetName()) +" + 1.0 ) *"  
                         + TString(_DPiXFrac->GetName()) + " )" ,
                      RooArgList(*_totBBNumEvts, *_DKXFrac, *_DPiXFrac));

  //set up connection between _numEvt between all these above numbers
  _numEvt[0] = _sigBadD0NumEvts; 
  _numEvt[1] = _sigGoodD0NumEvts;
  _numEvt[2] = _DpiBadD0NumEvts;
  _numEvt[3] = _DpiGoodD0NumEvts;
  _numEvt[4] = _DPiXNumEvts;
  _numEvt[5] = _DKXNumEvts;
  _numEvt[6] = _BBBadD0NumEvts;
  _numEvt[7] = _BBGoodD0NumEvts;
  _numEvt[8] = _qqBadD0NumEvts;
  _numEvt[9] = _qqGoodD0NumEvts; 

  // initialize _globalAsym:
  _globalAsym = new 
    RooRealVar(TString(GetName()) + ".globalAsym", 
	       TString(GetTitle()) + " global asymmetry", 
	       0, -1, 1);


  // for each evt type
  for (int t = 0; t < BdkEvtTypes::NTYPES; ++t){
    // Initialize the asymmetries:
    _typeAsym[t] = new 
      RooRealVar(TString(GetName()) + TString(".asym") + BdkEvtTypes::name(t),
		 TString(GetTitle()) + TString(" asymmetry of ") 
		 + BdkEvtTypes::name(t),
		 0, -1, 1);
  }

  // link bad D asym to good D asym if requested (for systematics studies)
  if (opt && BdkPdfOnRes::linkGoodBadSigAsym) 
    asymLink(BdkEvtTypes::SIG_BAD_D,BdkEvtTypes::SIG_GOOD_D);


  for (int t = 0; t < BdkEvtTypes::NTYPES; ++t){
    // The actual asymmetries are typeAsym + globalAsym:

    _asym[t] = new 
      RooFormulaVar(TString(GetName()) + TString(".finalAsym") 
		      + BdkEvtTypes::name(t),
		    TString(GetTitle()) + TString(" final asymmetry of ") 
		      + BdkEvtTypes::name(t),
		    TString(_typeAsym[t]->GetName()) + "+"
		      + TString(_globalAsym->GetName()),
		    RooArgList(*_typeAsym[t], *_globalAsym));
    // The # of events is N+ = N/2(1-A),  N- = N/2(1+A),  
    _numEvtPos[t] = new 
      RooFormulaVar(TString(GetName()) + TString(".numEvtPos")
		      + BdkEvtTypes::name(t),
		    TString(GetTitle()) + TString(" numEvts positive ")
		      + BdkEvtTypes::name(t),
		    TString(_numEvt[t]->GetName())
		      + TString("/2*(1-") 
		      + TString(_asym[t]->GetName())
		      + TString(")"),
		    RooArgList(*_numEvt[t], *_asym[t]));

    _numEvtNeg[t] = new 
      RooFormulaVar(TString(GetName()) + TString(".numEvtNeg")
		      + BdkEvtTypes::name(t),
		    TString(GetTitle()) + TString(" numEvts negative ")
		      + BdkEvtTypes::name(t),
		    TString(_numEvt[t]->GetName())
		      + TString("/2*(1+") 
		      + TString(_asym[t]->GetName())
		      + TString(")"),
		    RooArgList(*_numEvt[t], *_asym[t]));
  }

  // add the components:
  removeComponents();
  for (int c = 0; c < BdkEvtTypes::NTYPES; ++c) {
    addComponent(*prodP(c));
    addComponent(*prodN(c));
  }

  // number of BB events
  _nBB = new RooRealVar(TString(GetName()) + ".nBB",
                        TString(GetTitle()) + " nBB",
                        0);

  // absolute (cut&count) efficiency
  _absoluteEff = new RooRealVar(TString(GetName()) + ".absoluteEff",
                                TString(GetTitle()) + " absoluteEff",
                                0);

  // B -> DK branching fraction
  _brBtoDK = new RooRealVar(TString(GetName()) + ".brBtoDK",
                            TString(GetTitle()) + " brBtoDK",
                            0);

  // D -> K pi pi branching fraction
  _brDtoK2pi = new RooRealVar(TString(GetName()) + ".brDtoK2pi",
                              TString(GetTitle()) + " brDtoK2pi",
                              0);

  // (D -> 3pi) / (D -> K 2pi) branching fraction ratio
  _brDto3piOverK2pi = new RooRealVar(TString(GetName()) + ".brDto3piOverK2pi",
                                     TString(GetTitle()) + " brDto3piOverK2pi",
                                     0);

  setIsValid(kFALSE);
}

//------------------------------------------
RooArgSet BdkPdfOnRes::dependents()
{
  // make sure PDF is valid (e.g. after changing the components)
  if (getIsValid()==kFALSE) getPdf();

  return BdkPdfComposite::dependents();
}


//------------------------------------------
void BdkPdfOnRes::createPdf() {  
  int t;

  // init extended PDFs:
  for (t = 0; t < BdkEvtTypes::NTYPES; ++t) {
    _extPdfPos[t] = new 
      RooExtendPdf(TString(GetName()) 
		     + TString(".extPdfPos") 
		     + BdkEvtTypes::name(t),
		   TString(GetTitle()) 
		     + " Extended Pdf for positive "
		     + BdkEvtTypes::name(t),
		   *(prodP(t)->getPdf()), 
		   *_numEvtPos[t]); 
  
    _extPdfNeg[t] = new 
      RooExtendPdf(TString(GetName()) 
		     + TString(".extPdfNeg") 
		     + BdkEvtTypes::name(t),
		   TString(GetTitle()) 
		     + " Extended Pdf for negative "
		     + BdkEvtTypes::name(t),
		   *(prodN(t)->getPdf()), 
		   *_numEvtNeg[t]); 
  }  

  // Put the # of events variables and extended PDFs into ArgSets:
  RooArgList extPdfListPos;
  RooArgList evtNumsPos;

  RooArgList extPdfListNeg;
  RooArgList evtNumsNeg;

  for (t = 0; t < BdkEvtTypes::NTYPES; ++t) {    
    evtNumsPos.add(*_numEvtPos[t]);
    extPdfListPos.add(*_extPdfPos[t]);

    evtNumsNeg.add(*_numEvtNeg[t]);
    extPdfListNeg.add(*_extPdfNeg[t]);
  }
  
  // Make the 2 sum PDFs:
  _addPdfPos = 
    new RooAddPdf(TString(GetName())+".pdfPos", 
		  TString(GetTitle())+" pdf positive", 
		  extPdfListPos, evtNumsPos);

  _addPdfNeg = 
    new RooAddPdf(TString(GetName())+".pdfNeg", 
		  TString(GetTitle())+" pdf negative", 
		  extPdfListNeg, evtNumsNeg);

  // Make the simultaneous PDF:
  RooArgList simuPdfList(*_addPdfPos, *_addPdfNeg);

  RooSimultaneous * thePdf = new
    RooSimultaneous(TString(GetName()) + ".thePdf",
		    TString(GetTitle()) + " simultaneous pdf ",
		    simuPdfList,
		    *_charge);

  setPdf(*thePdf, kTRUE);


  // link additonal parameters to the pdf:
  RooListProxy * proxyList =  new RooListProxy(TString(GetName()) + ".proxyList",
                                               TString(GetTitle()) + " proxyList",
                                               thePdf);
  proxyList->add(*_nBB);
  proxyList->add(*_absoluteEff);
  proxyList->add(*_brBtoDK);
  proxyList->add(*_brDtoK2pi);
  proxyList->add(*_brDto3piOverK2pi);

  setIsValid(kTRUE);
}

//------------------------------------------
void BdkPdfOnRes::setMaxRange(int maxRange) {
  _maxRange = maxRange;
} 


//------------------------------------------
RooAbsReal * BdkPdfOnRes::numEvt(int index) {
  if (index < 0 || index >= BdkEvtTypes::NTYPES) {
    return 0;
  }

  return _numEvt[index];
}

//------------------------------------------
RooArgSet BdkPdfOnRes::numEvt() 
{
  RooArgSet set;
  for (Int_t i = 0; i < BdkEvtTypes::NTYPES; i++) 
    set.add(*_numEvt[i]);

  return set;
}

//------------------------------------------
RooAbsReal * BdkPdfOnRes::asym(int index) {
  if (index < 0 || index >= BdkEvtTypes::NTYPES) {
    return 0;
  }

  return _asym[index];
}
//----------------------------------------------
RooRealVar * BdkPdfOnRes::typeAsym(int index) {
  if (index < 0 || index >= BdkEvtTypes::NTYPES) {
    return 0;
  }

  return _typeAsym[index];
}

//------------------------------------------
RooFormulaVar * BdkPdfOnRes::numEvtPos(int index) {
  if (index < 0 || index >= BdkEvtTypes::NTYPES) {
    return 0;
  }

  return _numEvtPos[index];
}

//------------------------------------------
RooArgSet BdkPdfOnRes::numEvtPos() 
{
  RooArgSet set;
  for (Int_t i = 0; i < BdkEvtTypes::NTYPES; i++) 
    set.add(*_numEvtPos[i]);

  return set;
}

//------------------------------------------
RooFormulaVar * BdkPdfOnRes::numEvtNeg(int index) {
  if (index < 0 || index >= BdkEvtTypes::NTYPES) {
    return 0;
  }

  return _numEvtNeg[index];
}

//------------------------------------------
RooArgSet BdkPdfOnRes::numEvtNeg() 
{
  RooArgSet set;
  for (Int_t i = 0; i < BdkEvtTypes::NTYPES; i++) 
    set.add(*_numEvtNeg[i]);

  return set;
}

//------------------------------------------
RooExtendPdf * BdkPdfOnRes::extPdfPos(int index) {
  getPdf(); // wake up
   if (index < 0 || index >= BdkEvtTypes::NTYPES) {
    return 0;
  }

  return _extPdfPos[index];
}

//------------------------------------------
RooExtendPdf * BdkPdfOnRes::extPdfNeg(int index) {
  getPdf(); // wake up
   if (index < 0 || index >= BdkEvtTypes::NTYPES) {
    return 0;
  }

  return _extPdfNeg[index];
}

//------------------------------------------
RooDataSet * BdkPdfOnRes::generate(RooDataSet* protoDataSet) {
  return generate((int)totalNumEvts(), protoDataSet);
}

//------------------------------------------
RooDataSet * BdkPdfOnRes::generate(Int_t nEvents, RooDataSet* protoDataSet) {
  return generate(nEvents,false);
}

//------------------------------------------
RooDataSet * BdkPdfOnRes::generate(Int_t nEvents, Bool_t fluctuateNevt) {
  getPdf(); // make sure all are initialized
   RooDataSet * data = 0;
   RooDataSet * data1 = 0;
   Double_t SumPos = 0;
   Double_t SumNeg = 0;

   cout << GetName() << "::generate(): generating " << endl;
   for(Int_t i=0; i<BdkEvtTypes::NTYPES; i++) {
     cout << "  " << _numEvtPos[i]->getVal() 
	  << " positive and " << _numEvtNeg[i]->getVal() 
	  << " negative " 
	  << BdkEvtTypes::name(i) << " events." << endl;

     SumPos += _numEvtPos[i]->getVal();
     SumNeg += _numEvtNeg[i]->getVal();
   }

  Int_t nEvts1 = (Int_t) (nEvents*SumPos/(SumPos+SumNeg));
  Int_t nEvts2 = (Int_t) (nEvents*SumNeg/(SumPos+SumNeg));

  if (fluctuateNevt) {
    cout << GetName() <<".generate(): Requested number of events: B+ = "
         << nEvts1 << ", B- = " << nEvts2;

    nEvts1 = RooRandom::randomGenerator()->Poisson(nEvts1) ;
    nEvts2 = RooRandom::randomGenerator()->Poisson(nEvts2) ;

    cout << ". Generated: B+ = " << nEvts1 << ", B- = " << nEvts2 
         << ", asym = "<<(double)(nEvts2-nEvts1)/(nEvts1+nEvts2)      
         << endl;
  }

  RooArgSet deps = dependents(); 

  Bool_t genVerbose = kFALSE;

  _charge->setIndex(1);
  data = _addPdfPos->generate(deps, nEvts1, genVerbose);
  data->addColumn(*_charge);
  
  _charge->setIndex(-1); 
  data1 = _addPdfNeg->generate(deps, nEvts2, genVerbose);
  data1->addColumn(*_charge);
  data->append(*data1);

  delete data1;
  return data;
}

//------------------------------------------
double BdkPdfOnRes::totalNumEvts() const {
  double result = 0;
   for(Int_t i=0; i<BdkEvtTypes::NTYPES; ++i) {
     result += _numEvt[i]->getVal();
   }
   return result;
}

//---------------------------------------------------------
void BdkPdfOnRes::useNnCont(Bool_t use) {
  for (int i = 0; i < BdkEvtTypes::NTYPES; ++i) {
    prodP(i)->includeNnCont(use);
    prodN(i)->includeNnCont(use);
  }
}

//---------------------------------------------------------
void BdkPdfOnRes::useNnComb(Bool_t use) {
  for (int i = 0; i < BdkEvtTypes::NTYPES; ++i) {
    prodP(i)->includeNnComb(use);
    prodN(i)->includeNnComb(use);
  }
}

//---------------------------------------------------------
void BdkPdfOnRes::useDE(Bool_t use) {
  for (int i = 0; i < BdkEvtTypes::NTYPES; ++i) {
    prodP(i)->includeDeltaE(use);
    prodN(i)->includeDeltaE(use);
  }
}

//---------------------------------------------------------
void BdkPdfOnRes::useMd(Bool_t use) {
  for (int i = 0; i < BdkEvtTypes::NTYPES; ++i) {
    prodP(i)->includeMD(use);
    prodN(i)->includeMD(use);
  }
}

//---------------------------------------------------------
void BdkPdfOnRes::useMes(Bool_t use) {
  for (int i = 0; i < BdkEvtTypes::NTYPES; ++i) {
    prodP(i)->includeMes(use);
    prodN(i)->includeMes(use);
  }
}

//---------------------------------------------------------
void BdkPdfOnRes::useDalitz(Bool_t use) {
  for (int i = 0; i < BdkEvtTypes::NTYPES; ++i) {
    prodP(i)->includeDalitz(use);
    prodN(i)->includeDalitz(use);
  }
}

//---------------------------------------------------------
void BdkPdfOnRes::useVar(BdkPdfProdAll::Var var, Bool_t use) {
  switch (var) {
  case BdkPdfProdAll::DELTAE: 
    useDE(use);
    break;
  case BdkPdfProdAll::MD:
    useMd(use);
    break;
  case BdkPdfProdAll::MES:
    useMes(use);
    break;
  case BdkPdfProdAll::DALITZ:
    useDalitz(use);
    break;
  case BdkPdfProdAll::NNCONT:
    useNnCont(use);
    break;
  case BdkPdfProdAll::NNCOMB:
    useNnComb(use);
    break;
  }
}



//---------------------------------------------------------
void BdkPdfOnRes::printVariables() {
  // Print variables for arbitrary event type
  prodP(0)->printVariables();
}


//-----------------------------------------------
void BdkPdfOnRes::asymLink(int typeIndex, int sourceIndex) {

  if (typeIndex<0 || typeIndex>=BdkEvtTypes::NTYPES ||
      sourceIndex<0 || sourceIndex>=BdkEvtTypes::NTYPES) {
    cerr << GetName()<< ".asymLink(): Index " << typeIndex << " or "<<sourceIndex 
         << " is out of allowed range(0-9). Not linking asymmetries."
         << endl;
  }
  else {
    setIsValid(kFALSE);  
    _typeAsym[typeIndex] = (RooRealVar *) typeAsym(sourceIndex);
  }
}


//------------------------------------------------
BdkPdfProdAll * BdkPdfOnRes::prodP(int typeIndex) {
  switch(typeIndex) {
  case BdkEvtTypes::SIG_BAD_D   : return &sigBadD0P()   ; break;
  case BdkEvtTypes::SIG_GOOD_D  : return &sigGoodD0P()  ; break;
  case BdkEvtTypes::DPi_BAD_D   : return &DpiBadD0P()   ; break;
  case BdkEvtTypes::DPi_GOOD_D  : return &DpiGoodD0P()  ; break;
  case BdkEvtTypes::DPiX        : return &DPiXP()       ; break;
  case BdkEvtTypes::DKX         : return &DKXP()  ; break;
  case BdkEvtTypes::BB_BAD_D    : return &BBBadD0P() ; break;
  case BdkEvtTypes::BB_GOOD_D   : return &BBGoodD0P(); break;
  case BdkEvtTypes::QQ_BAD_D    : return &qqBadD0P()  ; break;
  case BdkEvtTypes::QQ_GOOD_D   : return &qqGoodD0P() ; break;

  }
  if( typeIndex<0 || typeIndex>=BdkEvtTypes::NTYPES ) {
   cerr<< " Index " << typeIndex << " is out of allwoed range(0-9)."
       << " Unrecognized argument \"" << typeIndex << "\". Returning 0."
       << endl;
 
  }

  return 0;
}

//------------------------------------------------
BdkPdfProdAll * BdkPdfOnRes::prodN(int typeIndex) {
  switch(typeIndex) {
  case BdkEvtTypes::SIG_BAD_D   : return &sigBadD0N()   ; break;
  case BdkEvtTypes::SIG_GOOD_D  : return &sigGoodD0N()  ; break;
  case BdkEvtTypes::DPi_BAD_D   : return &DpiBadD0N()   ; break;
  case BdkEvtTypes::DPi_GOOD_D  : return &DpiGoodD0N()  ; break;
  case BdkEvtTypes::DPiX        : return &DPiXN()       ; break;
  case BdkEvtTypes::DKX         : return &DKXN()  ; break;
  case BdkEvtTypes::BB_BAD_D    : return &BBBadD0N() ; break;
  case BdkEvtTypes::BB_GOOD_D   : return &BBGoodD0N(); break;
  case BdkEvtTypes::QQ_BAD_D    : return &qqBadD0N()  ; break;
  case BdkEvtTypes::QQ_GOOD_D   : return &qqGoodD0N() ; break;

  }
  if( typeIndex<0 || typeIndex>=BdkEvtTypes::NTYPES ) {
   cerr<< " Index " << typeIndex << " is out of allwoed range(0-9)."
       << " Unrecognized argument \"" << typeIndex << "\". Returning 0."
       << endl;
 
  }

  return 0;
}

//------------------------------------------------
// Set the values of the nsig and asym using the CP parameters:
void BdkPdfOnRes::setNsigAsymFromXY(Bool_t recalcErrors) 
{
  BdkOnResNLLYields yields(*this, 0);
  _sigGoodD0NumEvts->setVal(yields.expectedNsig());
  _typeAsym[BdkEvtTypes::SIG_GOOD_D]->setVal(yields.expectedAsym());

  if (recalcErrors) {
    double nsig = _sigGoodD0NumEvts->getVal();
    double asym = _typeAsym[BdkEvtTypes::SIG_GOOD_D]->getVal();

    double nsigErr = sqrt(nsig);
    double asymErr = sqrt(1. / nsig / (1 - asym * asym));

    _sigGoodD0NumEvts->setError(nsigErr);
    _typeAsym[BdkEvtTypes::SIG_GOOD_D]->setError(asymErr);
  }    

  cout << GetName() << "::setNsigAsymFromXY(): set"
       << " nBB=" << _nBB->getVal()
       << ", nSig=" << _sigGoodD0NumEvts->getVal()
       << " +/- " << _sigGoodD0NumEvts->getError()
       << ", asym=" << _typeAsym[BdkEvtTypes::SIG_GOOD_D]->getVal()
       << " +/- " << _typeAsym[BdkEvtTypes::SIG_GOOD_D]->getError()
       << " from"<<endl;
  cpParams().Print("v");
}


//------------------------------------------------
Bool_t BdkPdfOnRes::fit(RooAbsData& data, Bool_t usePenalty, Bool_t usePenaltyOnly)
{
  // Do the yield fit
  // Use N2 but no Dalitz
  useDE();
  useNnCont();
  useNnComb();
  useDalitz(false);

  cout << "----------------------------------------------------------"<<endl;
  cout << GetName() << ".fit(): Doing yield fit:" << endl;
  printVariables();
  parametersFree().Print("v");
  
  _yieldFitResult = getPdf()->fitTo(data,Extended(true),Save(true),
                                    Minos(_useMinos),PrintLevel(_minuitPrintLevel));
    
  if (_yieldFitResult==0) {
    cout << GetName() <<".fit(): No fit result from first fit."<<endl;
    _xyFitResult = 0;
    return kFALSE;
  }

  // Do the x/y fit
  return fit(data, *_yieldFitResult, usePenalty, usePenaltyOnly);
}


//------------------------------------------------
Bool_t BdkPdfOnRes::fit(RooAbsData& data, const RooFitResult& yieldFitResult,
                        Bool_t usePenalty, Bool_t usePenaltyOnly)
{
  // Get correlation between signal yield and asymmetry
  Double_t corrNsigAsym = yieldFitResult.correlation(*sigGoodD0NumEvts(),
                                                     *typeAsym(BdkEvtTypes::SIG_GOOD_D));

  return fit(data, corrNsigAsym, usePenalty, usePenaltyOnly);
}


//------------------------------------------------
Bool_t BdkPdfOnRes::fit(RooAbsData& data, Double_t corrNsigAsym,
                        Bool_t usePenalty, Bool_t usePenaltyOnly)
{
  TString fitOption = "er";
  if (!_useMinos) fitOption += "m";
  
  // Use Dalitz but no N2
  useNnComb(false);
  useDalitz();

  // Fix everything except x and y
  RooArgSet paramsFree = parametersFree();  // save to restore later
  fixAll();

  // If x/y (rho/theta) was free before make it floating again  
  RooArgList cpList(cpParams());

  for (int i=0; i<cpList.getSize(); i++) {
    RooRealVar& r = (RooRealVar&)cpList[i];
    if (paramsFree.contains(r)) r.setConstant(kFALSE);
  }

  // RooMinuit.cc line 833 gives error status if the NLL returns 0. So
  // move out of the minimum a little:
  ((RooRealVar&)cpList[0]).setVal(((RooRealVar&)cpList[0]).getVal() + 0.0001);

  RooNLLVar* nllVar = new RooNLLVar("nll","-log(likelihood)",*getPdf(),data,
                                    RooArgSet(),true);
  
  
  RooAbsReal* finalNLL = nllVar;
  BdkOnResNLLYields* penalty = 0;

  cout << "----------------------------------------------------------"<<endl;
  cout << GetName() << ".fit(): Doing x/y fit:" << endl;
  printVariables();
  nBB()->Print();
  cout << "corrNsigAsym = " << corrNsigAsym << endl;
  parametersFree().Print("v");


  if (usePenalty) {
    
    // Create penalty term
    if (_nllYieldsSystBit<0) penalty = new BdkOnResNLLYields(*this, corrNsigAsym);
    else penalty = new BdkOnResNLLYields(*this, corrNsigAsym, _nllYieldsSystBit);

    if (penalty==0) {
      cout << GetName() << ".fit(): Cannot create penalty term."<<endl;
      _yieldFitResult = 0;
      _xyFitResult = 0;
      return kFALSE;
    }

    penalty->Print();

    if (usePenaltyOnly) {    
      finalNLL = penalty;
      cout << GetName() << ".fit(): Fitting with penalty only" << endl;
    }    
    else {         
      finalNLL = new RooFormulaVar("nll+penalty", "nll+penalty", "@0+@1",
                                   RooArgList(*nllVar, *penalty));      
      cout << GetName() <<".fit(): Fitting with penalty+nll" << endl;
    }
  }
  else cout << GetName() <<".fit(): Fitting without penalty" << endl;
    
  _minuit = new RooMinuit(*finalNLL);
  _minuit->optimizeConst(1);   // const optimizer
  _minuit->setPrintLevel(_minuitPrintLevel);
  _minuit->fit(fitOption) ;
  
  _xyFitResult = _minuit->save() ;


  // Restore floating parameters
  TIterator* iter = paramsFree.createIterator();
  while (RooRealVar* r = (RooRealVar*)iter->Next()) {
    r->setConstant(kFALSE);
  }
  delete iter;
  
  return kTRUE;
}

// Return final fit parameters from both fits (caller owns)
RooArgSet* BdkPdfOnRes::fitResult()
{
  if (_yieldFitResult==0 && _xyFitResult==0) return 0;

  RooArgSet* result = new RooArgSet("fitResult");

  if (_yieldFitResult) {
    result->addClone(_yieldFitResult->floatParsFinal());
    result->addOwned(*(new RooRealVar("nllVar","nllVar",_yieldFitResult->minNll())));
  }
  if (_xyFitResult) {
    result->addClone(_xyFitResult->floatParsFinal());
    result->addOwned(*(new RooRealVar("xyNllVar","xyNllVar",_xyFitResult->minNll())));
  }

  return result;
}  


//------------------------------------------------
RooRealVar * BdkPdfOnRes::xPlus() {
  BdkPdfDKDalitz* pdf = (BdkPdfDKDalitz*)sigGoodD0P().getDalitzPdf();
  if (pdf->coord()==BdkPdfDKDalitz::CART) return (RooRealVar*)pdf->x();
  else return 0;
}

//------------------------------------------------
RooRealVar * BdkPdfOnRes::yPlus() {
  BdkPdfDKDalitz* pdf = (BdkPdfDKDalitz*)sigGoodD0P().getDalitzPdf();
  if (pdf->coord()==BdkPdfDKDalitz::CART) return (RooRealVar*)pdf->y();
  else return 0;
}

//------------------------------------------------
RooRealVar * BdkPdfOnRes::xMinus() {
  BdkPdfDKDalitz* pdf = (BdkPdfDKDalitz*)sigGoodD0N().getDalitzPdf();
  if (pdf->coord()==BdkPdfDKDalitz::CART) return (RooRealVar*)pdf->x();
  else return 0;
}

//------------------------------------------------
RooRealVar * BdkPdfOnRes::yMinus() {
  BdkPdfDKDalitz* pdf = (BdkPdfDKDalitz*)sigGoodD0N().getDalitzPdf();
  if (pdf->coord()==BdkPdfDKDalitz::CART) return (RooRealVar*)pdf->y();
  else return 0;
}



//------------------------------------------------
RooRealVar * BdkPdfOnRes::rhoPlus() {
  BdkPdfDKDalitz* pdf = (BdkPdfDKDalitz*)sigGoodD0P().getDalitzPdf();
  if (pdf->coord()==BdkPdfDKDalitz::POLAR) return (RooRealVar*)pdf->rho();
  else return 0;
}

//------------------------------------------------
RooRealVar * BdkPdfOnRes::thetaPlus() {
  BdkPdfDKDalitz* pdf = (BdkPdfDKDalitz*)sigGoodD0P().getDalitzPdf();
  if (pdf->coord()==BdkPdfDKDalitz::POLAR) return (RooRealVar*)pdf->theta();
  else return 0;
}

//------------------------------------------------
RooRealVar * BdkPdfOnRes::rhoMinus() {
  BdkPdfDKDalitz* pdf = (BdkPdfDKDalitz*)sigGoodD0N().getDalitzPdf();
  if (pdf->coord()==BdkPdfDKDalitz::POLAR) return (RooRealVar*)pdf->rho();
  else return 0;
}

//------------------------------------------------
RooRealVar * BdkPdfOnRes::thetaMinus() {
  BdkPdfDKDalitz* pdf = (BdkPdfDKDalitz*)sigGoodD0N().getDalitzPdf();
  if (pdf->coord()==BdkPdfDKDalitz::POLAR) return (RooRealVar*)pdf->theta();
  else return 0;
}



//------------------------------------------------
// Return the set of CP parameters
RooArgSet BdkPdfOnRes::cpParams()
{
  const BdkPdfDKDalitz& pdfN = *(BdkPdfDKDalitz*)sigGoodD0N().getDalitzPdf();
  const BdkPdfDKDalitz& pdfP = *(BdkPdfDKDalitz*)sigGoodD0P().getDalitzPdf();
  
  RooArgSet cpSet(GetName()+TString(".cpParams"));

  if (pdfN.coord()==BdkPdfDKDalitz::CART) {
    cpSet.add(*pdfN.x());
    cpSet.add(*pdfN.y());
  }
  else {
    cpSet.add(*pdfN.rho());
    cpSet.add(*pdfN.theta());
  }

  if (pdfP.coord()==BdkPdfDKDalitz::CART) {
    cpSet.add(*pdfP.x());
    cpSet.add(*pdfP.y());
  }
  else {
    cpSet.add(*pdfP.rho());
    cpSet.add(*pdfP.theta());
  }

  return cpSet;
}
