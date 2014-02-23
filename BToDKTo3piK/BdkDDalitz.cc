/*****************************************************************************
 * Project: BTodKTo3piK
 *    File: $Id: BdkDDalitz.cc,v 1.13 2007/04/18 11:57:04 fwinkl Exp $
 *
 * RooAbsPdf subclass that fits the Dalitz distribution of D decays.
 *
 * History: Abi Soffer, Oct 31, 2005, adapted from code by Ben Lau and 
 *    Kalanand Mishra
 *                                                                           *
 *****************************************************************************/

#include <iostream>
#include <math.h>
#include "BToDKTo3piK/BdkDDalitz.hh"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooRandom.hh"
#include "RooFitCore/RooComplex.hh"
#include "RooFitCore/RooRealProxy.hh"
#include "RooFitCore/RooListProxy.hh"
#include "RooFitCore/RooArgList.hh"


using namespace std;


ClassImp(BdkDDalitz)

/// Define our constructor, we need m12^2 and m13^2, as our dalitz variable
BdkDDalitz::BdkDDalitz(const char *name, const char *title,
                       RooAbsReal& m12, RooAbsReal& m13,
                       BdkDalitzBase::Flavor flavor, 
                       BdkAbsDDalitzAmp * amp,
                       Int_t comps,
                       BdkDalitzBase::Mode DdecMode,
                       Int_t spinResComp):
  BdkDalitzBase(name, title, flavor, DdecMode),
  _m12("m12","Invariant Mass square of M12",this,m12),
  _m13("m13","Invariant Mass square of M13",this,m13)
{
  setComponents(comps, spinResComp, amp);
}


/// Copy Constructor
BdkDDalitz::BdkDDalitz(const BdkDDalitz& other, const char* name) :
  BdkDalitzBase(other, name),
  _m12("m12",this,other._m12),
  _m13("m13",this,other._m13),
  _dalitzAmp(other._dalitzAmp)
{
  // update the dalitzAmp's inDalitz calculator to this:
  _dalitzAmp->registerParams(this);  
}

/// Destructor
BdkDDalitz::~BdkDDalitz()
{
  // Can't delete the _dalitzAmp since the copy constructor doesn't properly copy
  //delete _dalitzAmp;
}

/// reinitialize the _dalitzAmp with new components:
void BdkDDalitz::setComponents(Int_t comp, Int_t spinResComp,
                               BdkAbsDDalitzAmp * amp)
{
  if (0 == amp) { // make our own:
    _dalitzAmp = new BdkDDalitzAmp(TString(GetName()) + ".dalitzAmp",
				   TString(GetTitle()) + " dalitzAmp",
				   this,
				   comp, spinResComp);
  }
  else { // use the one given:
    _dalitzAmp = amp;
    _dalitzAmp->registerParams(this); // register the parameters for this PDF
  }
}

Double_t BdkDDalitz::evaluate() const
{
  RooComplex Dbaramp;

  if (BdkDalitzBase::D0BAR == flavor()) {
    Dbaramp = _dalitzAmp->getamp(_m13,_m12); // flip for D0bar
  }
  else {
    Dbaramp = _dalitzAmp->getamp(_m12,_m13); // for D0
  }
  double result = Dbaramp.abs2();
  return result;
}


Int_t BdkDDalitz::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  if (matchArgs(allVars,analVars,_m12,_m13)) return 1 ;
  return 0 ;
}


///the normalization is calculated in _dalitzAmp
Double_t BdkDDalitz::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1);
  Double_t norm = _dalitzAmp->getNormalization();
  return norm;
} 




