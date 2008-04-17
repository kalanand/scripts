 /****************************************************************************
 * File: $Id: BdkDKDalitz.cc,v 1.20 2007/04/18 11:57:05 fwinkl Exp $
 *
 * History:
 * Sep 19 2005, Abi sOffer, adapted from Ben Lau's RooDKminus
 *                                                                           *
 *****************************************************************************/
          
// -- CLASS DESCRIPTION [PDF] --
// Here we declared the B- -> D0 K- pdf
// read hep-ex/0308043 equation (2)
// M- = f(m-^2,m+^2)+rb Exp(i(-\gamma_+\delta_) * f(m+^2,m-^2)

//#include "BaBar/BaBar.hh"

#include <iostream>
#include <fstream>
#include <math.h>


#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooRandom.hh"
#include "RooFitCore/RooComplex.hh"
#include "RooFitCore/RooRealProxy.hh"
#include "BToDKTo3piK/BdkDKDalitz.hh"
#include "BToDKTo3piK/BdkDDalitzAmp.hh"

using namespace std;



ClassImp(BdkDKDalitz)

// constructor
BdkDKDalitz::BdkDKDalitz(const char *theName, const char *theTitle,
                         RooAbsReal& m12, RooAbsReal& m13,
                         RooAbsReal& x, RooAbsReal& y,
                         BdkDalitzBase::Flavor flavor, 
                         BdkDalitzBase::Mode DdecMode,
                         BdkAbsDDalitzAmp* amp):
  BdkDalitzBase(theName, theTitle, flavor, DdecMode),
  _m12("m12","Invariant Mass square of pi0 pi+",this,m12),
  _m13("m13","Invariant Mass square of pi0 pi-",this,m13),
  _x("x","Real part of rB*phase",this,x),
  _y("y","Imaginary part of rB*phase",this,y)
{
  if (amp==0) _dalitzAmp = new BdkDDalitzAmp(TString(GetName()) + ".dalitzAmp",
                                             TString(GetTitle()) + " dalitzAmp",
                                             this);
  else {
    _dalitzAmp = amp;
    // register all parameters with us
    _dalitzAmp->registerParams(this);
  }
}

//Copy Constructor

BdkDKDalitz::BdkDKDalitz(const BdkDKDalitz& other, const char* name) :
  BdkDalitzBase(other, name),
  _m12("m12",this,other._m12),
  _m13("m13",this,other._m13),
  _x("x",this,other._x),
  _y("y",this,other._y),
  _dalitzAmp(other._dalitzAmp)
{
  // register all parameters with us
  _dalitzAmp->registerParams(this);
}

BdkDKDalitz::~BdkDKDalitz()
{
}

Double_t BdkDKDalitz::evaluate() const
{
  // these are m12 and m13 for B-:
  double m12temp = _m12;
  double m13temp = _m13;

  // flip for B+:
  if (BdkDalitzBase::D0BAR == flavor()) {
    m12temp = _m13;
    m13temp = _m12;
  }

  RooComplex Damp    =  _dalitzAmp->getamp(m12temp, m13temp);
  RooComplex Dbaramp =  _dalitzAmp->getamp(m13temp, m12temp);

  //  RooComplex c(_rB*cos(phaseRad()),_rB*sin(phaseRad()));
  RooComplex c(_x,_y);

  _lastAmp = Damp + c*(Dbaramp);
  return _lastAmp.abs2();
}


Int_t BdkDKDalitz::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  if (matchArgs(allVars,analVars,_m12,_m13)) return 1 ;
  return 0 ;
}



Double_t BdkDKDalitz::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1) ;
  return getNormalization();
} 


Double_t BdkDKDalitz::getNormalization() const
{
  //  RooComplex c(_x,_y);
  //  Double_t norm = _p3.getVal()+(c.abs2())*_p2.getVal()+
  //                   2.0*(_p0.getVal()*(c.re())+_p1.getVal()*(c.im()));

  Double_t norm = _dalitzAmp->normDSqr() +
                  (_x*_x + _y*_y) * _dalitzAmp->normDbarSqr() +
                  2*_x*_dalitzAmp->normReDDbar() + 
                  2*_y*_dalitzAmp->normImDDbar();

  return norm;
}

