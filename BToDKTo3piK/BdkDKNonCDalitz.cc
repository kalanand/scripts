/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkDKNonCDalitz.cc,v 1.1 2007/04/18 11:57:05 fwinkl Exp $
 * Description:
 *    A PDF with two non C eigenstate Dalitz amplitudes
 * History:
 *   11 Apr 2007, created, Frank Winklmeier
 *
 * Copyright (C) 2007 Colorado State University and SLAC
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>


#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooRandom.hh"
#include "RooFitCore/RooComplex.hh"
#include "RooFitCore/RooRealProxy.hh"
#include "BToDKTo3piK/BdkDKNonCDalitz.hh"
#include "BToDKTo3piK/BdkDDalitzAmp.hh"

using namespace std;



ClassImp(BdkDKNonCDalitz)

// constructor
BdkDKNonCDalitz::BdkDKNonCDalitz(const char *theName, const char *theTitle,
                                 RooAbsReal& m12, RooAbsReal& m13,
                                 RooAbsReal& x, RooAbsReal& y, RooAbsReal& deltaD,
                                 BdkDalitzBase::Flavor flavor, 
                                 BdkDalitzBase::Mode DdecMode,
                                 BdkAbsDDalitzAmp* amp,
                                 BdkAbsDDalitzAmp* ampBar):
  BdkDalitzBase(theName, theTitle, flavor, DdecMode),
  _m12("m12","Invariant Mass square of pi0 pi+",this,m12),
  _m13("m13","Invariant Mass square of pi0 pi-",this,m13),
  _x("x","Real part of rB*phase",this,x),
  _y("y","Imaginary part of rB*phase",this,y),
  _deltaD("deltaD","Phase between D/Dbar amplitudes",this,deltaD),
  _dalitzAmp(amp),
  _dalitzAmpBar(ampBar)
{
  if (_dalitzAmp==0 || _dalitzAmpBar==0) {
    cout << GetName() << ": You need to supply two non-zero BdkAbsDDaltizAmp objects!"
         << endl;
    return;
  }
      
  // register all parameters with us
  _dalitzAmp->registerParams(this);
  _dalitzAmpBar->registerParams(this);  
}

//Copy Constructor

BdkDKNonCDalitz::BdkDKNonCDalitz(const BdkDKNonCDalitz& other, const char* name) :
  BdkDalitzBase(other, name),
  _m12("m12",this,other._m12),
  _m13("m13",this,other._m13),
  _x("x",this,other._x),
  _y("y",this,other._y),
  _deltaD("deltaD",this,other._deltaD),
  _dalitzAmp(other._dalitzAmp),
  _dalitzAmpBar(other._dalitzAmpBar)
{
  // register all parameters with us
  if (_dalitzAmp) _dalitzAmp->registerParams(this);
  if (_dalitzAmpBar) _dalitzAmpBar->registerParams(this);
}

BdkDKNonCDalitz::~BdkDKNonCDalitz()
{
}

Double_t BdkDKNonCDalitz::evaluate() const
{
  // this is for B-:
  BdkAbsDDalitzAmp* amp = _dalitzAmp;
  BdkAbsDDalitzAmp* ampBar = _dalitzAmpBar;

  // flip for B+:
  if (BdkDalitzBase::D0BAR == flavor()) {
    amp = _dalitzAmpBar;
    ampBar = _dalitzAmp;
  }
  
  RooComplex Damp    =  amp->getamp(_m12, _m13);
  RooComplex Dbaramp =  ampBar->getamp(_m12, _m13);

  RooComplex c(_x,_y+_deltaD);
  RooComplex result = Damp + c*(Dbaramp);

  return result.abs2();
}


Int_t BdkDKNonCDalitz::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  if (matchArgs(allVars,analVars,_m12,_m13)) return 1 ;
  return 0 ;
}


Double_t BdkDKNonCDalitz::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1) ;

  return getNormalization();
} 


Double_t BdkDKNonCDalitz::getNormalization() const
{
  Double_t norm = _dalitzAmp->normDSqr() +
                  (_x*_x + _y*_y) * _dalitzAmp->normDbarSqr() +
                  2*_x*_dalitzAmp->normReDDbar() + 
                  2*_y*_dalitzAmp->normImDDbar();
   
  return norm;
}


void BdkDKNonCDalitz::calDDbarNorm(int nEvents)
{
  cout << GetName() 
       << ": performing MC integration for D/Dbar interference" << endl;

  TVectorD p = _dalitzAmp->calDIntNorm(_dalitzAmpBar, nEvents);
  
  cout << "Precision of MC integration:"<<endl;
  cout << "p1/p0      = "<< p[1]/p[0] << endl;
  cout << "(p2-p3)/p2 = "<< (p[2]-p[3])/p[2] << endl;
  
  // Copy to RooRealVars of _dalitzAmp (we are a friend of BdkAbsDDalitzAmp)
  _dalitzAmp->_normReDDbar->setVal(p[0]);      
  _dalitzAmp->_normImDDbar->setVal(p[1]);         // Theoretically, p1 = 0
  _dalitzAmp->_normDbarSqr->setVal(p[2]); // (p[2]+p[3])/2);   // Theoretically, p2 = p3
  _dalitzAmp->_normDSqr->setVal(p[3]); //(p[2]+p[3])/2);
}
