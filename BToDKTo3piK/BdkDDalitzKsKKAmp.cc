/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkDDalitzKsKKAmp.cc,v 1.3 2008/01/02 20:30:09 fwinkl Exp $
 * Description:
 *   Implementation of D0 -> K+K-Ks Dalitz amplitude
 * History:
 *   10 Apr 2007, created, Frank Winklmeier
 *
 * Copyright (C) 2007 Colorado State University and SLAC
 *****************************************************************************/


#include <iostream>

#include "BToDKTo3piK/BdkDDalitzKsKKAmp.hh"
#include "BToDKTo3piK/BdkDalitzBase.hh"
#include "BToDKTo3piK/BdkAmp.hh"

using namespace std;


ClassImp(BdkDDalitzKsKKAmp);


// Constructor
BdkDDalitzKsKKAmp::BdkDDalitzKsKKAmp(const char * name, 
                                     const char * title, 
                                     BdkDalitzBase * pdf,
                                     Int_t componentsBit) :
  BdkAbsDDalitzAmp(name, title, pdf, componentsBit)
{
  initResonance();
  createParams();
  registerParams(_pdf);
}


// Destructor
BdkDDalitzKsKKAmp::~BdkDDalitzKsKKAmp() 
{
}

Double_t BdkDDalitzKsKKAmp::efficiency(Double_t m12, Double_t m13) const
{
  //return _pdf->efficiency(m12, m13);
  return 1;
}


// initialize the components:
void BdkDDalitzKsKKAmp::initResonance() 
{

  if (nComps()!=0) {
    cout << "BdkDDalitzKsKKAmp::initResonance(): Cannot call initResonance more than once!" << endl;
    return;
  }
  
  // Define all RooRealVar floating observables.

  /* The trackinfo variable defines the cyclycal nature of the resonances
     and the final state pions, with the following order. 
     See enum ResDaughters.

     index   particle
     1       Ks
     2       K+
     3       K-
     --------------
     index   resonance   dtr particles   dtr indices
     1                   K+ K-         23
     2                   Ks K-         13
     3                   Ks K+         12
  */

  int a0PIndex = -1;

  // amplitudes and phases from hep-ex/0507026
  // We don't use a coupled channel model for the a0/f0!
  // Masses and widths are estimated from pdg listings.
  
  if (A0P & componentsBit()) {
    a0PIndex = nComps();  // store its index
    addComp("a0+", 0.46, 3.59/DEGTORAD, 0.981, 0.092, PI0_PIP, SPIN0, -1, A0RES);
  }

  if (A00 & componentsBit()) {
    addComp("a00", 1.0, 0.0, 0.981, 0.092, PIP_PIM, SPIN0, a0PIndex, A0RES);
  }

  if (F0_1370 & componentsBit()) {
    addComp("F0_1370", 0.435, -2.63/DEGTORAD, 1.400, 0.150, PIP_PIM, SPIN0);
  }

  if (PHI & componentsBit()) {
    addComp("Phi", 0.437, 1.91/DEGTORAD, 1.0195, 0.00426, PIP_PIM, SPIN1);
  }

  
  // Masses of the above daughters
  _mDaug[0] = 0.0;
  _mDaug[1] = BdkDalitz::K0MASS;
  _mDaug[2] = BdkDalitz::KMASS;
  _mDaug[3] = BdkDalitz::KMASS;

  // Parameters for a0 coupled channel BW
  _a0gEtaPi = new RooRealVar(TString(GetName()) + ".a0_gEtaPi",
                             TString(GetTitle()) + " a0_gEtaPi^2",
                             0.324);
  _params.addOwned(*_a0gEtaPi);

  _a0gKK = new RooRealVar(TString(GetName()) + ".a0_gKK",
                          TString(GetTitle()) + " a0_gKK^2",
                          0.329);
                          _params.addOwned(*_a0gKK);
}


RooComplex BdkDDalitzKsKKAmp::matrixElement(Double_t m12, Double_t m13, Int_t i) const
{   
  if (_typeRes[i]==A0RES) {
    Double_t s = kinematics(m12,m13,_trackinfo[i]);
    return BdkAmp::a0Amp(s, _massRes[i]->getVal(),
                         _a0gEtaPi->getVal(), _a0gKK->getVal());
  }
  else return BdkAbsDDalitzAmp::matrixElement(m12, m13, i);
}

