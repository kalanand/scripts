/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkDDalitzKsKPiAmp.cc,v 1.4 2008/01/03 18:11:16 fwinkl Exp $
 * Description:
 *   Implementation of D -> Ks K- pi+ Dalitz amplitude
 *   Can also handle D -> Ks K+ pi- (piCharge in constructor)
 * History:
 *   11 Apr 2007, created, Frank Winklmeier
 *
 * Copyright (C) 2007 Colorado State University and SLAC
 *****************************************************************************/

#include <iostream>

#include "BToDKTo3piK/BdkDDalitzKsKPiAmp.hh"
#include "BToDKTo3piK/BdkDalitzBase.hh"
#include "BToDKTo3piK/BdkAmp.hh"

using namespace std;


ClassImp(BdkDDalitzKsKPiAmp);


// Constructor
BdkDDalitzKsKPiAmp::BdkDDalitzKsKPiAmp(const char * name, 
                                       const char * title, 
                                       BdkDalitzBase * pdf,
                                       Int_t componentsBit,
                                       Int_t piCharge) :
  BdkAbsDDalitzAmp(name, title, pdf, componentsBit),
  _piCharge(piCharge)
{
  initResonance();
  createParams();
  registerParams(_pdf);
}


// Destructor
BdkDDalitzKsKPiAmp::~BdkDDalitzKsKPiAmp() 
{
}

Double_t BdkDDalitzKsKPiAmp::efficiency(Double_t m12, Double_t m13) const
{
  //return _pdf->efficiency(m12, m13);
  return 1;
}


void BdkDDalitzKsKPiAmp::calDDbarNorm(int nEvents)
{
  cout << GetName() << " does not implement calDDbarNorm() because it represents "
       << "a non C eigenstate. Use the normalization method of the appropriate "
       << "RooAbsPdf derived classs."
       << endl;
}

// initialize the components:
void BdkDDalitzKsKPiAmp::initResonance() 
{

  if (nComps()!=0) {
    cout << "BdkDDalitzKsKPiAmp::initResonance(): Cannot call initResonance more than once!" << endl;
    return;
  }
  
  // Define all RooRealVar floating observables.

  /* The trackinfo variable defines the cyclycal nature of the resonances
     and the final state pions, with the following order. 
     See enum ResDaughters.

     index   particle
     1       Ks
     2       pi+
     3       K-
     --------------
     index   resonance   dtr particles   dtr indices
     1                   pi+ K-          23
     2                   Ks  K-          13
     3                   Ks  pi+         12
  */

  int kstPIndex = -1;
  int kst0_1430PIndex = -1;
  int kst2_1430PIndex = -1;
  int kst1_1680PIndex = -1;

  // The same components are needed for D0 -> Ks K+ pi-
  TString sign, signBar;
  ResDaughters daug, daugBar;
  if (_piCharge<0) {
    sign = "-";
    signBar = "+";
    daug = PIM_PI0;
    daugBar = PI0_PIP;    
  }
  else {
    sign = "+";
    signBar = "-";
    daug = PI0_PIP;
    daugBar = PIM_PI0;
  }

  // Masses of the above daughters
  _mDaug[0] = 0.0;
  _mDaug[1] = BdkDalitz::K0MASS;
  _mDaug[2] = (_piCharge>=0 ? BdkDalitz::PIMASS : BdkDalitz::KMASS);
  _mDaug[3] = (_piCharge>=0 ? BdkDalitz::KMASS : BdkDalitz::PIMASS);
  
  // Ks pi+ resonances:
  if (KSTP & componentsBit()) {
    kstPIndex = nComps();
    addComp("Kst"+sign, 0.0, 0.0, 0.8917, 0.0508, daug, SPIN1);
  }

  if (KST0_1430P & componentsBit()) {
    kst0_1430PIndex = nComps();
    addComp("Kst0_1430"+sign, 0.0, 0.0, 1.414, 0.290, daug, SPIN0);
  }

  if (KST2_1430P & componentsBit()) {
    kst2_1430PIndex = nComps();
    addComp("Kst2_1430"+sign, 0.0, 0.0, 1.425, 0.099, daug, SPIN2);
  }

  if (KST1_1680P & componentsBit()) {
    kst1_1680PIndex = nComps();
    addComp("Kst1_1680"+sign, 0.0, 0.0, 1.717, 0.322, daug, SPIN1);
  }

  // K- pi+ resonances:
  if (KST & componentsBit()) {
    addComp("Kst", 0.0, 0.0, 0.8917, 0.0508, PIP_PIM, SPIN1, kstPIndex);
  }

  if (KST0_1430 & componentsBit()) {
    addComp("Kst0_1430", 0.0, 0.0, 1.414, 0.290, PIP_PIM, SPIN0, kst0_1430PIndex);
  }

  if (KST2_1430 & componentsBit()) {
    addComp("Kst2_1430", 0.0, 0.0, 1.425, 0.099, PIP_PIM, SPIN2, kst2_1430PIndex);
  }

  if (KST1_1680 & componentsBit()) {
    addComp("Kst1_1680", 0.0, 0.0, 1.717, 0.322, PIP_PIM, SPIN1, kst1_1680PIndex);
  }
  
  // Ks K- resonances:
  if (A0M & componentsBit()) {
    addComp("a0"+signBar, 0.0, 0.0, 0.9847, 0.092, daugBar, SPIN0, -1, A0RES);
  }

  if (A0_1450M & componentsBit()) {
    addComp("a0_1450"+signBar, 0.0, 0.0, 1.474, 0.265, daugBar, SPIN0);
  }

  if (A2_1310M & componentsBit()) {
    addComp("a2_1310"+signBar, 0.0, 0.0, 1.318, 0.107, daugBar, SPIN2);
  } 
    
   
  // Non-resonant:
  if (NONRES & componentsBit()) {
    addComp("Nonres", 1.03, 77, 0, 0, PIP_PIM_PI0, SPIN0);
  }

  // Parameters for a0 coupled channel BW
  _a0gEtaPi = new RooRealVar(TString(GetName()) + ".a0-_gEtaPi",
                             TString(GetTitle()) + " a0-_gEtaPi^2",
                             0.324);
  _params.addOwned(*_a0gEtaPi);

  _a0gKK = new RooRealVar(TString(GetName()) + ".a0-_gKK",
                          TString(GetTitle()) + " a0-_gKK^2",
                          0.329);
  _params.addOwned(*_a0gKK);
}


RooComplex BdkDDalitzKsKPiAmp::matrixElement(Double_t m12, Double_t m13, Int_t i) const
{   
  if (_typeRes[i]==A0RES) {
    Double_t s = kinematics(m12,m13,_trackinfo[i]);
    return BdkAmp::a0AmpAntimo(s, _massRes[i]->getVal(),
			       _a0gEtaPi->getVal(), _a0gKK->getVal());
  }
  else return BdkAbsDDalitzAmp::matrixElement(m12, m13, i);
}
