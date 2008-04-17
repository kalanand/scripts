/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkDDalitzKsPiPiAmp.cc,v 1.6 2008/01/02 20:30:10 fwinkl Exp $
 * Description:
 *   Implementation of D0 -> pi+pi-Ks Dalitz amplitude
 * History:
 *   04 Apr 2007, created, Frank Winklmeier
 *
 * Copyright (C) 2007 Colorado State University and SLAC
 *****************************************************************************/


#include <iostream>

#include "BToDKTo3piK/BdkDDalitzKsPiPiAmp.hh"
#include "BToDKTo3piK/BdkDalitzBase.hh"
#include "BToDKTo3piK/BdkAmp.hh"

using namespace std;


ClassImp(BdkDDalitzKsPiPiAmp);


// Constructor
BdkDDalitzKsPiPiAmp::BdkDDalitzKsPiPiAmp(const char * name, 
                                         const char * title, 
                                         BdkDalitzBase * pdf,
                                         Int_t componentsBit,
                                         Int_t defaultRhoSpin) :
  BdkAbsDDalitzAmp(name, title, pdf, componentsBit),
  _defaultRhoSpin(defaultRhoSpin)
{
  initResonance();
  createParams();
  registerParams(_pdf);
}


// Destructor
BdkDDalitzKsPiPiAmp::~BdkDDalitzKsPiPiAmp() 
{
}

Double_t BdkDDalitzKsPiPiAmp::efficiency(Double_t m12, Double_t m13) const
{
  //return _pdf->efficiency(m12, m13);
  return 1;
}


// initialize the components:
void BdkDDalitzKsPiPiAmp::initResonance() 
{

  if (nComps()!=0) {
    cout << "BdkDDalitzKsPiPiAmp::initResonance(): Cannot call initResonance more than once!" << endl;
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

  int kstPIndex = -1;
  int kst0_1430PIndex = -1;
  int kst2_1430PIndex = -1;
  int nonresIndex = -1;
  int kst_1680PIndex = -1;
  
  // pi+pi- resonances:
  if (RHO0 & componentsBit()) {
    addComp("Rho0", 0.59, 10.0, 0.7758, 0.1503, PIP_PIM, _defaultRhoSpin,-1,GSRES);
  }

  if (RHO2S0 & componentsBit()) {
    addComp("Rho2s0", 0.0, 0.0, 1.465, 0.400, PIP_PIM, _defaultRhoSpin,-1,GSRES);
  }

  if (SIGMA & componentsBit()) {
    addComp("Sigma", 0.0, 0.0, 0.500, 0.400, PIP_PIM, SPIN0);
  }

  if (SIGMA2 & componentsBit()) {
    addComp("Sigma2", 0.0, 0.0, 0.500, 0.400, PIP_PIM, SPIN0);
  }
 
  if (OMEGA & componentsBit()) {
    addComp("Omega", 0.0, 0.0, 0.78268, 0.00849, PIP_PIM, SPIN1);
  }

  if (F0 & componentsBit()) {
    addComp("F0", 0.0, 0.0, 0.980, 0.044, PIP_PIM, SPIN0);

    // width for the F0 is taken from PRL 86, 765 (E791, 2001).
    // which looked at Ds->pi+pi-pi+
  }

  if (F0_1370 & componentsBit()) {
    addComp("F0_1370", 0.0, 0.0, 1.434, 0.173, PIP_PIM, SPIN0);

    // parameters for the F0_1370 are taken from PRL 86, 765 (E791, 2001).
    // which looked at Ds->pi+pi-pi+
  }

  if (F2 & componentsBit()) {
    addComp("F2", 0.0, 0.0, 1.2754, 0.1851, PIP_PIM, SPIN2);
  }

  // Ks pi+ resonances:
  if (KSTP & componentsBit()) {
    kstPIndex = nComps();
    addComp("Kst+", 0.0, 0.0, 0.8917, 0.0508, PI0_PIP, SPIN1);
  }

  if (KST0_1430P & componentsBit()) {
    kst0_1430PIndex = nComps();
    addComp("Kst0_1430+", 0.0, 0.0, 1.414, 0.290, PI0_PIP, SPIN0);
  }

  if (KST2_1430P & componentsBit()) {
    kst2_1430PIndex = nComps();
    addComp("Kst2_1430+", 0.0, 0.0, 1.425, 0.099, PI0_PIP, SPIN2);
  }

  if (KST_1680P & componentsBit()) {
    kst_1680PIndex = nComps();
    addComp("Kst_1680+", 0.0, 0.0, 1.717, 0.322, PIM_PI0, SPIN1);
  }
  // Ks pi- resonances:
  if (KSTM & componentsBit()) {
    addComp("Kst-", 0.0, 0.0, 0.8917, 0.0508, PIM_PI0, SPIN1,
            kstPIndex);
  }

  if (KST0_1430M & componentsBit()) {
    addComp("Kst0_1430-", 0.0, 0.0, 1.414, 0.290, PIM_PI0, SPIN0,
            kst0_1430PIndex);
  }

  if (KST2_1430M & componentsBit()) {
    addComp("Kst2_1430-", 0.0, 0.0, 1.425, 0.099, PIM_PI0, SPIN2,
            kst2_1430PIndex);
  }

  if (KST_1410M & componentsBit()) {
    addComp("Kst_1410-", 0.0, 0.0, 1.414, 0.232, PIM_PI0, SPIN1);
  }

  if (KST_1680M & componentsBit()) {
    addComp("Kst_1680-", 0.0, 0.0, 1.717, 0.322, PIM_PI0, SPIN1,
            kst_1680PIndex);
  }

  
  // Non-resonant:
  if (NONRES & componentsBit()) {
    nonresIndex = nComps();  // store its index
    addComp("Nonres", 1.03, 77, 0, 0, PIP_PIM_PI0, SPIN0);
  }
   
  // Masses of the above daughters
  _mDaug[0] = 0.0;
  _mDaug[1] = BdkDalitz::K0MASS;
  _mDaug[2] = BdkDalitz::PIMASS;
  _mDaug[3] = BdkDalitz::PIMASS;
}


RooComplex BdkDDalitzKsPiPiAmp::matrixElement(Double_t m12, Double_t m13, Int_t i) const
{ 
  if (_typeRes[i]==GSRES) {   //Gounaris-Sakurai for the rhos

    Double_t m0 = _massRes[i]->getVal();
    Double_t k0 = pionCMmom(m0*m0,_trackinfo[i]);
    Double_t s = kinematics(m12,m13,_trackinfo[i]);
    Double_t k = pionCMmom(s,_trackinfo[i]);
    Double_t g = runningWidth(m12,m13,m0,
			      _gammaRes[i]->getVal(),
			      _spinRes[i],
			      _trackinfo[i]);

    return BdkAmp::GS_Rho(m0, _gammaRes[i]->getVal(), k0, 
			  s, g, k) 
      * SpinFactor(m12, m13, _spinRes[i], _trackinfo[i], i);

  }
  else  // regular Breit-Wigner
    return BdkAbsDDalitzAmp::matrixElement(m12, m13, i);
}
