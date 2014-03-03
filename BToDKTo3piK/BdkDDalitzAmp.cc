/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkDDalitzAmp.cc,v 1.42 2007/04/05 15:47:34 fwinkl Exp $
 * Description:
 *   Implementation of D0 -> pi+pi-pi0 Dalitz amplitude
 * History:
 *   04 Apr 2007, created, Frank Winklmeier
 *
 * Copyright (C) 2007 Colorado State University and SLAC
 *****************************************************************************/

#include <iostream>

#include "BToDKTo3piK/BdkDDalitzAmp.hh"
#include "BToDKTo3piK/BdkDalitzBase.hh"

using namespace std;


ClassImp(BdkDDalitzAmp);


// Constructor
BdkDDalitzAmp::BdkDDalitzAmp(const char * name, 
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
BdkDDalitzAmp::~BdkDDalitzAmp() 
{
}

Double_t BdkDDalitzAmp::efficiency(Double_t m12, Double_t m13) const
{
  return _pdf->efficiency(m12, m13);
}


// initialize the components:
void BdkDDalitzAmp::initResonance() 
{

  if (nComps()!=0) {
    cout << "BdkDDalitzAmp::initResonance(): Cannot call initResonance more than once!" << endl;
    return;
  }
  
  // Define all RooRealVar floating observables.

  /* The trackinfo variable defines the cyclycal nature of the resonances
     and the final state pions, with the following order. 
     See enum ResDaughters.

     index   particle
     1       pi0
     2       pi+
     3       pi-
     --------------
     index   resonance   dtr particles   dtr indices
     1       rho0        pi+ pi-         23
     2       rho-        pi0 pi-         13
     3       rho+        pi0 pi+         12
  */

  int rhoPIndex = -1;
  int rho2sPIndex = -1;
  int rho1700PIndex = -1;
  int nonresIndex = -1;

  // pi+pi0 resonances:
  if (RHOP & componentsBit()) {
    rhoPIndex = nComps();  // store its index
    addComp("Rho+", 1.0, 0.0, 0.7758, 0.1503, PI0_PIP, _defaultRhoSpin);
  }

  if (RHO2SP & componentsBit()) {
    rho2sPIndex = nComps(); // store its index
    addComp("Rho2s+", 0.0, 0.0, 1.465, 0.400, PI0_PIP, _defaultRhoSpin);
  }

  if (RHO1700P & componentsBit()) {
    rho1700PIndex = nComps(); // store its index
    addComp("Rho1700+", 0.0, 0.0, 1.720, 0.250, PI0_PIP, _defaultRhoSpin);
  }

  // pi+pi- resonances:
  if (RHO0 & componentsBit()) {
    addComp("Rho0", 0.59, 10.0, 0.7758, 0.1503, PIP_PIM, _defaultRhoSpin, 
	    rhoPIndex);
  }

  if (RHO2S0 & componentsBit()) {
    addComp("Rho2s0", 0.0, 0.0, 1.465, 0.400, PIP_PIM, _defaultRhoSpin,
	    rho2sPIndex);
  }

  if (RHO17000 & componentsBit()) {
    addComp("Rho17000", 0.0, 0.0, 1.720, 0.250, PIP_PIM, _defaultRhoSpin,
	    rho1700PIndex);
  }

  if (SIGMA & componentsBit()) {
    addComp("Sigma", 0.0, 0.0, 0.500, 0.400, PIP_PIM, SPIN0);
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

  if (F0_1500 & componentsBit()) {
    addComp("F0_1500", 0.0, 0.0, 1.507, 0.109, PIP_PIM, SPIN0);
  }

  if (F0_1710 & componentsBit()) {
    addComp("F0_1710", 0.0, 0.0, 1.714, 0.140, PIP_PIM, SPIN0);
  }

  if (F2 & componentsBit()) {
    addComp("F2", 0.0, 0.0, 1.2754, 0.1851, PIP_PIM, SPIN2);
  }

  if (F2P1525 & componentsBit()) {
    addComp("F2P1525", 0.0, 0.0, 1.525, 0.073, PIP_PIM, SPIN2);
  }

  // pi-pi0 resonances:
  if (RHOM & componentsBit()) {
    addComp("Rho-", 0.65, -4.0, 0.7758, 0.1503, PIM_PI0, _defaultRhoSpin,
	    rhoPIndex);
  }

  if (RHO2SM & componentsBit()) {
    addComp("Rho2s-", 0.0, 0.0, 1.465, 0.400, PIM_PI0, _defaultRhoSpin,
	    rho2sPIndex);
  }

  if (RHO1700M & componentsBit()) {
    addComp("Rho1700-", 0.0, 0.0, 1.720, 0.250, PIM_PI0, _defaultRhoSpin,
	    rho1700PIndex);
  }

  if (NONRES & componentsBit()) {
    nonresIndex = nComps();  // store its index
    addComp("Nonres", 1.03, 77, 0, 0, PIP_PIM_PI0, SPIN0);
  }
 
  if (NRPW_PM & componentsBit()) {
    addComp("NRPW_PM", 0.0, 0.0, 0, 0, PIP_PIM, SPIN1, nonresIndex);
  }
 
  if (NRPW_M0 & componentsBit()) {
    addComp("NRPW_M0", 0.0, 0.0, 0, 0, PIM_PI0, SPIN1, nonresIndex);
  }
 
  if (NRPW_0P & componentsBit()) {
    addComp("NRPW_0P", 0.0, 0.0, 0, 0, PI0_PIP, SPIN1, nonresIndex);
  }
 
  // Masses of the above daughters
  _mDaug[0] = 0.0;
  _mDaug[1] = BdkDalitz::PI0MASS;
  _mDaug[2] = BdkDalitz::PIMASS;
  _mDaug[3] = BdkDalitz::PIMASS;  
}
