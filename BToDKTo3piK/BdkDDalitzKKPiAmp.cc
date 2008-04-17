/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkDDalitzKKPiAmp.cc,v 1.4 2008/01/02 20:30:09 fwinkl Exp $
 * Description:
 *   Implementation of D0 -> K+K-pi0 Dalitz amplitude
 * History:
 *   04 Apr 2007, created, Frank Winklmeier
 *                code by Kalanand Mishra
 *
 * Copyright (C) 2007 Colorado State University and SLAC
 *****************************************************************************/

#include <iostream>

#include "BToDKTo3piK/BdkDDalitzKKPiAmp.hh"
#include "BToDKTo3piK/BdkDalitzBase.hh"
#include "BToDKTo3piK/BdkMath.hh"
#include "BToDKTo3piK/BdkAmp.hh"

using namespace std;


ClassImp(BdkDDalitzKKPiAmp);

// Constructor
BdkDDalitzKKPiAmp::BdkDDalitzKKPiAmp(const char * name, 
                                     const char * title, 
                                     BdkDalitzBase * pdf,
                                     Int_t componentsBit,
                                     Int_t defaultKstarSpin) :
  BdkAbsDDalitzAmp(name, title, pdf, componentsBit),
  _defaultKstarSpin(defaultKstarSpin)
{
  initResonance();
  createParams();
  registerParams(_pdf);
}

// Destructor
BdkDDalitzKKPiAmp::~BdkDDalitzKKPiAmp() 
{
}

Double_t BdkDDalitzKKPiAmp::efficiency(Double_t m12, Double_t m13) const
{
//      if(m12>1.0 && m13>1.0) return 0.0;
//   Double_t m23 =  1.8645*1.8645+ 0.1349766*0.1349766 
//     + 2.0*0.493677*0.493677 - m12 -m13;
//   if(!(m23>1.0 && m23<1.2)) return 0.0;

  return 1.0;

  /*  
  //Here I set efficiency 
  Double_t c0 = -2.39917,
    s1 = 4.30997,
    s2 = -2.35996,
    s3 = 0.431582,
    s4 = 0.421878,
    s5 = -2.07148,
    Normalization = 60.03,
    value = Normalization*(c0 + s1*(m12+m13) + s2*(m12*m12+m13*m13) + 
			   s3*(m13*m13*m13+m12*m12*m12) + 
			   s4*(m12*m12*m13+m12*m13*m13) + s5*m12*m13); 
  
  return value;
  */
  //  return _pdf->efficiency(m12, m13);
}


// initialize the components:
void BdkDDalitzKKPiAmp::initResonance() 
{
  if (nComps()!=0) {
    cout << "BdkDDalitzKKPiAmp::initResonance(): Cannot call initResonance more than once!" << endl;
    return;
  }
    
  // Define all RooRealVar floating observables.

  /* The trackinfo variable defines the cyclycal nature of the resonances
     and the final state pions, with the following order. 
     See enum ResDaughters.

     index   particle
     1       pi0
     2       K+
     3       K-
     --------------
     index   resonance   dtr particles   dtr indices
     1       phi        K+ K-         23
     2       K*-        pi0 K-         13
     3       K*+        pi0 K+         12
  */

  int KstarIndex = -1;
  int Kstar1410PIndex = -1;
  int nonresIndex = -1;

  // K+pi0 resonances:
  if (KSTARP & componentsBit()) {
    KstarIndex = nComps();  // store its index
    addComp("K*+", 1.0, 0.0, 0.89166, 0.0508, PI0_PIP, _defaultKstarSpin);
  }

  if (KSTAR1410P & componentsBit()) {
    Kstar1410PIndex = nComps(); // store its index
    addComp("K*1410+", 1.0, 61.0, 1.414, 0.232, PI0_PIP, _defaultKstarSpin);
  }

  if (NONRESP & componentsBit()) {
    nonresIndex = nComps();  // store its index
    addComp("K+pi0_SW", 1.5, 120.0, 1.414, 0.290, PI0_PIP, SPIN0, -1, KPI_SW);
  }

  // K+K- resonances:
  if (PHI & componentsBit()) {
    addComp("Phi", 0.7, 43.0, 1.0195, 0.00426, PIP_PIM, SPIN1);
  }

  if (F2P1525 & componentsBit()) {
    addComp("F2P1525", 0.6, -40.0, 1.525, 0.073, PIP_PIM, SPIN2);
  }

  if (F0 & componentsBit()) {
    addComp("F0", 0.9, 170, 0.965, 0.044, PIP_PIM, SPIN0, -1, F0);

    // width for the F0 is taken from Flatte' formula and BES parameterization
    // which looked at J/psi-->phi pi+pi- and J/psi-->phi K+K-
  }

  if (A0 & componentsBit()) {
    addComp("A0", 0.9, 0, 0.998, 0.01, PIP_PIM, SPIN0, -1, A0);
  }

  // K-pi0 resonances:
  if (KSTARM & componentsBit()) {
    addComp("K*-", 0.67, -31.0, 0.89166, 0.0508, PIM_PI0, _defaultKstarSpin,
	    KstarIndex);
  }

  if (KSTAR1410M & componentsBit()) {
    addComp("K*1410-", 1.0, 148.0, 1.414, 0.232, PIM_PI0, _defaultKstarSpin,
	    Kstar1410PIndex);
  }

   if (NONRESM & componentsBit()) {
    addComp("K-pi0_SW", 1.5, -45.0, 1.414, 0.290, PIM_PI0, SPIN0,nonresIndex);
  }
 

  // Masses of the above daughters
  _mDaug[0] = 0.0;
  _mDaug[1] = BdkDalitz::PI0MASS;
  _mDaug[2] = BdkDalitz::KMASS;
  _mDaug[3] = BdkDalitz::KMASS;


  _r = new RooRealVar(TString(GetName()) + ".r",
                              TString(GetTitle()) + " Effective Range",
                              3.32);
  _params.addOwned(*_r);

  _a = new RooRealVar(TString(GetName()) + ".a",
                              TString(GetTitle()) + " Scattering Length",
                              2.07);
  _params.addOwned(*_a);

  _R = new RooRealVar(TString(GetName()) + ".R",
                              TString(GetTitle()) + " Resonant Amplitude",
                              1.0);
  _params.addOwned(*_R);

  _B = new RooRealVar(TString(GetName()) + ".B",
                              TString(GetTitle()) + " Background Amplitude",
                              1.0);
  _params.addOwned(*_B);

  _phiR = new RooRealVar(TString(GetName()) + ".phiR",
                              TString(GetTitle()) + " Resonant Phase",
                              0.0);
  _params.addOwned(*_phiR);

  _phiB = new RooRealVar(TString(GetName()) + ".phiB",
                              TString(GetTitle()) + " Background Phase",
                              0.0);
  _params.addOwned(*_phiB);

  _E791 = new RooRealVar(TString(GetName()) + ".E791",
                              TString(GetTitle()) + " whether using E791 parameters",
                              0.0);
  _params.addOwned(*_E791);
}



RooComplex BdkDDalitzKKPiAmp::matrixElement(Double_t m12, Double_t m13, Int_t i) const
{  
  if (_typeRes[i]==A0RES) {
    Double_t s = kinematics(m12,m13,_trackinfo[i]);
    return BdkAmp::a0Amp(s, _massRes[i]->getVal(), 0.105, 0.102); // this is a0(980)
  }
  else if (_typeRes[i]==F0RES) {
    Double_t s = kinematics(m12,m13,_trackinfo[i]);
    return BdkAmp::f0Amp(s, _massRes[i]->getVal(), 0.165, 0.165*4.21); // this is f0(980)
  }
  else if (_typeRes[i]==KPI_SW) {
    Double_t s = kinematics(m12,m13,_trackinfo[i]);
    /*********  LASS parameters:  return LASS(s, i);      *************/ 
    /*********  E791 parameters:  return E791Table(s);    *************/ 
    
    if(_E791->getVal()>0.0) return E791Table(s);
    else return LASS(s, i);
  }
  // Standard Breit-Wigner
  else return BdkAbsDDalitzAmp::matrixElement(m12, m13, i);
}




//replace Breit-Wigner with LASS

// RooComplex BdkDDalitzKKPiAmp::LASS(Double_t s, Int_t i) const 
// {

//   Double_t _mass = _massRes[i]->getVal();  //K*(1430) mass
//   Double_t _width = _gammaRes[i]->getVal(); //K*(1430) width
//   Double_t _m1=0.493677;       //This is K
//   Double_t _m2=0.1349766;      //This is pi0 
    
//   Double_t q      = pionCMmom(s, _m1, _m2);  //pion momentum in C.M. frame
//   Double_t q0     = pionCMmom(_mass*_mass, _m1, _m2);  
//   Double_t GammaM = _width * _mass/sqrt(s) * q/q0;  //mass dependent width


//   //calculate the background phase motion
//   Double_t cot_deltaB = 1.0/(_a->getVal()*q) + 0.5*_r->getVal()*q;
//   Double_t _deltaB = atan( 1.0/cot_deltaB);
//   Double_t totalB = (_deltaB + _phiB->getVal()*DEGTORAD) ;
  
//   //calculate the resonant phase motion
//   Double_t deltaR = atan((_mass*GammaM/(_mass*_mass - s)));
//   Double_t totalR = deltaR + _phiR->getVal()*DEGTORAD;
  
//   //sum them up
//   RooComplex  bkgB,resT;
//   bkgB = RooComplex(_B->getVal()*sin(totalB),0)*RooComplex(cos(totalB),sin(totalB));
//   resT = RooComplex(_R->getVal()*sin(deltaR),0)*RooComplex(cos(totalR),sin(totalR))*
//     RooComplex(cos(2*totalB),sin(2*totalB));
//   RooComplex T = bkgB + resT;  
//   T = T*RooComplex(sqrt(s)/q,0);  
//   return T;
// }


RooComplex BdkDDalitzKKPiAmp::LASS(Double_t s, Int_t i) const
{
  Double_t _resMass = _massRes[i]->getVal();   //K*(1430) mass
  Double_t _resWidth = _gammaRes[i]->getVal(); //K*(1430) width
  Int_t n = _trackinfo[i];
  Double_t _mDaugSum  = _mDaug[n] + _mDaug[1];
  Double_t _mDaugSumSq = _mDaugSum*_mDaugSum;

  Double_t _mDaugDiff = _mDaug[n] - _mDaug[1];
  Double_t _mDaugDiffSq = _mDaugDiff*_mDaugDiff;

  RooComplex resAmplitude = RooComplex(0.0, 0.0);
  RooComplex bkgAmplitude = RooComplex(0.0, 0.0);
  RooComplex totAmplitude = RooComplex(0.0, 0.0);
  
  if (s < 1e-10) {
    cout<<"Warning in LASS ::amplitude. Mass < 1e-10."<<endl;
    return RooComplex(0.0, 0.0);
  } else if (s > 2.2) {
    return RooComplex(0.0, 0.0);
  }

  //---------------------------
  // First do the resonant part
  //---------------------------

  // Calculate the width of the resonance (as a function of mass)
  // q is the momentum of either daughter in the resonance rest-frame
  Double_t q(0.0);
  if ((s - _mDaugSumSq)>0.0) { // protect against negative sqrt due to rounding errors
    q = sqrt((s - _mDaugSumSq)*(s - _mDaugDiffSq))/(2.0*sqrt(s));
  }
  Double_t _q0 = sqrt((_resMass*_resMass - _mDaugSumSq)*
	   (_resMass*_resMass - _mDaugDiffSq))/(2.0*_resMass);
  Double_t qRatio = q/_q0;

  Double_t totWidth = _resWidth*qRatio*(_resMass/sqrt(s));

  Double_t massSqTerm = _resMass*_resMass - s;

  // Compute the complex amplitude
  resAmplitude = RooComplex(massSqTerm, _resMass*totWidth);

  // Scale by the denominator factor
  resAmplitude = resAmplitude*
    ((_resMass*_resMass*_resWidth/_q0)/(massSqTerm*massSqTerm +
				 _resMass*_resMass*totWidth*totWidth));

  // Calculate the phase shift term
  Double_t deltaB = TMath::ATan((2.0*_a->getVal()*q)/
				(2.0 + _a->getVal()*_r->getVal()*q*q));
  Double_t cos2PhaseShift = TMath::Cos(2.0*(deltaB + _phiB->getVal()*DEGTORAD));
  Double_t sin2PhaseShift = TMath::Sin(2.0*(deltaB + _phiB->getVal()*DEGTORAD));
  RooComplex phaseShift = RooComplex(cos2PhaseShift, sin2PhaseShift);

  // Add in the R e^{i phiR} term
  Double_t reR = _R->getVal() * TMath::Cos(_phiR->getVal()*DEGTORAD);
  Double_t imR = _R->getVal() * TMath::Sin(_phiR->getVal()*DEGTORAD);
  RooComplex R = RooComplex(reR, imR);

  // Multiply by the phase shift and R e^{i phiR}
  resAmplitude = resAmplitude * phaseShift * R;


  //--------------------------------
  // Now do the effective range part
  //--------------------------------

  // Form the real and imaginary parts
  Double_t realTerm = q/TMath::Tan(deltaB + _phiB->getVal()*DEGTORAD);
  Double_t imagTerm = q;

  // Compute the complex amplitude
  bkgAmplitude = RooComplex(realTerm, imagTerm);
  bkgAmplitude = bkgAmplitude*(sqrt(s)*_B->getVal());

  // Scale by the denominator factor
  bkgAmplitude = bkgAmplitude*(1.0/(realTerm*realTerm 
				    + imagTerm*imagTerm));


  //------------------
  // Add them together
  //------------------

  totAmplitude = bkgAmplitude + resAmplitude;

  return totAmplitude;
}



//replace Breit-Wigner with Brian's Table III (E791 MIPWA D+ --> Kpipi paper

RooComplex BdkDDalitzKKPiAmp::E791Table(Double_t s) const 
{
  if( sqrt(s)<0.672 || sqrt(s)>1.707 ) return RooComplex(0.0,0.0);

  double mass[38], FD[38], Amp[38], Phas[38];
  mass[0] = 0.672;        FD[0] = 0.26;       Amp[0] = 8.37;      Phas[0] = -102;
  mass[1] = 0.719;        FD[1] = 0.27;       Amp[1] = 9.04;      Phas[1] = -96;
  mass[2] = 0.764;        FD[2] = 0.29;       Amp[2] = 7.82;      Phas[2] = -73;
  mass[3] = 0.807;        FD[3] = 0.31;       Amp[3] = 7.42;      Phas[3] = -77;
  mass[4] = 0.847;        FD[4] = 0.33;       Amp[4] = 6.47;      Phas[4] = -60;
  mass[5] = 0.885;        FD[5] = 0.34;       Amp[5] = 5.57;      Phas[5] = -54;
  mass[6] = 0.922;        FD[6] = 0.36;       Amp[6] = 5.90;      Phas[6] = -68;
  mass[7] = 0.958;        FD[7] = 0.38;       Amp[7] = 6.17;      Phas[7] = -72;
  mass[8] = 0.992;        FD[8] = 0.40;       Amp[8] = 4.87;      Phas[8] = -41;
  mass[9] = 1.025;        FD[9] = 0.42;       Amp[9] = 4.42;      Phas[9] = -43;
  mass[10] = 1.057;       FD[10] = 0.44;      Amp[10] = 4.02;     Phas[10] = -38;
  mass[11] = 1.088;       FD[11] = 0.46;      Amp[11] = 3.74;     Phas[11] = -22;
  mass[12] = 1.118;       FD[12] = 0.49;      Amp[12] = 3.81;     Phas[12] = -29;
  mass[13] = 1.147;       FD[13] = 0.51;      Amp[13] = 3.16;     Phas[13] = -3;
  mass[14] = 1.176;       FD[14] = 0.53;      Amp[14] = 3.21;     Phas[14] = -11;
  mass[15] = 1.204;       FD[15] = 0.55;      Amp[15] = 2.86;     Phas[15] = -3;
  mass[16] = 1.231;       FD[16] = 0.58;      Amp[16] = 3.11;     Phas[16] = -3;
  mass[17] = 1.258;       FD[17] = 0.60;      Amp[17] = 2.92;     Phas[17] = 8;
  mass[18] = 1.284;       FD[18] = 0.62;      Amp[18] = 2.80;     Phas[18] = 11;
  mass[19] = 1.310;       FD[19] = 0.65;      Amp[19] = 2.77;     Phas[19] = 11;
  mass[20] = 1.335;       FD[20] = 0.67;      Amp[20] = 2.83;     Phas[20] = 22;
  mass[21] = 1.360;       FD[21] = 0.69;      Amp[21] = 2.73;     Phas[21] = 31;
  mass[22] = 1.384;       FD[22] = 0.71;      Amp[22] = 2.29;     Phas[22] = 30;
  mass[23] = 1.408;       FD[23] = 0.74;      Amp[23] = 2.38;     Phas[23] = 46;
  mass[24] = 1.431;       FD[24] = 0.76;      Amp[24] = 2.05;     Phas[24] = 55;
  mass[25] = 1.454;       FD[25] = 0.78;      Amp[25] = 1.59;     Phas[25] = 64;
  mass[26] = 1.477;       FD[26] = 0.80;      Amp[26] = 1.33;     Phas[26] = 80;
  mass[27] = 1.499;       FD[27] = 0.82;      Amp[27] = 1.23;     Phas[27] = 74;
  mass[28] = 1.522;       FD[28] = 0.84;      Amp[28] = 0.66;     Phas[28] = 34;
  mass[29] = 1.543;       FD[29] = 0.86;      Amp[29] = 0.57;     Phas[29] = 18;
  mass[30] = 1.565;       FD[30] = 0.88;      Amp[30] = 0.50;     Phas[30] = 22;
  mass[31] = 1.586;       FD[31] = 0.90;      Amp[31] = 1.18;     Phas[31] = 10;
  mass[32] = 1.607;       FD[32] = 0.92;      Amp[32] = 1.35;     Phas[32] = 11;
  mass[33] = 1.627;       FD[33] = 0.93;      Amp[33] = 1.11;     Phas[33] = 19;
  mass[34] = 1.648;       FD[34] = 0.95;      Amp[34] = 1.37;     Phas[34] = 2;
  mass[35] = 1.668;       FD[35] = 0.96;      Amp[35] = 1.82;     Phas[35] = 28;
  mass[36] = 1.687;       FD[36] = 0.98;      Amp[36] = 1.16;     Phas[36] = 8;
  mass[37] = 1.707;       FD[37] = 0.99;      Amp[37] = 1.47;     Phas[37] = 11;


  int n = 0;
  for(int i=0; i<36; i++) { if( sqrt(s)>mass[i] && sqrt(s)<mass[i+1] ) n = i; }
  

  double slope = (sqrt(s)-mass[n])/(mass[n+1]-mass[n]);
  double ReLo = FD[n]*Amp[n]*cos(Phas[n]*DEGTORAD);
  double ImLo = FD[n]*Amp[n]*sin(Phas[n]*DEGTORAD);
  double ReDiff = FD[n+1]*Amp[n+1]*cos(Phas[n+1]*DEGTORAD) - ReLo;
  double ImDiff = FD[n+1]*Amp[n+1]*sin(Phas[n+1]*DEGTORAD) - ImLo;
  RooComplex loV = RooComplex(ReLo, ImLo);
  RooComplex diff = RooComplex(ReDiff, ImDiff)*slope;


  RooComplex val = loV + diff;

  return val;
}














