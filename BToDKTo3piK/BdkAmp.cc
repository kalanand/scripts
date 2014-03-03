/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkAmp.cc,v 1.2 2008/01/03 18:11:16 fwinkl Exp $
 * Description:
 *   Collection of amplitudes
 * History:
 *   09 Nov 2007, created, Frank Winklmeier
 *
 * Copyright (C) 2007 Colorado State University and SLAC
 *****************************************************************************/

#include "BToDKTo3piK/BdkAmp.hh"

#include "BToDKTo3piK/BdkMath.hh"
#include "BToDKTo3piK/BdkDalitz.hh"

ClassImp(BdkAmp);
  
// Coupled channel Breit-Wigner function for a0(980)
// Code from Kalanad Mishra
// gXX are the squared coupling constants!
RooComplex BdkAmp::a0Amp(Double_t s, Double_t m0,
                         Double_t gEtaPi, Double_t gKK)
{

  //  Double_t m0 = 0.998;
  //  Double_t gEtaPi = 0.105;        
  //  Double_t gKK   = 0.102;   

  Double_t mSumSq1_ = 4.0*BdkDalitz::ETAMASS*BdkDalitz::PI0MASS;
  Double_t mSumSq2_ = 4.0*sqr(BdkDalitz::KMASS);
  Double_t mSumSq3_ = 4.0*sqr(BdkDalitz::K0MASS);

  Double_t dMSq = m0*m0 - s;

  Double_t rho1(0.0), rho2(0.0);

  ////////////////////////////////////////
  gKK = sqr(gKK);
  gEtaPi = sqr(gEtaPi);
  ////////////////////////////////////////
  
  if (s > mSumSq1_) {
    rho1 += TMath::Sqrt(1.0 - mSumSq1_/s);
    if (s > mSumSq2_) {
      rho2 = 0.5*TMath::Sqrt(1.0 - mSumSq2_/s);
      if (s > mSumSq3_) {
        rho2 += 0.5*TMath::Sqrt(1.0 - mSumSq3_/s);
      } else {
        // Continue analytically below higher channel thresholds
        // This contributes to the real part of the amplitude denominator
        dMSq += gKK*m0*0.5*TMath::Sqrt(mSumSq3_/s - 1.0);
      }
    } else {
      // Continue analytically below higher channel thresholds
      // This contributes to the real part of the amplitude denominator
      rho2 = 0.0;
      dMSq += gKK*m0*(0.5*TMath::Sqrt(mSumSq2_/s - 1.0) + 0.5*TMath::Sqrt(mSumSq3_/s - 1.0));
    }
  } else {
    // Continue analytically below higher channel thresholds
    // This contributes to the real part of the amplitude denominator
    dMSq += gEtaPi*m0*2.0*TMath::Sqrt(mSumSq1_/s - 1.0);
  }
  
  //  Double_t widthTerm = gEtaPi*rho1*m0 + gKK*rho2*m0;
  Double_t widthTerm = gEtaPi*rho1 + gKK*rho2;
  RooComplex resAmplitude = RooComplex(dMSq, widthTerm);  
  Double_t denomFactor = dMSq*dMSq + widthTerm*widthTerm;  
  Double_t invDenomFactor = 0.0;
  if (denomFactor > 1e-10) {invDenomFactor = 1.0/denomFactor;}  
  resAmplitude = resAmplitude*(invDenomFactor);

  return resAmplitude;
}


// Coupled channel Breit-Wigner function for f0(980)
// Code from Kalanad Mishra
// gXX are the squared coupling constants!
RooComplex BdkAmp::f0Amp(Double_t s, Double_t m0,
                         Double_t gPiPi, Double_t gKK)
{
  //  Double_t gPiPi = 0.165;        // +/- 0.010 +/- 0.015 GeV/c^2
  //  Double_t gKK   = gPiPi*4.21;   // +/- 0.25 +/- 0.21

  Double_t mSumSq0_ = 4.0*sqr(BdkDalitz::PI0MASS);
  Double_t mSumSq1_ = 4.0*sqr(BdkDalitz::PIMASS);
  Double_t mSumSq2_ = 4.0*sqr(BdkDalitz::KMASS);
  Double_t mSumSq3_ = 4.0*sqr(BdkDalitz::K0MASS);

  Double_t dMSq = m0*m0 - s;

  ////////////////////////////////////////
  gPiPi = sqr(gPiPi);
  gKK = sqr(gKK);
  ////////////////////////////////////////
  
  Double_t rho1(0.0), rho2(0.0);
  if (s > mSumSq0_) {
    rho1 = TMath::Sqrt(1.0 - mSumSq0_/s)/3.0;
    if (s > mSumSq1_) {
      rho1 += 2.0*TMath::Sqrt(1.0 - mSumSq1_/s)/3.0;
      if (s > mSumSq2_) {
	rho2 = 0.5*TMath::Sqrt(1.0 - mSumSq2_/s);
	if (s > mSumSq3_) {
	  rho2 += 0.5*TMath::Sqrt(1.0 - mSumSq3_/s);
	} else {
	  // Continue analytically below higher channel thresholds
	  // This contributes to the real part of the amplitude denominator
	  dMSq += gKK*m0*0.5*TMath::Sqrt(mSumSq3_/s - 1.0);
	}
      } else {
	// Continue analytically below higher channel thresholds
	// This contributes to the real part of the amplitude denominator
	rho2 = 0.0;
	dMSq += gKK*m0*(0.5*TMath::Sqrt(mSumSq2_/s - 1.0) + 0.5*TMath::Sqrt(mSumSq3_/s - 1.0));
      }
    } else {
      // Continue analytically below higher channel thresholds
      // This contributes to the real part of the amplitude denominator
      dMSq += gPiPi*m0*2.0*TMath::Sqrt(mSumSq1_/s - 1.0)/3.0;
    }
  }
  // Double_t widthTerm = gPiPi*rho1*m0 + gKK*rho2*m0;
  Double_t widthTerm = gPiPi*rho1 + gKK*rho2;  
  RooComplex resAmplitude = RooComplex(dMSq, widthTerm);  
  Double_t denomFactor = dMSq*dMSq + widthTerm*widthTerm;  
  Double_t invDenomFactor = 0.0;
  if (denomFactor > 1e-10) {invDenomFactor = 1.0/denomFactor;}  
  resAmplitude = resAmplitude*(invDenomFactor);

  return resAmplitude;
}


// Coupled channel Breit-Wigner function for a0(980)
// According to Antimo's BAD 120
RooComplex BdkAmp::a0AmpAntimo(Double_t s, Double_t m0,
                               Double_t gEtaPi, Double_t gKK)

{
  Double_t rhoEtaPi = sqrt((m0*m0 - sqr(BdkDalitz::ETAMASS+BdkDalitz::PI0MASS)) * 
			   (m0*m0 - sqr(BdkDalitz::ETAMASS-BdkDalitz::PI0MASS)))/(m0*m0);

  Double_t rhoKK = sqrt(m0*m0 - sqr(2*BdkDalitz::K0MASS))/m0;
   
  Double_t widthTerm = rhoEtaPi*sqr(gEtaPi) + rhoKK*sqr(gKK);

  RooComplex amp(gKK);
  amp = amp / RooComplex(m0*m0-s,-widthTerm);
  return amp;
}
 
// Gounaris-Sakurai (GS) parameterization of the rho
// PRL 21,244 (1968) http://prola.aps.org/abstract/PRL/v21/i4/p244_1
// Equations also in BAD 637 Version 4
// Code taken from RooDKDalitz/RooPto3PAmp.cc
// Parameters:
//     m0      Resonance (rho) mass
//     g0      Nominal width of resonance at m0
//     k0      Pion momentum in resonance rest frame at energy m0
//     s       CM energy squared
//     g       Energy dependent (running) width at energy s
//     k       Pion momentum in resonance rest frame at energy s
RooComplex BdkAmp::GS_Rho(Double_t m0, Double_t g0, Double_t k0,
			  Double_t s, Double_t g, Double_t k)
{
  Double_t mPi = BdkDalitz::PIMASS;
  Double_t mPiSqr = mPi*mPi;
  Double_t m = TMath::Sqrt(s);

  Double_t GS_f = g0*m0*m0/(k0*k0*k0)*(k*k*(GS_h(m,k)-GS_h(m0,k0))
				       + (m0*m0-s)*k0*k0*GS_dhds(m0,k0));

  Double_t GS_d = 3/M_PI*mPiSqr/(k0*k0)*log((m0+2*k0)/(2*mPi)) + 
    m0/(2*M_PI*k0) - mPiSqr*m0/(M_PI*k0*k0*k0);

  RooComplex amp(1+GS_d*g0/m0);
  amp = amp / RooComplex(m0*m0-s+GS_f,-m*g);
  return amp;
}

inline Double_t BdkAmp::GS_h(const Double_t& m, const Double_t& k) 
{
  return 2/M_PI * k/m*log((m+2*k)/(2*BdkDalitz::PIMASS));
}

// First derivative of GS_h
inline Double_t BdkAmp::GS_dhds(const Double_t& m0, const Double_t& k0) 
{
  return GS_h(m0,k0)*( 1/(8*k0*k0) - 1/(2*m0*m0) ) + 1/(2*M_PI*m0*m0);
}

