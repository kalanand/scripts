/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkDalitzHist.cc,v 1.4 2006/05/22 21:53:04 fwinkl Exp $
 * Description:
 *   Histogram PDF that knows about the Dalitz boundaries
 * History:
 *   Mar 03, 2006     Frank Winklmeier, created
 *
 * Copyright (C) 2006 Colorado State University and SLAC
 *****************************************************************************/

#include "BToDKTo3piK/BdkDalitzHist.hh"

#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooRealProxy.hh"

#include "TH2.h"

/// constructor
/// m12 = x-axis, m13 = y-axis of hist
BdkDalitzHist::BdkDalitzHist(const char* name, const char* title,
                             BdkDalitzBase::Flavor flavor, BdkDalitzBase::Mode DdecMode,
                             RooAbsReal &m12, RooAbsReal &m13, const TH2& hist) :
  BdkDalitzBase(name,title,flavor,DdecMode),
  _m12("m12","Invariant Mass square of pi0 pi+",this,m12),
  _m13("m13","Invariant Mass square of pi0 pi-",this,m13),
  _hist(hist),
  _area(0)
{
}

/// copy constructor
BdkDalitzHist::BdkDalitzHist(const BdkDalitzHist& other, const char* name) :
  BdkDalitzBase(other, name),
  _m12("m12",this,other._m12),
  _m13("m13",this,other._m13),
  _hist(other._hist)
{
  if (other._area) _area = (TH2*)other._area->Clone(name+TString(".area"));
  else _area = 0;
}

/// destructor
BdkDalitzHist::~BdkDalitzHist()
{
  delete _area;
}

/// returns 0 for events outside the Dalitz plot
Double_t BdkDalitzHist::evaluate() const
{
  if (inDalitz(_m12,_m13))
    return _hist.GetBinContent(((TH2&)_hist).FindBin(_m12,_m13));
  else
    return 0;
}


Int_t BdkDalitzHist::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                                           const char* rangeName) const
{
  // Only analytical integrals over the full range are defined
  if (rangeName!=0) return 0;

  if (matchArgs(allVars,analVars,_m12,_m13)) return 1;
  if (matchArgs(allVars,analVars,_m12)) return 2;
  if (matchArgs(allVars,analVars,_m13)) return 3;

  return 0;
}


Double_t BdkDalitzHist::analyticalIntegral(Int_t code, const char* /*rangeName*/) const 
{
  Double_t r = 0;

  // Simplest scenario, integration over all dependents
  if (code==1) {
    if (!_area) {
      // Calculate area of each bin and store in separate histogram
      _area = (TH2*)_hist.Clone(GetName()+TString(".area"));
      weightBins(_area, kFALSE);
    }
    for (Int_t xbin=1; xbin<=_hist.GetNbinsX(); xbin++) {
      for (Int_t ybin=1; ybin<=_hist.GetNbinsY(); ybin++) {        
        r += _hist.GetBinContent(xbin,ybin) * _area->GetBinContent(xbin,ybin);
      }
    }
  }
  else if (code==2) {
      // Dalitz boundaries
    Double_t max = m12Max(_m13);
    Double_t min = m12Min(_m13);
    if (max==0 || min==0) return 0;  // not in Dalitz region

    Double_t m12VetoMin = -1;
    Double_t m12VetoMax = -1;
    if (hasM23Veto()) {
      // Veto window boundaries in m13
      m12VetoMin = mOtherVal(_m13, m23VetoMax());
      m12VetoMax = mOtherVal(_m13, m23VetoMin());
      // Is veto window part of Dalitz region?
      if (m12VetoMin>max || m12VetoMax<min) m12VetoMin = m12VetoMax = -1;
      else {
	// Only consider veto window inside Dalitz region
	m12VetoMin = TMath::Max(m12VetoMin,min);
	m12VetoMax = TMath::Min(m12VetoMax,max);
      }
    }

    // Loop over m12 bins
    Int_t ybin = _hist.GetYaxis()->FindBin(_m13);
    for (Int_t xbin=1; xbin<=_hist.GetNbinsX(); xbin++) {

      // Calculate length of bin inside Dality region
      Double_t x1 = _hist.GetXaxis()->GetBinLowEdge(xbin);
      Double_t x2 = _hist.GetXaxis()->GetBinUpEdge(xbin);
      Double_t len = TMath::Min(x2,max)-TMath::Max(x1,min);
      if (len<0) len = 0;

      // Subtract veto window if neccessary
      if (len>0 && m12VetoMin>0 && m12VetoMax>0) {
	Double_t veto = TMath::Min(x2,m12VetoMax)-TMath::Max(x1,m12VetoMin);
	if (veto>0) len -= veto;
      }

      Double_t w = _hist.GetBinContent(xbin,ybin);
      r += w * len;
    }
  }
  else if (code==3) {
     // Dalitz boundaries
    Double_t max = m13Max(_m12);
    Double_t min = m13Min(_m12);
    //cout << "Dalitz boundaries at m12 = "<<_m12 << " are "<<min<<" to "<<max<<endl;
    if (max==0 || min==0) return 0;  // not in Dalitz region

    Double_t m13VetoMin = -1;
    Double_t m13VetoMax = -1;
    if (hasM23Veto()) {
      // Veto window boundaries in m13
      m13VetoMin = mOtherVal(_m12, m23VetoMax());
      m13VetoMax = mOtherVal(_m12, m23VetoMin());
      // Is veto window part of Dalitz region?
      if (m13VetoMin>max || m13VetoMax<min) m13VetoMin = m13VetoMax = -1;
      else {
	// Only consider veto window inside Dalitz region
	m13VetoMin = TMath::Max(m13VetoMin,min);
	m13VetoMax = TMath::Min(m13VetoMax,max);
      }
      //cout << "Veto window is "<<m13VetoMin<<" to "<<m13VetoMax<<endl;
    }

    // Loop over m13 bins
    Int_t xbin = _hist.GetXaxis()->FindBin(_m12);
    for (Int_t ybin=1; ybin<=_hist.GetNbinsY(); ybin++) {

      // Calculate length of bin inside Dality region
      Double_t y1 = _hist.GetYaxis()->GetBinLowEdge(ybin);
      Double_t y2 = _hist.GetYaxis()->GetBinUpEdge(ybin);
      Double_t len = TMath::Min(y2,max)-TMath::Max(y1,min);
      if (len<0) len = 0;

      // Subtract veto window if neccessary
      if (len>0 && m13VetoMin>0 && m13VetoMax>0) {
	Double_t veto = TMath::Min(y2,m13VetoMax)-TMath::Max(y1,m13VetoMin);
	if (veto>0) len -= veto;
      }

      Double_t w = _hist.GetBinContent(xbin,ybin);
      //cout << "xbin = "<<xbin<<"  ybin = "<<ybin<<"("<<y1<<","<<y2<<")"<<"  w = "<<w<<"  len = "<<len<<endl;
      r += w * len;
    }
  }
  else assert(0);
  
  //cout <<" = "<<r<<endl;
  return r;
}
