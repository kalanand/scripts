/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkDalitz.cc,v 1.8 2007/04/18 11:57:05 fwinkl Exp $
 * Description:
 *   Base class for Dalitz PDFs
 * History:
 *   26 Oct 2005, created, Frank Winklmeier
 *
 * Copyright (C) 2005 Colorado State University and SLAC
 *****************************************************************************/

#include <math.h>
#include "BToDKTo3piK/BdkMath.hh"
#include "BToDKTo3piK/BdkDalitz.hh"
#include "BToDKTo3piK/BdkDalitzCfg.hh"

#include "RooFitCore/RooAbsFunc.hh"
#include "RooFitCore/RooArgSet.hh"
#include "TH2.h"

ClassImp(BdkDalitz)

// Static members:
RooRealVar* BdkDalitz::_m23VetoMass = 0;
const Double_t BdkDalitz:: D0MASS  = 1.8646;
const Double_t BdkDalitz:: PI0MASS = 0.1349766;
const Double_t BdkDalitz:: PIMASS  = 0.13957;
const Double_t BdkDalitz:: KMASS   = 0.493677;
const Double_t BdkDalitz:: K0MASS  = 0.497648;
const Double_t BdkDalitz:: ETAMASS  = 0.54751;

// constructors:
BdkDalitz::BdkDalitz(BdkDalitz::Flavor flavor, 
                     BdkDalitz::Mode DdecMode)
{
  _flavor = flavor;
  setDdecMode(DdecMode);
}

BdkDalitz::BdkDalitz(const BdkDalitz& other) :
  _flavor(other._flavor),
  _M(other._M),
  _m1(other._m1),
  _m2(other._m2),
  _m3(other._m3),
  _DdecMode(other._DdecMode)
{
}

// destructor:  
BdkDalitz::~BdkDalitz() 
{
}


/// set the D decay mode
/// BE CAREFULL: The daughter masses must be in the same order as
/// the _mDaug array in the dalitz amp. class!
void BdkDalitz::setDdecMode(BdkDalitz::Mode DdecMode)
{
  _DdecMode = DdecMode;
  
  // D0 -> pi+ pi- pi0
  if (DdecMode==BdkDalitz::PPP0)       setMasses(D0MASS,PI0MASS,PIMASS,PIMASS);
  // D0 -> K+ pi- pi0
  else if (DdecMode==BdkDalitz::KPP0)  setMasses(D0MASS,PI0MASS,KMASS,PIMASS);
  // D0 -> K+ K- pi0
  else if (DdecMode==BdkDalitz::KKP0)  setMasses(D0MASS,PI0MASS,KMASS,KMASS);
  // D0 -> pi+ pi- Ks
  else if (DdecMode==BdkDalitz::PPK0)  setMasses(D0MASS,K0MASS,PIMASS,PIMASS);
  // D0 -> K+ K- Ks
  else if (DdecMode==BdkDalitz::KKK0)  setMasses(D0MASS,K0MASS,KMASS,KMASS);
  // D0 -> K- pi+ Ks
  else if (DdecMode==BdkDalitz::KPK0)  setMasses(D0MASS,K0MASS,PIMASS,KMASS);
  
  else                                 setMasses(-1,-1,-1,-1);
}


/// Returns maximum of m13 at given m12
Double_t BdkDalitz::m13Max(Double_t m12, Double_t M, 
                               Double_t m1, Double_t m2, Double_t m3)
{
  if ((m12<sqr(m1+m2)) || (m12>sqr(M-m3))) return 0;

  // Energies of particles 1 and 3 in m12 restframe
  Double_t e1star = (m12 - m2*m2 + m1*m1) / (2.0*sqrt(m12));
  Double_t e3star = (M*M - m12 - m3*m3  ) / (2.0*sqrt(m12));

  return sqr(e1star+e3star)-sqr(sqrt(e1star*e1star-m1*m1)-sqrt(e3star*e3star - m3*m3));
}

/// Returns minimum of m13 at given m12
Double_t BdkDalitz::m13Min(Double_t m12, Double_t M, 
                               Double_t m1, Double_t m2, Double_t m3)
{
  if ((m12<sqr(m1+m2)) || (m12>sqr(M-m3))) return 0;

  // Energies of particles 1 and 3 in m12 restframe
  Double_t e1star = (m12 - m2*m2 + m1*m1) / (2.0*sqrt(m12));
  Double_t e3star = (M*M - m12 - m3*m3  ) / (2.0*sqrt(m12));

  return sqr(e1star+e3star)-sqr(sqrt(e1star*e1star-m1*m1)+sqrt(e3star*e3star - m3*m3));
}

/// Returns maximum of m12 at given m13
Double_t BdkDalitz::m12Max(Double_t m13, Double_t M, 
                               Double_t m1, Double_t m2, Double_t m3)
{
  if ((m13<sqr(m1+m3)) || (m13>sqr(M-m2))) return 0;

  // Energies of particles 1 and 2 in m13 restframe
  Double_t e1prime = (m13 - m3*m3 + m1*m1) / (2.0*sqrt(m13));
  Double_t e2prime = (M*M - m13 - m2*m2) / (2.0*sqrt(m13));

  return sqr(e1prime+e2prime)-sqr(sqrt(e1prime*e1prime-m1*m1)-sqrt(e2prime*e2prime - m2*m2));
}

/// Returns minumum of m12 at given m13
Double_t BdkDalitz::m12Min(Double_t m13, Double_t M, 
                               Double_t m1, Double_t m2, Double_t m3)
{
  if ((m13<sqr(m1+m3)) || (m13>sqr(M-m2))) return 0;

  // Energies of particles 1 and 2 in m13 restframe
  Double_t e1prime = (m13 - m3*m3 + m1*m1) / (2.0*sqrt(m13));
  Double_t e2prime = (M*M - m13 - m2*m2) / (2.0*sqrt(m13));

  return sqr(e1prime+e2prime)-sqr(sqrt(e1prime*e1prime-m1*m1)+sqrt(e2prime*e2prime - m2*m2));
}


/// Decide if the Dalitz point (m12,m13) is within the kinematic region
Bool_t BdkDalitz::inDalitzBounds(Double_t m12, Double_t m13, 
                                 Double_t M, Double_t m1, Double_t m2, Double_t m3)
{ 
  // kinematic limits
  Double_t m13High = m13Max(m12,M,m1,m2,m3);
  if (m13High<=0) return kFALSE;

  Double_t m13Low = m13Min(m12,M,m1,m2,m3);
  if (m13Low<=0) return kFALSE; 

  // decide
  if ( (m13 < m13Low) || (m13 > m13High)) return kFALSE;

  return kTRUE;
}

/// Decide if the Dalitz point (m12,m13) is within the Dalitz 
/// boundaries of the decay M -> m1 m2 m3
/// Also check the veto window.
Bool_t BdkDalitz::inDalitz(Double_t m12, Double_t m13, 
                           Double_t M, Double_t m1, Double_t m2, Double_t m3)
{
  if (inDalitzBounds(m12,m13,M,m1,m2,m3) &&
      !inM23Veto(m12,m13,M,m1,m2,m3)) return kTRUE;
  else return kFALSE;
}


/// check if Dalitz point (m12,m13) is inside the m23 mass veto window
Bool_t BdkDalitz::inM23Veto(Double_t m12, Double_t m13, 
				Double_t M, Double_t m1, Double_t m2, Double_t m3)
{
  if (!hasM23Veto()) return kFALSE;

  Double_t m23 = mOtherVal(m12, m13, M, m1, m2, m3);
  if (m23>m23VetoMin() && m23<m23VetoMax()) return kTRUE;
  else return kFALSE;
}


/// Find intersection of func with Dalitz boundary
/// Result in m12Left/m13Right. -1 if no intersection found.
/// Returns number of intersections found.
Int_t BdkDalitz::findIntersection(const RooAbsFunc& func, 
                                  Double_t& m12Left, Double_t& m12Right) const
{
  const Double_t dx = 0.001;   // step in m12
  const Double_t dy = 0.01;    // accuracy in m13

  Double_t max = m12MaxLimit();
  m12Left = -1;     // Default value indicating error
  m12Right = -1;
  
  // Loop for upper and lower Dality boundary
  for (int j=0; j<2; j++) {
    Double_t m12 = m12MinLimit();    // start on left side of Dalitz plot  
    while (m12<=max) {
      Double_t y = (j==0 ? m13Max(m12) : m13Min(m12));   // upper or lower bound?
      if (abs(func(&m12)-y)<dy) {   // found intersection?
	if (m12Left<0) m12Left = m12;
	else m12Right = m12;
      }
      m12 += dx;
    }
  }
  /*
  if (_verboseEval) cout << "Intersections with Dalitz at: "
			 << m12Left << "," << func(&m12Left) << " and "
			 << m12Right << "," << func(&m12Right) << endl;
  */
  if (m12Left<0 && m12Right<0) return 0;
  else if (m12Left*m12Right<0) return 1;
  else return 2;
}



/// Return a TGraph with the Dalitz plot boundaries
/// Caller owns object
TGraph* BdkDalitz::drawBoundary(Double_t M, Double_t m1, Double_t m2, Double_t m3,
                                    Int_t points)
{
  if (points<=0) return 0;

  TGraph *g = new TGraph(2*points);
  if (!g) return 0;

  Double_t x1 = m12MinLimit(M,m1,m2,m3);     // determine min and max of m12
  Double_t x2 = m12MaxLimit(M,m1,m2,m3);
  
  Double_t h = (x2-x1)/points;      // spacing between points
  Double_t x = x1;
  Double_t y = 0;

  // draw lower boundary
  for (Int_t i=0; i<points; i++) {         
    y = m13Min(x,M,m1,m2,m3);
    if (y>0) g->SetPoint(i,x,y);
    x += h;
  }
  
  // draw upper boundary
  x = x2;
  for (Int_t i=points; i<2*points; i++) {  
    y = m13Max(x,M,m1,m2,m3);
    if (y>0) g->SetPoint(i,x,y);
    x -= h;
  }

  // remove points with invalid coordindates
  for (Int_t i=0; i<2*points; i++) {
    g->GetPoint(i,x,y);
    if ((x<=0) || (y<=0)) g->RemovePoint(i);
  }

  // add one more point that closes the boundary
  g->Set(g->GetN()+1);
  g->GetPoint(0,x,y);
  g->SetPoint(g->GetN()-1,x,y);

  return g;
}


/// Calculate weight of each bin in histogram "h" with respect to Dalitz boundaries.
/// Set weightContent to kTRUE if you want the original content of h to be
/// reweighted. Otherwise the weights overwrite the content of h and the resulting
/// histogram is normalized to its initial integral.
/// "step" is the step size for each axis
/// IT IS ASSUMED THAT m12 IS PLOTTED ON THE X-AXIS !!!!
Bool_t BdkDalitz::weightBins(Double_t M, Double_t m1, Double_t m2, Double_t m3,
                             TH2* h, Bool_t weightContent, Double_t step)
{
  if (!h) return kFALSE;
  if (h->GetDimension()!=2) return kFALSE;

  TAxis *ax = h->GetXaxis();
  TAxis *ay = h->GetYaxis();
  if (!ax || !ay) return kFALSE;

  Double_t integral = h->Integral();   // for normalization at the end
  // Loop over all bins
  for (Int_t i=ax->GetFirst(); i<=ax->GetLast(); i++) {
    for (Int_t j=ay->GetFirst(); j<=ay->GetLast(); j++) {
   
      // Coordinates of bin corners
      Double_t x1 = ax->GetBinLowEdge(i);
      Double_t x2 = ax->GetBinUpEdge(i);
      Double_t y1 = ay->GetBinLowEdge(j);
      Double_t y2 = ay->GetBinUpEdge(j);

      Int_t NinDalitz = 0;
      Int_t N = 0;
      for (Double_t x=x1; x<=x2; x+=step) {
	for (Double_t y=y1; y<=y2; y+=step) {
	  N++;
	  if (inDalitz(x,y,M,m1,m2,m3)) NinDalitz++;
	}
      }
      Double_t weight = (Double_t)NinDalitz/N*(x2-x1)*(y2-y1);

      if (weightContent) {
	if (weight>0) {
	  if (h->GetBinContent(i,j)==0)
	    cout << "BdkDalitz::weightBins: WARNING: Bin ("
		 << i << "," << j << ") has no entries."<<endl;
	  else h->SetBinContent(i,j,h->GetBinContent(i,j)/weight);
	}
        else h->SetBinContent(i,j,0);
      }
      else h->SetBinContent(i,j,weight);
    }  
  }
  // Normalize to the original integral
  if (weightContent) h->Scale(integral/h->Integral());
  return kTRUE;
}


/// Returns the number of events in data that are OUTSIDE the Daltiz plot.
/// Returns -1 on error.
/// Set useVetoWindow if you want to include the m23 veto window in the test
Int_t BdkDalitz::inDalitzDataSet(Double_t M, Double_t m1, Double_t m2, Double_t m3,
                                 const RooAbsData& data, 
                                 const RooRealVar& m12, const RooRealVar& m13,
                                 Bool_t useVetoWindow, Bool_t verbose)
{
  if (data.numEntries()==0) {
    if (verbose) cout << "BdkDalitz::inDalitzDataSet(): empty dataset."<<endl;
    return -1;
  }
  if (data.get(0)->find(m12.GetName())==0) {
    if (verbose) cout <<"BdkDalitz::inDalitzDataSet(): "
                      << m12.GetName() << " is not in dataset." << endl;
    return -1;
  }
  if (data.get(0)->find(m13.GetName())==0) {
    if (verbose) cout << "BdkDalitz::inDalitzDataSet(): "
                      << m13.GetName() << " is not in dataset." << endl;
    return -1;
  }


  Int_t outOfDalitz = 0;

  for (int i=0; i<data.numEntries(); i++) {
    const RooArgSet* vars = data.get(i);
    Double_t m12Val = ((RooRealVar*)vars->find(m12.GetName()))->getVal();
    Double_t m13Val = ((RooRealVar*)vars->find(m13.GetName()))->getVal();
    
    Bool_t b;
    if (useVetoWindow) b = inDalitz(m12Val, m13Val, M, m1, m2, m3);
    else b = inDalitzBounds(m12Val, m13Val, M, m1, m2, m3);

    if (!b) {
      if (verbose) cout <<"BdkDalitz::inDalitzDataSet(): m12 = " << m12Val
                        << " m13 = "<<m13Val<<" is outside of Dalitz boundary."
                        << endl;
      outOfDalitz++;
    }
  }
  return outOfDalitz;
}
