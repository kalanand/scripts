/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: Bdk2DpolyDalitz.cc,v 1.22 2007/04/04 08:58:01 fwinkl Exp $
 * Description:
 *   2D polynomial
 * History:
 *   25 Oct 2005, created, Frank Winklmeier, 
 *                adapted from Kalanand Mishra's Roo2Dpoly_hhPi0
 *
 * Copyright (C) 2005 Colorado State University and SLAC
 *****************************************************************************/

#include <iostream>
#include <fstream>

#include "BToDKTo3piK/Bdk2DpolyDalitz.hh"
#include "BToDKTo3piK/BdkMath.hh"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooAbsFunc.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooRandom.hh"
#include "RooFitCore/RooComplex.hh"
#include "RooFitCore/RooRealProxy.hh"
#include "RooFitCore/RooFormulaVar.hh"


using namespace std;
using std::abs;

ClassImp(Bdk2DpolyDalitz)
;


// Linear function representing m13 mass veto window
class BdkLinFunc : public RooAbsFunc {
public:
  BdkLinFunc(Double_t slope = 1, Double_t c = 0) :
    RooAbsFunc(1), _slope(slope), _c(c) {}

  virtual Double_t operator()(const Double_t xvector[]) const {
    return (_slope*xvector[0] + _c);
  }
  virtual Double_t getMinLimit(UInt_t dimension) const { return -1e6; }
  virtual Double_t getMaxLimit(UInt_t dimension) const { return +1e6; }

private:
  Double_t _slope;
  Double_t _c;
};


/// Constructor for polynomial with normal parametrization
/**
 The polynomial is defined over the two variables x = m12, y = m13 as:

        p(x,y) = c0 + c1*x + c2*y + c3*x^2 + c4*y^2 + 
                 c5*x*y + c6*x^3 + c7*y^3 + c8*x^2*y + c9*x*y^2
*/
Bdk2DpolyDalitz::Bdk2DpolyDalitz(const char *name, const char *title, 
                                 BdkDalitzBase::Flavor flavor, BdkDalitzBase::Mode DdecMode, 
                                 RooAbsReal& m12, RooAbsReal& m13, 
                                 RooAbsReal& c0, RooAbsReal& c1, RooAbsReal& c2, 
                                 RooAbsReal& c3, RooAbsReal& c4, RooAbsReal& c5, 
                                 RooAbsReal& c6, RooAbsReal& c7, RooAbsReal& c8, 
                                 RooAbsReal& c9) :

  BdkDalitzEff(name, title, flavor, DdecMode),
  _m12("m12","Invariant Mass square of pi0 pi+",this,m12),
  _m13("m13","Invariant Mass square of pi0 pi-",this,m13), 
  _s1(TString(GetName())+"s1","0",RooArgList()),   // those are not needed in this parametrization
  _s2(TString(GetName())+"s2","0",RooArgList()),
  _s3(TString(GetName())+"s3","0",RooArgList()),
  _s4(TString(GetName())+"s4","0",RooArgList()),
  _a1(TString(GetName())+"a1","0",RooArgList()),
  _a2(TString(GetName())+"a2","0",RooArgList()),
  _a3(TString(GetName())+"a3","0",RooArgList()),
  _a4(TString(GetName())+"a4","0",RooArgList()),
  _c0("c0","c0",this,c0),
  _c1("c1","c1",this,c1),
  _c2("c2","c2",this,c2),
  _c3("c3","c3",this,c3),
  _c4("c4","c4",this,c4),
  _c5("c5","c5",this,c5),
  _c6("c6","c6",this,c6),
  _c7("c7","c7",this,c7),
  _c8("c8","c8",this,c8),
  _c9("c9","c9",this,c9),
  _customInt(0),
  _epsRel(1e-4)
{
  setDdecMode(DdecMode);
}


/// Constructor for polynomial with a symmetric/asymmetric parametrization
/**
   The polynomial is defined over the two variables x = m12, y = m13 as:

        p(x,y) = c0 + s1(x+y) + s2(x^2+y^2) + s3(x^3+y^3) + s4(yx^2+xy^2) + s5(xy)
                      a1(x-y) + a2(x^2-y^2) + a3(x^3-y^3) + a4(yx^2-xy^2)

   Must be called with any value for saParams to distinguish it from the
   default constructor.
*/
Bdk2DpolyDalitz::Bdk2DpolyDalitz(const char *name, const char *title, 
                                 BdkDalitzBase::Flavor flavor, BdkDalitzBase::Mode DdecMode, 
                                 RooAbsReal& m12, RooAbsReal& m13, 
                                 RooAbsReal& c0, RooAbsReal& s1, RooAbsReal& s2, 
                                 RooAbsReal& s3, RooAbsReal& s4, RooAbsReal& s5, 
                                 RooAbsReal& a1, RooAbsReal& a2, RooAbsReal& a3, 
                                 RooAbsReal& a4, Bool_t saParams) :
  BdkDalitzEff(name, title, flavor, DdecMode),
  _m12("m12","Invariant Mass square of pi0 pi+",this,m12),
  _m13("m13","Invariant Mass square of pi0 pi-",this,m13),
  _s1(TString(GetName())+"s1","@0+@1",RooArgList(s1,a1)),
  _s2(TString(GetName())+"s2","@0-@1",RooArgList(s1,a1)),
  _s3(TString(GetName())+"s3","@0+@1",RooArgList(s2,a2)),
  _s4(TString(GetName())+"s4","@0-@1",RooArgList(s2,a2)),
  _a1(TString(GetName())+"a1","@0+@1",RooArgList(s3,a3)),
  _a2(TString(GetName())+"a2","@0-@1",RooArgList(s3,a3)),
  _a3(TString(GetName())+"a3","@0+@1",RooArgList(s4,a4)),
  _a4(TString(GetName())+"a4","@0-@1",RooArgList(s4,a4)),
  _c0("c0","c0",this,c0),
  _c1("c1","c1",this,_s1),
  _c2("c2","c2",this,_s2),
  _c3("c3","c3",this,_s3),
  _c4("c4","c4",this,_s4),
  _c5("c5","c5",this,s5),
  _c6("c6","c6",this,_a1),
  _c7("c7","c7",this,_a2),
  _c8("c8","c8",this,_a3),
  _c9("c9","c9",this,_a4),
  _customInt(0),
  _epsRel(1e-4)
{
  setDdecMode(DdecMode);
}


/// Copy Constructor
Bdk2DpolyDalitz::Bdk2DpolyDalitz(const Bdk2DpolyDalitz& other, const char* name) :
  BdkDalitzEff(other, name),
  _m12("m12",this,other._m12),
  _m13("m13",this,other._m13),
  _s1(other._s1),
  _s2(other._s2),
  _s3(other._s3),
  _s4(other._s4),
  _a1(other._a1),
  _a2(other._a2),
  _a3(other._a3),
  _a4(other._a4),
  _c0("c0",this,other._c0),
  _c1("c1",this,other._c1),
  _c2("c2",this,other._c2),
  _c3("c3",this,other._c3),
  _c4("c4",this,other._c4),
  _c5("c5",this,other._c5),
  _c6("c6",this,other._c6),
  _c7("c7",this,other._c7),
  _c8("c8",this,other._c8),
  _c9("c9",this,other._c9),
  _customInt(other._customInt),
  _epsRel(other._epsRel)
{
}

/// Evaluate polynomial
Double_t Bdk2DpolyDalitz::evaluate() const
{
  if (BdkDalitzBase::D0 == _flavor) {
    return evaluateAt(_m12, _m13);
  }
  // otherwise, return flipped coordinates:
  return evaluateAt(_m13, _m12);
}

/// Evaluate polynomial
Double_t Bdk2DpolyDalitz::evaluateAt(Double_t m12, Double_t m13) const
{
  // Check if inside Dalitz boundaries
  if (!inDalitz(m12,m13)) { 
    if (_verboseEval>1) cout << "("<<m12<<","<<m13<<") is outside the Dalitz plot."
                             << " Returning 0."<<endl;
    return 0.0;
  }

  // Evaluate
  Double_t value = _c0 + _c1*m12 + _c2*m13 + _c3*m12*m12 + _c4*m13*m13 + 
    _c5*m12*m13 + _c6*pow3(m12) + _c7*pow3(m13) + _c8*m12*m12*m13 + 
    _c9*m12*m13*m13 ;
    
  return value;
}


/// Decide which integral to use
Int_t Bdk2DpolyDalitz::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,const char* rangeName) const 
{
  if (_forceNumInt) return 0;                      // numerical integration

  if (matchArgs(allVars,analVars,_m12,_m13)) {     // 2D integration
    return (_customInt>0 ? _customInt : 2);        // 2 = 2D hybrid, other = custom
  }

  if (matchArgs(allVars,analVars,_m12)) return 3;  // 1D analytical
  if (matchArgs(allVars,analVars,_m13)) return 4;  // 1D analytical
   
  return 0 ;
}


/// Trapezoidal integration (from Numerical Recipes)
Double_t trapzd(const Bdk2DpolyDalitz *poly, Double_t (Bdk2DpolyDalitz::*f)(Double_t) const, Double_t a, Double_t b, Int_t n)
{
  Double_t x    ;
  static Double_t s;
  Int_t it,j;
  if (n == 1) {
    return (s = 0.5*(b-a)*((poly->*f)(a)+(poly->*f)(b)));
  } 
  else {
    it = 1;
    for (j=1; j<n-1; j++) {
      it <<= 1;    // multiply by 2
    }
    Double_t del = (b-a)/it;
    x = a + 0.5*del;
    Double_t sum = 0;
    for (Int_t j=1; j<=it; j++) {
      x += del;
      sum += (poly->*f)(x);
    }
    s = 0.5*(s+(b-a)*sum/it );
    return s;
  }
}

/// Adaptive numerical integration (from Numerical Recipes)
Double_t qtrap(const Bdk2DpolyDalitz *poly, Double_t (Bdk2DpolyDalitz::*f)(Double_t) const, Double_t a, Double_t b)
{
  Double_t s;
  Double_t olds = 0.0;

  const Int_t JMIN = 5;      // minumum number of refinements
  const Int_t JMAX = 100;    // maximum number of refinements
  const Double_t EPS = poly->epsRel();

  for (Int_t j=1;j<=JMAX;j++) {
    s = trapzd(poly,f,a,b,j);
    if (j > JMIN)
      if (abs(s-olds) < EPS*abs(olds) || (s == 0.0 && olds == 0.0)) {
        //        cout << j <<" integrations steps needed."<<endl;
        return s;
      }
    olds=s;
  }
  cout << "Too many steps in routine qtrap" << endl;
  return 0.0;
}


/// Integral of polynomial over m12 evaluated at (m12,m13)
Double_t Bdk2DpolyDalitz::intPolym12(Double_t m12, Double_t m13) const
{
  return (_c0 + _c2*m13 + _c4*m13*m13 + _c7*m13*m13*m13)*m12 +
         (_c1 + _c5*m13 + _c9*m13*m13)*m12*m12/2 + 
         (_c3 + _c8*m13)*pow3(m12)/3 +
         (_c6)*pow4(m12)/4;
}


/// Integral of polynomial over m13 evaluated at (m12,m13)
Double_t Bdk2DpolyDalitz::intPolym13(Double_t m12, Double_t m13) const
{
  return (_c0 + _c1*m12 + _c3*m12*m12 + _c6*m12*m12*m12)*m13 +
         (_c2 + _c5*m12 + _c8*m12*m12)*m13*m13/2 + 
         (_c4 + _c9*m12)*pow3(m13)/3 +
         (_c7)*pow4(m13)/4;
}


/// Integrate polynomial over m13 within the Dalitz region
Double_t Bdk2DpolyDalitz::intPolym13Dalitz(Double_t m12) const
{
  // Dalitz boundaries
  Double_t max = m13Max(m12);
  Double_t min = m13Min(m12);  
  if (max==0 || min==0) return 0;  // not in Dalitz region

  Double_t result = intPolym13(m12,max) - intPolym13(m12,min);

  if (hasM23Veto()) {
    // Veto window boundaries in m13
    Double_t m13VetoMin = mOtherVal(m12, m23VetoMax());
    Double_t m13VetoMax = mOtherVal(m12, m23VetoMin());

    // Is veto window part of Dalitz region?
    if (m13VetoMin>max || m13VetoMax<min) return result;

     // Adjust integration range if necessary
    if (min<m13VetoMin) min = m13VetoMin;
    if (max>m13VetoMax) max = m13VetoMax;

    // Integrate veto window
    Double_t vetoInt = intPolym13(m12,max) - intPolym13(m12,min);
    result -= vetoInt;  
    if (_verboseEval>1) cout << "Veto window integration over m13 at m12 = "<<m12
                             <<" from "<< min <<" to " << max << " = "<<vetoInt<<"."
                             << " Final result = "<<result<<endl;
  }
  return result;
}


/// Integrate polynomial over m12 within the Dalitz region
Double_t Bdk2DpolyDalitz::intPolym12Dalitz(Double_t m13) const
{
  // Dalitz boundaries
  Double_t max = m12Max(m13);
  Double_t min = m12Min(m13);  
  if (max==0 || min==0) return 0;  // not in Dalitz region

  Double_t result = intPolym12(max,m13) - intPolym12(min,m13);

  if (hasM23Veto()) {
    // Veto window boundaries in m13
    Double_t m12VetoMin = mOtherVal(m13, m23VetoMax());
    Double_t m12VetoMax = mOtherVal(m13, m23VetoMin());

    // Is veto window part of Dalitz region?
    if (m12VetoMin>max || m12VetoMax<min) return result;

     // Adjust integration range if necessary
    if (min<m12VetoMin) min = m12VetoMin;
    if (max>m12VetoMax) max = m12VetoMax;

    // Integrate veto window
    Double_t vetoInt = intPolym12(max,m13) - intPolym12(min,m13);
    result -= vetoInt;  
    if (_verboseEval>1) cout << "Veto window integration over m12 at m13 = "<<m13
                             <<" from "<< min <<" to " << max << " = "<<vetoInt<<"."
                             << " Final result = "<<result<<endl;
  }
  return result;
}


/// Alternative purely anlytical 2D integration of mass veto window
/// Currently not used.
Double_t Bdk2DpolyDalitz::intM23Veto(Double_t m12Min, Double_t m12Max) const
{
  Double_t k = _M*_M + _m1*_m1 + _m2*_m2 + _m3*_m3;

  Double_t h = k - m23VetoMin();
  Double_t g = k - m23VetoMax();

  Double_t a = -1;
  Double_t b = -1;
  BdkLinFunc linFunc(-1, h);  
  Int_t i = findIntersection(linFunc,a,b);
  if (i!=2) {
    cout << "WARNING: Cannot integrate m23 mass veto window because only "
	 << i << "intersection(s) with Dalitz boundary were found."<<endl;
    return -1;
  }
  if (m12Min>a) a = m12Min;
  if (m12Max<b) b = m12Max;

  Double_t result = 1/4*(_c7+_c8-_c9-_c6)*(g-h)*(pow4(b)-pow4(a))
    + 1/3*((3/2*_c7+_c8/2-_c9)*(h*h-g*g)+(_c3+_c4-_c5)*(h-g))*(pow3(b)-pow3(a))
    + 1/2*((_c7-_c9/3)*(pow3(g)-pow3(h))+(_c4-_c5/2)*(g*g-h*h)+(_c2-_c1)*(g-h))*(b*b-a*a)
    + (_c7/4*(pow4(h)-pow4(g))+_c4/3*(pow3(h)-pow3(g))+_c2/2*(h*h-g*g)+_c0*(h-g))*(b-a);

  return result;
}


Double_t Bdk2DpolyDalitz::analyticalIntegral(Int_t code,const char* rangeName) const 
{
  Double_t result = 0;
  
  Double_t max1 = _m12.max(rangeName);
  Double_t min1 = _m12.min(rangeName);

  if (code==1) {       // 2D analytical integration (with some approximations)
    
    // Here used to be Kalanand's aproximated analytical integral
    // I removed it since it wasn't sufficient. 
    // Go back to CVS revision 1.14 to retrieve it.

    cout << "Analytical integral not supported."<<endl;

    if (_verboseEval) cout << "Doing analytical integration from "
                           <<min1<<" to "<<max1 << " = "<< result<<endl;
  }

  else if (code==2) {   // 2D hybrid integration (m13 analytical, m12 numerical)

    result = qtrap(this, &Bdk2DpolyDalitz::intPolym13Dalitz,min1,max1);

    if (_verboseEval) cout << "Doing custom numerical integration from "
                           <<min1<<" to "<<max1 <<" = "<<result<<endl;
  }
  else if (code==3) {  // 1D analytical integration over m12

    result = intPolym12Dalitz(_m13);

    if (_verboseEval) cout << "Doing 1D analytical integration over m12 from "
                           << min1 <<" to " << max1 << " = " <<result<<endl;
 
  }
  else if (code==4) {  // 1D analytical integration over m13

    result = intPolym13Dalitz(_m12);
  
    if (_verboseEval) cout << "Doing 1D analytical integration over m13 from "
                           << min1 <<" to " << max1 << " = " <<result<<endl;

  }
  else assert(0);

  return result;
} 

