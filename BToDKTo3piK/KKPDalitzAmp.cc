/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: RooFitDalitz
 * Authors:
 *   BL, Ben Lau, Princeton, yanpan@slac.stanford.edu
 *   Abi Soffer, Colorado State, abi@slac.stanford.edu
 *   Kalanand Mishra, U. Cincinnati, kalanand@slac.stanford.edu
 *  NOTE: Only works when particles 2 and 3 are CP-conjugates of each other
 *****************************************************************************/

// -- CLASS DESCRIPTION --
// Dalitz plot amplitude for D0->K+K-pi0
//
//Note: The KKPDalitzAmp class consturct the pdf for D0->pi0 K+ K-
//Here I specific the phase convection:
//s12 -> pi0 K+ 
//s13 -> pi0 K- 


#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>

#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooListProxy.hh"
#include "RooFitCore/RooComplex.hh"
#include "BToDKTo3piK/KKPDalitzAmp.hh"
#include "RooFitCore/RooRandom.hh"
#include "RooFitCore/RooRealVar.hh"
#include "BToDKTo3piK/BdkDalitzBase.hh"

using namespace std;


ClassImp(KKPDalitzAmp);

// Static members:
std::vector<KKPDalitzAmp *> KKPDalitzAmp::_registry;
const Double_t KKPDalitzAmp::F_D = 1.0;
const Double_t KKPDalitzAmp::DEGTORAD = M_PI/180.0;
const Bool_t KKPDalitzAmp::_enforceTransversality = kTRUE;


void KKPDalitzAmp::normalizeAll(int nEvents){ 
  // Normalize all objects
  cout << "----- KKPDalitzAmp::normalizeAll: normalizing " << _registry.size()
       << " objects with " << nEvents << " events each." << endl
       << "Note: Only objects whose PDFs have been created are normalized."
       << endl
       << "To force the creation of a PDF call its wrapper's parameters()." 
       << endl;

  for (int i = 0; i < _registry.size(); ++i){
    cout << "Normalizing " << _registry[i]->GetName() << endl;
    _registry[i]->calNorm(nEvents);
    _registry[i]->normParams().Print("V");
  }
}


void KKPDalitzAmp::registerObject(KKPDalitzAmp * amp) {
  cout << "KKPDalitzAmp::Registry: registering " << amp->GetName() << endl;
  _registry.push_back(amp);
}

void KKPDalitzAmp::deregisterObject(KKPDalitzAmp * amp) {
  remove(_registry.begin(),_registry.end(),amp);
}

void KKPDalitzAmp::cleanRegistry() {
  _registry.clear();
}


// Constructor
KKPDalitzAmp::KKPDalitzAmp(const char * name, 
                             const char * title, 
                             BdkDalitzBase * pdf,
                             Int_t componentsBit,
                             Int_t defaultKstarSpin) :
  TNamed(name, title),
  _componentsBit(componentsBit),
  _defaultKstarSpin(defaultKstarSpin),
  _pdf(pdf)
{
  registerObject(this);
  initResonance();
  registerParams(_pdf);
}

// link our parameters to the pdf:
void KKPDalitzAmp::registerParams(BdkDalitzBase * pdf) {
  RooListProxy * proxyList = 
    new RooListProxy(TString(GetName()) + ".proxyList",
		     TString(GetTitle()) + " proxyList",
		     pdf);
  
  TIterator*  tIter = _params.createIterator();
  RooAbsArg*  coef;
  while(coef = (RooAbsArg*)tIter->Next()) {
    if (!dynamic_cast<RooAbsReal*>(coef)) assert(0);
    proxyList->add(*coef);
  }
  delete tIter;

  // also link the normalization parameters:
  tIter = _normParams.createIterator();
  while(coef = (RooAbsArg*)tIter->Next()) {
    if (!dynamic_cast<RooAbsReal*>(coef)) assert(0);
    proxyList->add(*coef);
  }
  delete tIter;
}


// Destructor
KKPDalitzAmp::~KKPDalitzAmp() 
{
  deregisterObject(this);
}

Double_t KKPDalitzAmp::efficiency(Double_t m12, Double_t m13) const
{
//      if(m12>1.0 && m13>1.0) return 0.0;
//   Double_t m23 =  1.8645*1.8645+ 0.1349766*0.1349766 
//     + 2.0*0.493677*0.493677 - m12 -m13;
//   if(!(m23>1.0 && m23<1.2)) return 0.0;

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
  //  return _pdf->efficiency(m12, m13);
}


//Return the normalization
Double_t KKPDalitzAmp::getNormalization() const 
{  
  Double_t norm = 0;

  for(int j=0;j<_nComps;j++) {

    Double_t amp_j = _ampRes[j]->getVal();
    if (amp_j==0) continue;

    RooComplex coeff_j(amp_j * cos(_phaseRes[j]->getVal()*DEGTORAD),
                       amp_j * sin(_phaseRes[j]->getVal()*DEGTORAD));

    for(int k=0;k<_nComps;k++){
   
      Double_t amp_k = _ampRes[k]->getVal();
      if (amp_k==0) continue;

      RooComplex coeff_k(amp_k * cos(_phaseRes[k]->getVal()*DEGTORAD),
                         amp_k * sin(_phaseRes[k]->getVal()*DEGTORAD));

      RooComplex normCoeff(_normReal[j][k]->getVal(), 
			   _normImag[j][k]->getVal());

      norm += ((coeff_j*coeff_k.conj()) * normCoeff).re();
    }
  }  
  return norm;
}

//Return the normalization
Double_t KKPDalitzAmp::unNormalizedFitFraction(int j) const 
{  
  RooComplex coeff_j(_ampRes[j]->getVal() * cos(_phaseRes[j]->getVal()*DEGTORAD),
                     _ampRes[j]->getVal() * sin(_phaseRes[j]->getVal()*DEGTORAD));
  
  RooComplex normCoeff(_normReal[j][j]->getVal(), 
		       _normImag[j][j]->getVal());
  
  Double_t fraction = ((coeff_j * coeff_j.conj()) * normCoeff).re();
  return fraction;
}

// Returns the fit fractions. Causes a memory leak, but who cares... 
RooArgSet KKPDalitzAmp::fitFractions() const {
  Double_t norm = getNormalization();
  RooArgSet result;
  for(int j=0;j<_nComps;j++) {  
    TString name("FitFraction_");
    name += _nameRes[j];
    RooRealVar * var = new RooRealVar(name, name, 
				      unNormalizedFitFraction(j) / norm);
    result.add(*var);
  }
  return result;
}

void KKPDalitzAmp::calNorm(int events, Bool_t enforceNormalizeAll) 
{
  if (enforceNormalizeAll) {
    // Caller is enforcing normalization of all components
    setNeedToNormalizeAll();

    cout << GetName() << ": KKPDalitzAmp::calNorm(): performing MC integration on all components" << endl;
  }

  if (kFALSE == needToNormalizeAny()) {
    // No need to normalize. Note that this won't happen if enforceNormalizeall
    // is true.
    return;
  }

  if (verbose().Contains("n")) {
    cout << GetName() << ": start normalization" << endl;;
  }

  Double_t m12gen, m13gen;
  RooComplex normArray[SIZE][SIZE]; 

  for(int i=1;i<events;i++) {    
    m12gen = RooRandom::uniform()*2.0;
    m13gen = RooRandom::uniform()*2.0;
    
    if (inDalitz(m12gen,m13gen)) {

      for(int j=0;j<_nComps;j++) {
	// Get the Breit-Wigner and spin terms
	Double_t Spin_j = SpinFactor(m12gen, m13gen, _spinRes[j], 
				     _trackinfo[j],j);
	
	RooComplex BW_j = RelativisticBW(m12gen, m13gen, 
					 _massRes[j]->getVal(), 
					 _gammaRes[j]->getVal(),j);

	// use symmetry of normalization constants: N_jk = (N_kj)*
	for(int k=j;k<_nComps;k++){
	  //check if component k or j needs to be normalized:
	  if (needToNormalize(k) || needToNormalize(j)) {
	    
	    Double_t Spin_k = SpinFactor(m12gen, m13gen, _spinRes[k], 
					 _trackinfo[k],k);
	    
	    RooComplex BW_k = RelativisticBW(m12gen, m13gen, 
					     _massRes[k]->getVal(), 
					     _gammaRes[k]->getVal(),k);
	    
	    RooComplex integral = (BW_j*BW_k.conj()*Spin_j*Spin_k) * 
	      efficiency(m12gen,m13gen);
	    
	    normArray[j][k] = normArray[j][k] + integral;
	  }
	}
      }
    }  //end the inDalitz if loop

    //give a status report when user asked for normalization:
    if(enforceNormalizeAll && fmod(i,events/10.0)==0) {
      cout << "Done with " << i << " events" << endl;
    }
  } //end the for loop

  // Convert the temporary normArray to RRV's that can be read in like
  // any other PDF parameter:
  for(int j=0;j<_nComps;j++) {
    for(int k=j;k<_nComps;k++){
      if (needToNormalize(k) || needToNormalize(j)) {
	_normReal[j][k]->setVal((normArray[j][k]*4.0/events).re());
	_normImag[j][k]->setVal((normArray[j][k]*4.0/events).im());

	// use symmetry of normalization constants: N_jk = (N_kj)*
        _normReal[k][j]->setVal(_normReal[j][k]->getVal());
        _normImag[k][j]->setVal((-1)*_normImag[j][k]->getVal());
        
      }
    }
  }

  if (verbose().Contains("n")) {
    cout << GetName() << ": end normalization" << endl;
  }

  setIsNormalized();
}





RooComplex KKPDalitzAmp::getamp(Double_t m12, Double_t m13) const
{

  //This return the f(m+^2,m-^2)
  if (kFALSE == inDalitz(m12, m13)) {
    return RooComplex(0.0, 0.0);
  }

  // normalize (mutable operation) if need to:
  ((KKPDalitzAmp*)this)->calNorm(1000000, kFALSE); 

  RooComplex dalitzamplitude(0,0);

  for (int i=0; i<_nComps; i++) {

    Double_t amp = _ampRes[i]->getVal();
    if (amp==0) continue;

    RooComplex coeff = RooComplex(amp * cos(_phaseRes[i]->getVal()*DEGTORAD),
                                  amp * sin(_phaseRes[i]->getVal()*DEGTORAD));    

    RooComplex matrixelement = RelativisticBW(m12, m13, 
					      _massRes[i]->getVal(), 
					      _gammaRes[i]->getVal(), i)
                               * SpinFactor(m12, m13, _spinRes[i], _trackinfo[i], i);

    dalitzamplitude = dalitzamplitude + coeff*matrixelement;
  }
  
  // In the PDF we will take abs2() for the amplitude, so need to
  // take sqrt() for efficiency map now:
  return (dalitzamplitude * sqrt(efficiency(m12,m13)));
}

// initialize one component:
void KKPDalitzAmp::addComp(const char * name, double amp, double phase, 
			    double mass, double width, 
			    ResDaughters trackInfo, 
			    Int_t spin, int parSource) {
  
  _nameRes[_nComps] = TString(name); 
  
  _ampRes[_nComps] = new RooRealVar(TString(GetName()) + "."
				    + _nameRes[_nComps] + "_amp", 
				    TString(GetTitle()) + "."
				    + _nameRes[_nComps] + "_amp", 
				    amp);
  
  _phaseRes[_nComps] = new RooRealVar(TString(GetName()) + "."
				      + _nameRes[_nComps] + "_phase", 
				      TString(GetTitle()) + "."
				      + _nameRes[_nComps] + "_phase", 
				      phase);
  
  if (0 > parSource) { // make new mass and width parameters:
    _massRes[_nComps] = new RooRealVar(TString(GetName()) + "."
				       + _nameRes[_nComps] + "_mass", 
				       TString(GetTitle()) + "."
				       + _nameRes[_nComps] + "_mass", mass);  
    
    _gammaRes[_nComps] = new RooRealVar(TString(GetName()) + "."
					+ _nameRes[_nComps] + "_width",
					TString(GetTitle()) + "."
					+ _nameRes[_nComps] + "_width", width);
  }
  else { // the mass and width for this resonance come from parSource:
    _massRes[_nComps] = _massRes[parSource];
    _gammaRes[_nComps] = _gammaRes[parSource];
  }

  // put these on the parameters list:
  _params.addOwned(*_ampRes[_nComps]);
  _params.addOwned(*_phaseRes[_nComps]);
  _params.addOwned(*_massRes[_nComps]);
  _params.addOwned(*_gammaRes[_nComps]);

  // User's responsibility to read normalization parameters or
  // normalize after initialization. Hence:
  _massResLast[_nComps] = _massRes[_nComps]->getVal();
  _gammaResLast[_nComps] = _gammaRes[_nComps]->getVal();

  // set spin and daughters code:
  _spinRes[_nComps] = spin;
  _trackinfo[_nComps] = trackInfo;

  // update # of components:
  ++_nComps;  
}
  

// initialize the components:
void KKPDalitzAmp::initResonance() 
{
  _nComps = 0;
  
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
    KstarIndex = _nComps;  // store its index
    addComp("K*+", 1.0, 0.0, 0.89166, 0.0508, PI0_KP, _defaultKstarSpin);
  }

  if (KSTAR1410P & componentsBit()) {
    Kstar1410PIndex = _nComps; // store its index
    addComp("K*1410+", 1.0, 61.0, 1.414, 0.232, PI0_KP, _defaultKstarSpin);
  }

  if (NONRESP & componentsBit()) {
    nonresIndex = _nComps;  // store its index
    addComp("K+pi0_SW", 1.5, 120.0, 1.414, 0.290, PI0_KP, SPIN0);
  }

  // K+K- resonances:
  if (PHI & componentsBit()) {
    addComp("Phi", 0.7, 43.0, 1.0195, 0.00426, KP_KM, SPIN1);
  }

  if (F2P1525 & componentsBit()) {
    addComp("F2P1525", 0.6, -40.0, 1.525, 0.073, KP_KM, SPIN2);
  }

  if (F0 & componentsBit()) {
    addComp("F0", 0.9, 170, 0.965, 0.0, KP_KM, SPIN0);

    // width for the F0 is taken from Flatte' formula and BES parameterization
    // which looked at J/psi-->phi pi+pi- and J/psi-->phi K+K-
  }

  if (A0 & componentsBit()) {
    addComp("A0", 0.9, 0, 0.998, 0.0, KP_KM, SPIN0);
  }

  // K-pi0 resonances:
  if (KSTARM & componentsBit()) {
    addComp("K*-", 0.67, -31.0, 0.89166, 0.0508, KM_PI0, _defaultKstarSpin,
	    KstarIndex);
  }

  if (KSTAR1410M & componentsBit()) {
    addComp("K*1410-", 1.0, 148.0, 1.414, 0.232, KM_PI0, _defaultKstarSpin,
	    Kstar1410PIndex);
  }

   if (NONRESM & componentsBit()) {
    addComp("K-pi0_SW", 1.5, -45.0, 1.414, 0.290, KM_PI0, SPIN0,nonresIndex);
  }
 



  // Masses of the above daughters
  _mDaug[0]=0.0;
  _mDaug[1]=0.1349766;
  _mDaug[2] = 0.493677;
  _mDaug[3] = 0.493677;
  _mD = 1.8645; //Mass of D0-meson


  // create and store all RooRealVar normalization parameters:
  for (int j = 0; j < _nComps; ++j) {
    for (int k = 0; k < _nComps; ++k) {

      TString suffix;
      suffix += j;
      suffix += "_";
      suffix += k;

      TString nameReal(".normReal_");
      nameReal += suffix;

      TString nameImag(".normImag_");
      nameImag += suffix;

      _normReal[j][k] = new RooRealVar(GetName() + nameReal, 
				       GetTitle() + nameReal,
				       1.0);

      _normImag[j][k] = new RooRealVar(GetName() + nameImag, 
				       GetTitle() + nameImag,
				       1.0);

      _normParams.addOwned(*_normReal[j][k]);  
      _normParams.addOwned(*_normImag[j][k]); 
    }
  }

  // create the radial parameter for the Blatt-Weisskopf factor
  _resRadius = new RooRealVar(TString(GetName()) + ".resRadius",
                              TString(GetTitle()) + " resonance radius",
                              0.0);
  _params.addOwned(*_resRadius);


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


Double_t KKPDalitzAmp::kinematics(Double_t m12, Double_t m13, Int_t n) const {
  Double_t value=0;
  
  //compute m23(invarient mass square) using energy and momentum conservation
  Double_t m23 = _mD*_mD + _mDaug[1]*_mDaug[1] + _mDaug[2]*_mDaug[2] + 
    _mDaug[3]*_mDaug[3] - m12 -m13;
  
  //Now need to determine which variables to return using the info from Int_t n,(bachler tracks)
  //convention n=1 -> m23, n=2 -> m13, n=3 -> m12
  
  switch (n) {
  case 1:
    value=m23;
    break;
  case 2:
    value=m13;
    break;
  case 3:
    value=m12;
    break;    
  }
  return value;
}


Double_t KKPDalitzAmp::dampingFactorSquared(Double_t s, Int_t n, Int_t spin) const
{ 
  Double_t p = pionCMmom(s,n);
  Double_t R = _resRadius->getVal();

  switch (spin) {
  case 0: return 1;
  case 1: return (1 + R*R*p*p);
  case 2: return (9 + 3*R*R*p*p + pow(R,4)*pow(p,4));
  default: return 1;
  }
}

Double_t KKPDalitzAmp::FrEval(Double_t s, Double_t m0, Int_t n, Int_t spin) const
{
  //  if (!aboveThreshold( s )) return 1;

  Double_t Fr2_s = dampingFactorSquared(s, n, spin);
  Double_t Fr2_m = dampingFactorSquared(m0*m0, n, spin);

  return sqrt(Fr2_s/Fr2_m);
}



RooComplex KKPDalitzAmp::RelativisticBW(Double_t m12, Double_t m13, 
					 Double_t m0, Double_t width, 
                                         Int_t i) const 
{
  //notice m12 and m13 are invarient mass "square" of the decay product

  //s is the invariant mass "square" M^2(ab)
  //m0 is the mass of the resonace

 Double_t s=kinematics(m12,m13,_trackinfo[i]);

 if( fabs(_massRes[i]->getVal()-0.970)<0.015 ) return Flatte(m0,s); // this is f0(980)
 else if( (fabs(_massRes[i]->getVal()-0.998)<0.010) && (0 == _spinRes[i]) ) return a0(s); // this is a0(980)
 else if (0 == _spinRes[i]) {

    /*********  Breit_Wigner:     return RooComplex(1.0,0.0); *********/  
    /*********  LASS parameters:  return LASS(s, i);      *************/ 
    /*********  E791 parameters:  return E791Table(s);    *************/ 
    
    if(_E791->getVal()>0.0) return E791Table(s);
    else return LASS(s, i);
  }
  else{    
    Double_t A = (m0*m0 - s);
    Double_t B = m0*runningWidth(m12,m13,m0,_gammaRes[i]->getVal(),_spinRes[i],_trackinfo[i]);
    Double_t C = A*A+B*B;
    
  assert (C != 0);

  Double_t Fr = FrEval(s, m0 ,_trackinfo[i], _spinRes[i]);

  return RooComplex(A/C, B/C)*Fr*F_D;
  }
}



//replace Breit-Wigner with Flatte for f0(980)

RooComplex KKPDalitzAmp::Flatte(Double_t m0, Double_t s) const
{
  Double_t gPiPi = 0.165;        // +/- 0.010 +/- 0.015 GeV/c^2
  Double_t gKK   = gPiPi*4.21;   // +/- 0.25 +/- 0.21
  Double_t mK    = 0.493677;     //This is K
  Double_t mPi0  = 0.1349766;    //This is pi0 
  Double_t mPi   = 0.13957;      //This is pi 
  Double_t mK0s  = 0.497648;     //This is K0 

  Double_t mSumSq0_ = 4.0*mPi0*mPi0;
  Double_t mSumSq1_ = 4.0*mPi*mPi;
  Double_t mSumSq2_ = 4.0*mK*mK;
  Double_t mSumSq3_ = 4.0*mK0s*mK0s;

  Double_t dMSq = m0*m0 - s;

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
  

  Double_t widthTerm = gPiPi*rho1*m0 + gKK*rho2*m0;
  
  RooComplex resAmplitude = RooComplex(dMSq, widthTerm);
  
  Double_t denomFactor = dMSq*dMSq + widthTerm*widthTerm;
  
  Double_t invDenomFactor = 0.0;
  if (denomFactor > 1e-10) {invDenomFactor = 1.0/denomFactor;}
  
  resAmplitude = resAmplitude*(invDenomFactor);
  
  return resAmplitude;
}



RooComplex KKPDalitzAmp::a0(Double_t s) const
{

  Double_t m0 = 0.999;
  Double_t gEtaPi = 0.105;        
  Double_t gKK   = 0.102;   
  Double_t mK    = 0.493677;     //This is K
  Double_t mPi0  = 0.1349766;    //This is pi0 
  Double_t mEta   = 0.54775;      //This is eta 
  Double_t mK0s  = 0.497648;     //This is K0 

  Double_t mSumSq1_ = 4.0*mEta*mPi0;
  Double_t mSumSq2_ = 4.0*mK*mK;
  Double_t mSumSq3_ = 4.0*mK0s*mK0s;

  Double_t dMSq = m0*m0 - s;

  Double_t rho1(0.0), rho2(0.0);

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


  Double_t widthTerm = gEtaPi*rho1*m0 + gKK*rho2*m0;
  
  RooComplex resAmplitude = RooComplex(dMSq, widthTerm);
  
  Double_t denomFactor = dMSq*dMSq + widthTerm*widthTerm;
  
  Double_t invDenomFactor = 0.0;
  if (denomFactor > 1e-10) {invDenomFactor = 1.0/denomFactor;}
  
  resAmplitude = resAmplitude*(invDenomFactor);
  
  return resAmplitude;
}




//replace Breit-Wigner with LASS

// RooComplex KKPDalitzAmp::LASS(Double_t s, Int_t i) const 
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


RooComplex KKPDalitzAmp::LASS(Double_t s, Int_t i) const
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

RooComplex KKPDalitzAmp::E791Table(Double_t s) const 
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




//Try to introduce Relativistic Breit-Wigner function

inline Double_t KKPDalitzAmp::betaKin( Double_t s, Double_t m ) const
{
  Double_t rad = 1 - m*m/s;
  return (rad >= 0.0) ? sqrt( rad ) : 1;  
}


// For A -> B + C calculate momentum of B and C in rest frame of A 
// s is center of mass energy squared (mass of A)
// see (38.16) in PDG04
inline Double_t KKPDalitzAmp::pionCMmom( Double_t s, Double_t m1, Double_t m2 ) const
{
  return 0.5 * sqrt(s) * betaKin(s, m1+m2) * betaKin(s, m1-m2);
}

inline Double_t KKPDalitzAmp::pionCMmom(Double_t s, Int_t n) const
{
 /*  
   index   resonance   dtr particles   dtr indices
     1       phi        K+ K-         23
     2       K*-        pi0 K-         13
     3       K*+        pi0 K+         12
  */
  switch (n) {
  case 1: return pionCMmom(s, _mDaug[2], _mDaug[3]);
  case 2: return pionCMmom(s, _mDaug[3], _mDaug[1]);
  case 3: return pionCMmom(s, _mDaug[1], _mDaug[2]);
  default: return 0;
  }
}

Double_t KKPDalitzAmp::runningWidth(Double_t m12, Double_t m13, Double_t m0,
                                     Double_t width, Int_t spin, Int_t n) const
{
  //  if (!aboveThreshold( s )) return width;

  Double_t s = kinematics(m12,m13,n);  // invariant mass squared of daugthers

  Double_t k_s  = pionCMmom(s, n);
  Double_t k_m0 = pionCMmom(m0*m0, n);   // using nominal resonance mass

  Double_t Fr = FrEval(s, m0, n, spin);
  return width * m0/sqrt(s) * pow(k_s/k_m0, 2*spin+1) * Fr*Fr;
}




Double_t KKPDalitzAmp::SpinFactor(Double_t m12, Double_t m13, Int_t spin, 
                                   Int_t n, Int_t i) const
{
  Double_t value;
    
  if (_trackinfo[i]==0) {
    value = 1.0;   //take care the NonResonace Component is uniform over dalitz region
  }
  else {
    
    Double_t m23;
    Double_t a1, b1, c1;
    Double_t mAC=0, mBC=0, mAB=0;
    Double_t _mA=0, _mB=0, _mC=0;
    
    //Reminder:
    //_mDaug[1]=Pi0
    //_mDaug[2]=K+
    //_mDaug[3]=K-
    
    
    //compute m23(invarient mass square) using energy and momentum conservation
    m23 = _mD*_mD + _mDaug[1]*_mDaug[1] + _mDaug[2]*_mDaug[2] + _mDaug[3]*_mDaug[3] - m12 -m13;    
    switch (n) {  //n is the _trackinfo 
    case 1: //K+ K- is resonace pair
      _mA=_mDaug[2]; //K+
      _mB=_mDaug[3]; //K-
      _mC=_mDaug[1]; //Pi0
      mAC=m12;
      mBC=m13;
      mAB=m23;
      break;
    case 2: //Pi0 K- is the resonance pair (Cabibbo allow decay)
      _mA=_mDaug[3]; //K-
      _mB=_mDaug[1]; //pi0
      _mC=_mDaug[2]; //K+
      mAC=m23;
      mBC=m12;
      mAB=m13;
      break;
    case 3: //pi0 K+ is the resonace pair (Cabibbo suppress decay)
      _mA=_mDaug[1]; //Pi0
      _mB=_mDaug[2]; //K-
      _mC=_mDaug[3]; //K+
      mAC=m13;
      mBC=m23;
      mAB=m12;
      break;
    }

    double massFactor = 1./pow(_massRes[i]->getVal(),2);
    if (_enforceTransversality) {
      massFactor = 1./mAB;
    }
    
    switch (spin) {
    case 0: 
      value = 1;
      break;
    case 1:      //spin 1 case
      value = ((mAC-mBC) + (massFactor*(_mD*_mD - _mC*_mC)*(_mB*_mB-_mA*_mA)));
      break;
    case 2:
      a1 = (pow(((mBC-mAC) + (massFactor*(_mD*_mD - _mC*_mC)*(_mA*_mA-_mB*_mB))),2));
      b1 = ((mAB-(2*_mD*_mD)-(2*_mC*_mC))+massFactor*pow((_mD*_mD-_mC*_mC),2));
      c1 = ((mAB-(2*_mA*_mA)-(2*_mB*_mB))+massFactor*pow((_mA*_mA-_mB*_mB),2));      
      value =  a1-((1/3.0)*b1*c1);
      break;
    default:
      cout << "KKPDalitzAmp::SpinFactor() cannot handle spin>2 resonances."<<endl;
      value = 1;
      assert(spin<3);
    }
  }

  return value;
}



Bool_t KKPDalitzAmp::inDalitz(Double_t m12, Double_t m13) const
{
  return _pdf->inDalitz(m12, m13);
}

// Check if component j needs to be normalized:
Bool_t KKPDalitzAmp::needToNormalize(int j) const {
  if (_massRes[j]->getVal() != _massResLast[j] ||
      _gammaRes[j]->getVal() != _gammaResLast[j]) {
    return kTRUE;
  }
  return kFALSE;
}

Bool_t KKPDalitzAmp::needToNormalizeAny() const {
  Bool_t result = kFALSE;

  for(int j=0;j<_nComps;j++) {
    if (needToNormalize(j)) {
      result = kTRUE;
      
      if (verbose().Contains("n")) {
	cout << GetName() << ": need to normalize " << _nameRes[j] << endl;
      }
      else {
	// break out of the loop if not printing what needs to be
	// normalized. If printing, then finish the loop in order to
	// print everything.
	break;  
      }
    }
  }

  return result;
}	

// Reset normalization-needed flags: 
void KKPDalitzAmp::setIsNormalized() {
  for(int j=0;j<_nComps;j++) {
    _massResLast[j] = _massRes[j]->getVal();
    _gammaResLast[j] = _gammaRes[j]->getVal();
  }
}

// Change the last mass so that the component is flagged as needing to
// be normalized:
void KKPDalitzAmp::setNeedToNormalize(int component) {
  _massResLast[component] = _massRes[component]->getVal() + 1;
}

void KKPDalitzAmp::setNeedToNormalizeAll() {
  for(int j=0;j<_nComps;j++) {
    setNeedToNormalize(j);
  }
}









