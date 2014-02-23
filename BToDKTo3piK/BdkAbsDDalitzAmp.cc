/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: RooFitDalitz
 * Authors:
 *   BL, Ben Lau, Princeton, yanpan@slac.stanford.edu
 *   Abi Soffer, Colorado State, abi@slac.stanford.edu
 *  NOTE: Only works when particles 2 and 3 are CP-conjugates of each other
 *****************************************************************************/

// -- CLASS DESCRIPTION --
// Dalitz plot amplitude for D0->pi+pi-pi0
//
//Note: The BdkAbsDDalitzAmp class consturct the pdf for D0->pi0 pi+ pi-
//Here I specific the phase convection:
//s12 -> pi0 pi+ 
//s13 -> pi0 pi- 


#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>

#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooListProxy.hh"
#include "RooFitCore/RooComplex.hh"
#include "BToDKTo3piK/BdkAbsDDalitzAmp.hh"
#include "RooFitCore/RooRandom.hh"
#include "RooFitCore/RooRealVar.hh"
#include "BToDKTo3piK/BdkDalitzBase.hh"
#include "BToDKTo3piK/BdkMath.hh"

using namespace std;


ClassImp(BdkAbsDDalitzAmp);

// Static members:
std::vector<BdkAbsDDalitzAmp *> BdkAbsDDalitzAmp::_registry;
const Double_t BdkAbsDDalitzAmp::_mD = 1.8645; //Mass of D0-meson   
const Double_t BdkAbsDDalitzAmp::F_D = 1.0;
const Double_t BdkAbsDDalitzAmp::DEGTORAD = M_PI/180.0;
const Bool_t BdkAbsDDalitzAmp::_enforceTransversality = kTRUE;


void BdkAbsDDalitzAmp::normalizeAll(int nEvents){ 
  // Normalize all objects
  cout << "----- BdkAbsDDalitzAmp::normalizeAll: normalizing " << _registry.size()
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


void BdkAbsDDalitzAmp::registerObject(BdkAbsDDalitzAmp * amp) {
  cout << "BdkAbsDDalitzAmp::Registry: registering " << amp->GetName() << endl;
  _registry.push_back(amp);
}

void BdkAbsDDalitzAmp::deregisterObject(BdkAbsDDalitzAmp * amp) {
  remove(_registry.begin(),_registry.end(),amp);
}

void BdkAbsDDalitzAmp::cleanRegistry() {
  _registry.clear();
}


// Constructor
BdkAbsDDalitzAmp::BdkAbsDDalitzAmp(const char * name, 
                                   const char * title, 
                                   BdkDalitzBase * pdf,
                                   Int_t componentsBit) :  
  TNamed(name, title),
  _pdf(pdf),
  _componentsBit(componentsBit),
  _useFixedWidth(kFALSE),  
  _useAutoNorm(kFALSE)
{
  registerObject(this);
  _nComps = 0;
}


// Set the external pdf
// On request, register all parameters of this dalitzAmp with pdf
void BdkAbsDDalitzAmp::setPdf(BdkDalitzBase * pdf, Bool_t registerParamsWithPdf)
{
  _pdf = pdf;
  if (registerParamsWithPdf) registerParams(_pdf);
}


// link our parameters to the pdf:
void BdkAbsDDalitzAmp::registerParams(BdkDalitzBase * pdf)
{
  if (pdf==0) {
    //    cout << GetName()<< ".registerParams(): Cannot register parameters for NULL pdf."
    //         << endl;
    return;
  }
  
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
BdkAbsDDalitzAmp::~BdkAbsDDalitzAmp() 
{
  deregisterObject(this);
}


//Return the normalization
Double_t BdkAbsDDalitzAmp::getNormalization() const 
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
Double_t BdkAbsDDalitzAmp::unNormalizedFitFraction(int j) const 
{  
  RooComplex coeff_j(_ampRes[j]->getVal() * cos(_phaseRes[j]->getVal()*DEGTORAD),
                     _ampRes[j]->getVal() * sin(_phaseRes[j]->getVal()*DEGTORAD));
  
  RooComplex normCoeff(_normReal[j][j]->getVal(), 
		       _normImag[j][j]->getVal());
  
  Double_t fraction = ((coeff_j * coeff_j.conj()) * normCoeff).re();
  return fraction;
}


// Returns the fit fractions. Causes a memory leak, but who cares... 
RooArgSet BdkAbsDDalitzAmp::fitFractions() const {
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





//Return the normalization integral over the isospin amplitudes
// Valid only for D0-->pi+pi-pi0 decay.
 
RooArgSet BdkAbsDDalitzAmp::IntegralOverIsospin() const {    

  RooComplex ampcoeff(0.0, 0.0);
  RooComplex normcoeff(0.0, 0.0);
  RooComplex amp(0.0, 0.0);
  RooComplex amp_NR(0.0, 0.0);
  RooComplex amp_rhoP(0.0, 0.0);
  RooComplex amp_rho0(0.0, 0.0);
  RooComplex amp_rhoM(0.0, 0.0);
  RooComplex amp_F(0.0, 0.0);

  for(int j=0;j<_nComps;j++) {   

    ampcoeff = RooComplex(_ampRes[j]->getVal() * 
			  cos(_phaseRes[j]->getVal()*DEGTORAD),
			  _ampRes[j]->getVal() * 
			  sin(_phaseRes[j]->getVal()*DEGTORAD));  
    normcoeff = RooComplex(_normReal[j][j]->getVal(), 
			   _normImag[j][j]->getVal());

    amp = ampcoeff * sqrt(normcoeff.re());
  

    if (_nameRes[j].Contains("Nonres",TString::kExact)) { 
      amp_NR = amp_NR + amp; 
    }

    if(_nameRes[j].Contains("Rho+",TString::kExact) || 
       _nameRes[j].Contains("Rho2s+",TString::kExact) || 
       _nameRes[j].Contains("Rho1700+",TString::kExact)) { 
      amp_rhoP = amp_rhoP + amp; 
    }

    if(_nameRes[j].Contains("Rho-",TString::kExact) || 
       _nameRes[j].Contains("Rho2s-",TString::kExact) || 
       _nameRes[j].Contains("Rho1700-",TString::kExact)) { 
      amp_rhoM = amp_rhoM + amp; 
    }

    if(_nameRes[j].Contains("Rho0",TString::kExact) || 
       _nameRes[j].Contains("Rho2s0",TString::kExact) || 
       _nameRes[j].Contains("Rho17000",TString::kExact)) { 
      amp_rho0 = amp_rho0 + amp; 
    }

    if(_nameRes[j].Contains("F0",TString::kExact) || 
       _nameRes[j].Contains("F0_1370",TString::kExact) || 
       _nameRes[j].Contains("F0_1500",TString::kExact) ||
       _nameRes[j].Contains("F0_1710",TString::kExact) ||
       _nameRes[j].Contains("F2",TString::kExact) ||
       _nameRes[j].Contains("Sigma",TString::kExact)) { 
      amp_F = amp_F + amp; 
    }
  }

  RooArgSet normisospin;
  Double_t norm = getNormalization();

  RooRealVar* var1 = new RooRealVar("isospin-norm-NR", "isospin-norm-NR",
				    amp_NR.abs()/ sqrt(norm));
  normisospin.add(*var1);

  RooRealVar* var2 = new RooRealVar("isospin-norm-rhoP", "isospin-norm-rhoP", 
				   amp_rhoP.abs()/  sqrt(norm));
  normisospin.add(*var2);

  RooRealVar* var3 = new RooRealVar("isospin-norm-rho0", "isospin-norm-rho0", 
				    amp_rho0.abs()/ sqrt(norm));
  normisospin.add(*var3);

  RooRealVar* var4 = new RooRealVar("isospin-norm-rhoM", "isospin-norm-rhoM", 
				    amp_rhoM.abs()/ sqrt(norm));
  normisospin.add(*var4);

  RooRealVar* var5 = new RooRealVar("isospin-norm-F", "isospin-norm-F", 
				    amp_F.abs()/ sqrt(norm));
  normisospin.add(*var5);

  RooRealVar* pha1 = new RooRealVar("isospin-phase-NR", "isospin-phase-NR",
				    TMath::ATan(amp_NR.im()/amp_NR.re())/DEGTORAD);
  normisospin.add(*pha1);

  RooRealVar* pha2 = new RooRealVar("isospin-phase-rhoP","isospin-phase-rhoP", 
				   TMath::ATan(amp_rhoP.im()/amp_rhoP.re())/DEGTORAD);
  normisospin.add(*pha2);

  RooRealVar* pha3 = new RooRealVar("isospin-phase-rho0","isospin-phase-rho0", 
				    TMath::ATan(amp_rho0.im()/amp_rho0.re())/DEGTORAD);
  normisospin.add(*pha3);

  RooRealVar* pha4 = new RooRealVar("isospin-phase-rhoM","isospin-phase-rhoM", 
				    TMath::ATan(amp_rhoM.im()/amp_rhoM.re())/DEGTORAD);
  normisospin.add(*pha4);

  RooRealVar* pha5 = new RooRealVar("isospin-phase-F","isospin-phase-F", 
				    TMath::ATan(amp_F.im()/amp_F.re())/DEGTORAD);
  normisospin.add(*pha5);

  return normisospin;
}




//Return the unnormalization integral over the Breit-Wigner amplitude only:
// Integral_r =  \int |BW_r(s+, s-)|^2 ds+ ds- 
Double_t BdkAbsDDalitzAmp::IntegralOverBreitWigner(int j) const 
{    
  RooComplex coeff(_normReal[j][j]->getVal(), _normImag[j][j]->getVal());
  
  Double_t integral = coeff.abs();
  return integral;
}


// Returns the square-root of the normalization integral over the 
// Breit-Wigner amplitude: sqrt( Integral_r / norm )
RooArgSet BdkAbsDDalitzAmp::BreitWignerNormalizationCoefficients() const {
  Double_t norm = getNormalization();
  RooArgSet fraction;
  for(int j=0;j<_nComps;j++) {  
    Double_t val = sqrt(IntegralOverBreitWigner(j) / norm);
    RooRealVar * var = new RooRealVar( TString("BWFraction_")+ _nameRes[j], 
				       "", val);
    fraction.add(*var);
  }
  return fraction;
}




void BdkAbsDDalitzAmp::calNorm(int events, Bool_t enforceNormalizeAll) 
{
  if (enforceNormalizeAll) {
    // Caller is enforcing normalization of all components
    setNeedToNormalizeAll();

    cout << GetName() 
         << ": BdkAbsDDalitzAmp::calNorm(): performing MC integration on all components" 
         << endl;
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

  Double_t m12Range = m12NormRange();
  Double_t m13Range = m13NormRange();
  
  for(int i=0;i<events;i++) {    
    m12gen = RooRandom::uniform()*m12Range;
    m13gen = RooRandom::uniform()*m13Range;
    
    if (inDalitz(m12gen,m13gen)) {

      for(int j=0;j<_nComps;j++) {
	// Get the Breit-Wigner and spin terms
        RooComplex me_j = matrixElement(m12gen, m13gen, j);
        
	// use symmetry of normalization constants: N_jk = (N_kj)*
	for(int k=j;k<_nComps;k++){
	  //check if component k or j needs to be normalized:
	  if (needToNormalize(k) || needToNormalize(j)) {
            RooComplex me_k = matrixElement(m12gen, m13gen, k);
            RooComplex integral = (me_j*me_k.conj()) * efficiency(m12gen,m13gen);
                    
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
  
  const Double_t normArea = m12Range*m13Range;
  
  for(int j=0;j<_nComps;j++) {
    for(int k=j;k<_nComps;k++){
      if (needToNormalize(k) || needToNormalize(j)) {
	_normReal[j][k]->setVal((normArray[j][k]*normArea/events).re());
	_normImag[j][k]->setVal((normArray[j][k]*normArea/events).im());

	// use symmetry of normalization constants: N_jk = (N_kj)*
        _normReal[k][j]->setVal(_normReal[j][k]->getVal());
        _normImag[k][j]->setVal((-1)*_normImag[j][k]->getVal());
        
      }
    }
  }

  setIsNormalized();

  // Calculate normalization for D/Dbar interference
  calDDbarNorm(events);

  if (verbose().Contains("n")) {
    cout << GetName() << ": end normalization" << endl;
  }
}

// Helper function to calcualte the normalization integral of this
// amplitude with a second amplitude given in amp2.
// If amp==0 then use m12<->m13 flipped amplitude as amp2.
// Return vector of normalization coefficient p0,...,p3.
TVectorD BdkAbsDDalitzAmp::calDIntNorm(BdkAbsDDalitzAmp* amp2, int nEvents) const
{
  Double_t p0 = 0; 
  Double_t p1 = 0; 
  Double_t p2 = 0; 
  Double_t p3 = 0; 
 
  Double_t m12Range = m12NormRange();
  Double_t m13Range = m13NormRange();
  
  for(Double_t i=0;i<nEvents;i++) {       
    Double_t m12gen = RooRandom::uniform()*m12Range;
    Double_t m13gen = RooRandom::uniform()*m13Range;
    
    if(inDalitz(m12gen,m13gen)) {


      RooComplex Damp = getamp(m12gen, m13gen);
      RooComplex Dbaramp(0,0);

      // amplitude calculated with appropriately charge-flipped variables
      // or using second amplitude
      if (amp2) Dbaramp = amp2->getamp(m12gen, m13gen);
      else Dbaramp = getamp(m13gen,m12gen);

      p0 += ((Dbaramp)*(Damp.conj())).re();
      p1 += ((Dbaramp)*(Damp.conj())).im();
      p2 += Dbaramp.abs2();
      p3 += Damp.abs2();
      
    }  //end the if loop
  } //end the for loop

  const Double_t normArea = m12Range*m13Range;
    
  p0 = p0*normArea/nEvents;
  p1 = p1*normArea/nEvents;
  p2 = p2*normArea/nEvents;
  p3 = p3*normArea/nEvents;

  TVectorD v(4);
  v[0] = p0;
  v[1] = p1;
  v[2] = p2;
  v[3] = p3;
  return v;
}


// Calculates the four normalization coefficients of
// integral( |D + z Dbar|^2 ) where Dbar(m12,m13) = D(m13,m12)
// Therefore this only works for C eigenstates. This is a design flaw.
// This routine should be in BdkDKDalitz instead.
void BdkAbsDDalitzAmp::calDDbarNorm(int nEvents)
{
  cout << GetName() 
       << ": performing MC integration for D/Dbar interference" << endl;

  // Use ourself as intefering amplitude
  TVectorD p = calDIntNorm(0, nEvents);

  cout << "Precision of MC integration:"<<endl;
  cout << "p1/p0      = "<< p[1]/p[0] << endl;
  cout << "(p2-p3)/p2 = "<< (p[2]-p[3])/p[2] << endl;

  // Copy to RooRealVars
  _normReDDbar->setVal(p[0]);      
  _normImDDbar->setVal(0.0);         // Theoretically, p1 = 0
  _normDbarSqr->setVal((p[2]+p[3])/2);   // Theoretically, p2 = p3
  _normDSqr->setVal((p[2]+p[3])/2);
}


// Return matrix element of resonance i
RooComplex BdkAbsDDalitzAmp::matrixElement(Double_t m12, Double_t m13, Int_t i) const
{
  RooComplex me = RelativisticBW(m12, m13, 
				 _massRes[i]->getVal(), _gammaRes[i]->getVal(), i)
    * SpinFactor(m12, m13, _spinRes[i], _trackinfo[i], i);

  return me;
}


RooComplex BdkAbsDDalitzAmp::getamp(Double_t m12, Double_t m13) const
{
  //This return the f(m+^2,m-^2)
  if (kFALSE == inDalitz(m12, m13)) {
    return RooComplex(0.0, 0.0);
  }

  if (_useAutoNorm) {
    // normalize (mutable operation) if need to:
    ((BdkAbsDDalitzAmp*)this)->calNorm(1000000, kFALSE); 
  }
  
  RooComplex dalitzamplitude(0,0);

  for (int i=0; i<_nComps; i++) {

    Double_t amp = _ampRes[i]->getVal();
    if (amp==0) continue;

    RooComplex coeff = RooComplex(amp * cos(_phaseRes[i]->getVal()*DEGTORAD),
                                  amp * sin(_phaseRes[i]->getVal()*DEGTORAD));    

    dalitzamplitude = dalitzamplitude + coeff*matrixElement(m12,m13,i);
  }
  
  // In the PDF we will take abs2() for the amplitude, so need to
  // take sqrt() for efficiency map now:

  return (dalitzamplitude * sqrt(efficiency(m12,m13)));
}

// initialize one component:
void BdkAbsDDalitzAmp::addComp(const char * name, double amp, double phase, 
                               double mass, double width, 
                               ResDaughters trackInfo, 
                               Int_t spin, Int_t parSource, Int_t resType) {
  
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

  // type of resonance
  // can be used by derived classes to implement other matrix elements
  _typeRes[_nComps] = resType;

  // update # of components:
  ++_nComps;  
}

  
  
// initialize the components:
void BdkAbsDDalitzAmp::createParams() 
{
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

  // Add the D/Dbar inteference normalization parameters

  _normDSqr = new RooRealVar(TString(GetName())+".normDSqr", 
                             TString(GetTitle())+".normDSqr", 0);

  _normDbarSqr = new RooRealVar(TString(GetName())+".normDbarSqr", 
                                TString(GetTitle())+".normDbarSqr", 0);

  _normReDDbar = new RooRealVar(TString(GetName())+".normReDDbar", 
                                TString(GetTitle())+".normReDDbar", 0);

  _normImDDbar = new RooRealVar(TString(GetName())+".normImDDbar", 
                                TString(GetTitle())+".normImDDbar", 0);

  _normParams.addOwned(*_normDSqr);
  _normParams.addOwned(*_normDbarSqr);
  _normParams.addOwned(*_normReDDbar);
  _normParams.addOwned(*_normImDDbar);
  

  // create the radial parameter for the Blatt-Weisskopf factor
  _resRadius = new RooRealVar(TString(GetName()) + ".resRadius",
                              TString(GetTitle()) + " resonance radius",
                              1.5);
  _params.addOwned(*_resRadius);
  
}


Double_t BdkAbsDDalitzAmp::kinematics(Double_t m12, Double_t m13, Int_t n) const {
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


Double_t BdkAbsDDalitzAmp::dampingFactorSquared(Double_t s, Int_t n, Int_t spin) const
{ 
  Double_t p = pionCMmom(s,n);
  Double_t R = _resRadius->getVal();

  switch (spin) {
  case 0: return 1;
  case 1: return (1 + R*R*p*p);
  case 2: return (9 + 3*R*R*p*p + R*R*R*R*p*p*p*p);
  default: return 1;
  }
}

Double_t BdkAbsDDalitzAmp::FrEval(Double_t s, Double_t m0, Int_t n, Int_t spin) const
{
  //  if (!aboveThreshold( s )) return 1;

  Double_t Fr2_s = dampingFactorSquared(s, n, spin);
  Double_t Fr2_m = dampingFactorSquared(m0*m0, n, spin);

  return sqrt(Fr2_s/Fr2_m);
}


RooComplex BdkAbsDDalitzAmp::RelativisticBW(Double_t m12, Double_t m13, 
                                            Double_t m0, Double_t width, 
                                            Int_t i) const 
{
  //notice m12 and m13 are invarient mass "square" of the decay product

  //s is the invariant mass "square" M^2(ab)
  //m0 is the mass of the resonace

  Double_t s = kinematics(m12,m13,_trackinfo[i]);
  if (0 == _trackinfo[i] || 0 == _gammaRes[i]->getVal()) {
    // Take care the NonResonace Component is uniform over dalitz region.
    // Using the width to flag this works for pwave as well as swave.
    return RooComplex(1.0,0.0);   
  }
  else{    
    Double_t A = (m0*m0 - s);
    Double_t B = m0*runningWidth(m12,m13,m0,_gammaRes[i]->getVal(),_spinRes[i],_trackinfo[i]);
    Double_t C = A*A+B*B;
    
    assert (C != 0);

    Double_t Fr = FrEval(s, m0 ,_trackinfo[i], _spinRes[i]);

    return RooComplex(A/C, B/C)*Fr*F_D;
    //      ^^ = 1/(A-iB)
  }
}


//Try to introduce Relativistic Breit-Wigner function

inline Double_t BdkAbsDDalitzAmp::betaKin( Double_t s, Double_t m ) const
{
  Double_t rad = 1 - m*m/s;
  return (rad >= 0.0) ? sqrt( rad ) : 1;  
}


// For A -> B + C calculate momentum of B and C in rest frame of A 
// s is center of mass energy squared (mass of A)
// see (38.16) in PDG04
inline Double_t BdkAbsDDalitzAmp::pionCMmom( Double_t s, Double_t m1, Double_t m2 ) const
{
  return 0.5 * sqrt(s) * betaKin(s, m1+m2) * betaKin(s, m1-m2);
}

inline Double_t BdkAbsDDalitzAmp::pionCMmom(Double_t s, Int_t n) const
{
 /*  
   index   resonance   dtr particles   dtr indices
     1       rho0        pi+ pi-         23
     2       rho-        pi0 pi-         13
     3       rho+        pi0 pi+         12
  */
  switch (n) {
  case 1: return pionCMmom(s, _mDaug[2], _mDaug[3]);
  case 2: return pionCMmom(s, _mDaug[3], _mDaug[1]);
  case 3: return pionCMmom(s, _mDaug[1], _mDaug[2]);
  default: return 0;
  }
}

Double_t BdkAbsDDalitzAmp::runningWidth(Double_t m12, Double_t m13, Double_t m0,
                                        Double_t width, Int_t spin, Int_t n) const
{
  if (_useFixedWidth) return width;

  Double_t s = kinematics(m12,m13,n);  // invariant mass squared of daugthers

  Double_t k_s  = pionCMmom(s, n);
  Double_t k_m0 = pionCMmom(m0*m0, n);   // using nominal resonance mass

  Double_t Fr = FrEval(s, m0, n, spin);

  // The following loop replaces this: r = pow(k_s/k_m0, 2*spin+1)
  // Note: x^(2*n+1) = x*(x*x)^n
  Double_t x = k_s/k_m0;
  Double_t r = x;
  for (int i=0; i<spin; i++) r *= (x*x);
    
  return width * m0/sqrt(s) * r * Fr*Fr;
}


Double_t BdkAbsDDalitzAmp::SpinFactor(Double_t m12, Double_t m13, Int_t spin, 
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
    //_mDaug[2]=pi+
    //_mDaug[3]=pi-
    
    
    //compute m23(invarient mass square) using energy and momentum conservation
    m23 = _mD*_mD + _mDaug[1]*_mDaug[1] + _mDaug[2]*_mDaug[2] + _mDaug[3]*_mDaug[3] - m12 -m13;    
    switch (n) {  //n is the _trackinfo 
    case 1: //pi+ pi- is resonace pair
      _mA=_mDaug[2]; //pi+
      _mB=_mDaug[3]; //pi-
      _mC=_mDaug[1]; //Pi0
      mAC=m12;
      mBC=m13;
      mAB=m23;
      break;
    case 2: //Pi0 pi- is the resonance pair (Cabibbo allow decay)
      _mA=_mDaug[3]; //Pi-
      _mB=_mDaug[1]; //pi0
      _mC=_mDaug[2]; //pi+
      mAC=m23;
      mBC=m12;
      mAB=m13;
      break;
    case 3: //pi0 pi+ is the resonace pair (Cabibbo suppress decay)
      _mA=_mDaug[1]; //Pi0
      _mB=_mDaug[2]; //pi-
      _mC=_mDaug[3]; //pi+
      mAC=m13;
      mBC=m23;
      mAB=m12;
      break;
    }

    double massFactor = 1./(_massRes[i]->getVal()*_massRes[i]->getVal());
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
      a1 = sqr((mBC-mAC) + (massFactor*(_mD*_mD - _mC*_mC)*(_mA*_mA-_mB*_mB)));
      b1 = (mAB-(2*_mD*_mD)-(2*_mC*_mC))+massFactor*sqr(_mD*_mD-_mC*_mC);
      c1 = (mAB-(2*_mA*_mA)-(2*_mB*_mB))+massFactor*sqr(_mA*_mA-_mB*_mB);
      value =  a1-((1/3.0)*b1*c1);
      break;
    default:
      cout << "BdkAbsDDalitzAmp::SpinFactor() cannot handle spin>2 resonances."<<endl;
      value = 1;
      assert(spin<3);
    }
  }

  return value;
}


Bool_t BdkAbsDDalitzAmp::inDalitz(Double_t m12, Double_t m13) const
{
  return _pdf->inDalitz(m12, m13);
}

// Check if component j needs to be normalized:
Bool_t BdkAbsDDalitzAmp::needToNormalize(int j) const {

  double minDiff = 1e-4;
  if (_gammaRes[j]->getVal()>0) minDiff = _gammaRes[j]->getVal() / 20.;
    
  if (fabs(_massRes[j]->getVal() - _massResLast[j]) > minDiff ||
      fabs(_gammaRes[j]->getVal() - _gammaResLast[j]) > minDiff) {
    return kTRUE;
  }
  return kFALSE;
}

Bool_t BdkAbsDDalitzAmp::needToNormalizeAny() const {
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
void BdkAbsDDalitzAmp::setIsNormalized() {
  for(int j=0;j<_nComps;j++) {
    _massResLast[j] = _massRes[j]->getVal();
    _gammaResLast[j] = _gammaRes[j]->getVal();
  }
}

// Change the last mass so that the component is flagged as needing to
// be normalized:
void BdkAbsDDalitzAmp::setNeedToNormalize(int component) {
  _massResLast[component] = _massRes[component]->getVal() + 1;
}

void BdkAbsDDalitzAmp::setNeedToNormalizeAll() {
  for(int j=0;j<_nComps;j++) {
    setNeedToNormalize(j);
  }
}

void BdkAbsDDalitzAmp::printDDbarNorm() const
{
  RooArgSet set(*_normDSqr, *_normDbarSqr, 
                *_normReDDbar, *_normImDDbar);

  set.Print("v");
}

// Return list of amplitudes (caller owns)
RooArgList* BdkAbsDDalitzAmp::ampRes() const
{
  RooArgList* list = new RooArgList(TString(GetName())+"_amp");
  for (int i=0; i<_nComps; i++) list->add(*ampRes(i));
  return list;
}

// Return list of phases (caller owns)
RooArgList* BdkAbsDDalitzAmp::phaseRes() const
{
  RooArgList* list = new RooArgList(TString(GetName())+"_phase");
  for (int i=0; i<_nComps; i++) list->add(*phaseRes(i));
  return list;
}

// Return list of masses (caller owns)
RooArgList* BdkAbsDDalitzAmp::massRes() const
{
  RooArgList* list = new RooArgList(TString(GetName())+"_mass");
  for (int i=0; i<_nComps; i++) list->add(*massRes(i));
  return list;
}

// Return list of widths (caller owns)
RooArgList* BdkAbsDDalitzAmp::gammaRes() const
{
  RooArgList* list = new RooArgList(TString(GetName())+"_gamma");
  for (int i=0; i<_nComps; i++) list->add(*gammaRes(i));
  return list;
}

// Change from running to fixed widths:
void BdkAbsDDalitzAmp::setUseFixedWidth(Bool_t val) {
  _useFixedWidth = val;
  cout << "*** " << GetName() << ".setUseFixedWidth() called. ***" << endl
       << "*** May need to renormalize! *** " << endl;
}
