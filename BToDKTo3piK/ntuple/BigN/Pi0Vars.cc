//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: Pi0Vars.cc,v 1.1 2005/11/08 17:29:10 zhangjl Exp $
//
// Description:
//	Class Pi0Vars
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Abi Soffer                  Original Author
//
// Copyright Information:
//	Copyright (C) 2003
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------
#include "BetaMiniUser/Pi0Vars.hh"

//-------------
// C Headers --
//-------------
#include <assert.h>
#include <math.h>

//---------------
// C++ Headers --
//---------------
#include <iostream.h>
#include <math.h>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "Beta/BtaCandidate.hh"
#include "Beta/EventInfo.hh"
#include "ErrLogger/ErrLog.hh"
#include "HepTuple/Tuple.h"
#include "BetaCoreTools/BtaMcAssoc.hh"
#include "PDT/Pdt.hh"
#include "PDT/PdtEntry.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//----------------
// Constructors --
//----------------

// in general, a module constructor should not do much.  The begin(job) or
// begin(run) members are better places to put initialization
Pi0Vars::Pi0Vars(const BtaCandidate & pi0, 
	  const HepAList<BtaCandidate> & pi0List,
	  const EventInfo & eventInfo,
	  const HepAList<BtaCandidate> * gammaList,
	  const BtaMcAssoc * truthMap) : 
  _pi0(&pi0),
  _pi0List(&pi0List),
  _gammaList(gammaList),
  _truthMap(truthMap), 
  _boostVector(-(eventInfo.cmFrame().boostVector())),
  _cmFrame(-(eventInfo.cmFrame())),
  _mom(-100),
  _mass(-100),
  _helic(-100),
  _good(-100)
{
  for (int i = 0; i < NPAIRS; ++i) {
    _massH[i] = -100;
    _massS[i] = -100;
    _helicH[i] = -100;
    _helicS[i] = -100;
    _asymH[i] = -100;
    _asymS[i] = -100;
    _momH[i] = -100;
    _momS[i] = -100;
    _goodH[i] = -100;
    _goodS[i] = -100;
    _tMassH[i] = -100;
    _tMassS[i] = -100;
    _tHelicH[i] = -100;
    _tHelicS[i] = -100;
    _tGoodH[i] = -100;
    _tGoodS[i] = -100;
  }

  calc();
}  

//--------------
// Destructor --
//--------------

Pi0Vars::~Pi0Vars( )
{
}


//-----------
// Helpers --
//-----------


const BtaCandidate * 
Pi0Vars::gammaRecoFromMc(const BtaCandidate * mc) const {
  int i;
  if (0 != _gammaList && 0 != mc) {
    for (i = 0; i < _gammaList->length(); ++i) {
      const BtaCandidate * gamma = (*_gammaList)[i];
      const BtaCandidate * trueGamma = mcFromReco(gamma);
      if (0 != trueGamma && mc->p4() == trueGamma->p4()) {
	return gamma;
      }
    }
  }
    
  return 0;
}

//-------------------------------------------------------------------------
const BtaCandidate * 
Pi0Vars::mcFromReco(const BtaCandidate * reco) const {
  if (0 != _truthMap) {
    return _truthMap->mcFromReco(reco);
  }
  return 0;
}

//-------------------------------------------------------------------------
const BtaCandidate *  
Pi0Vars::good2Body(const BtaCandidate * mother, const PdtEntry * entry) const {

  const BtaCandidate * dtr1;
  const BtaCandidate * dtr2;

  const BtaCandidate * dtr = 0;    
  HepAListIterator<BtaCandidate> dtrIter = mother->daughterIterator();        
  int d = 0;
  while (dtr = dtrIter()){
    ++d;
    if (1 == d){
      dtr1 = dtr;
    } 
    else {
      dtr2 = dtr;
    }
  }
   
  return good2Body(dtr1, dtr2, entry);
}
 
//-------------------------------------------------------------------------
const BtaCandidate * 
Pi0Vars::good2Body(const BtaCandidate * g1, 
		   const BtaCandidate * g2, 
		   const PdtEntry * entry) const {
  
  const BtaCandidate * t1 = mcFromReco(g1);
  const BtaCandidate * t2 = mcFromReco(g2);

  if (0 != t1 && 0 != t2) {
    const BtaCandidate * m1 = t1->theMother();
    const BtaCandidate * m2 = t2->theMother();
    if (0 != m1 && m1 == m2 && m1->pdtEntry() == entry) {
      return m1;
    }
  }
  return 0;
}

//-------------------------------------------------------------------------
double Pi0Vars::calcHelicityCM(HepLorentzVector daughter,    // by value
			       const HepLorentzVector & mother) const {

  // Transform the grandmother and the daughter to the mother's frame:
  const Hep3Vector boostVector = -(mother.boostVector());
  daughter.boost(boostVector);
  HepLorentzVector grandmother = _cmFrame;
  grandmother.boost(boostVector);
  return cos(grandmother.angle(daughter));
}

//------------------------------------------------------------------------
double Pi0Vars::dtrMass(const BtaCandidate * mother) const {
  HepLorentzVector p4;
  const BtaCandidate * dtr = 0;    
  HepAListIterator<BtaCandidate> dtrIter = mother->daughterIterator();        
  while (dtr = dtrIter()){
    p4 += dtr->p4();
  }

  return  p4.mag();
}
 
//------------------------------------------------------------------------
void Pi0Vars::bestPair(const BtaCandidate * include,  // input
		       const BtaCandidate * exclude1,  // input
		       const BtaCandidate * exclude2,  // input
		       float & bestMass,
		       float & helic,
		       int & good) const {
		     
  if (0 == include) {
    return;
  }

  BtaCandidate * partner;
  float pairMom;
  float asym;

  bestPair(include, exclude1, exclude2, partner, 
	   bestMass, pairMom, helic, asym, good);
}

//------------------------------------------------------------------------
void Pi0Vars::bestPair(const BtaCandidate * include,  // input
		       const BtaCandidate * exclude1,  // input
		       const BtaCandidate * exclude2,  // input
		       BtaCandidate *& partner,
		       float & bestMass,
		       float & pairMom,
		       float & helic,
		       float & asym,
		       int & good) const {
		     
  static const double pi0Mass =  0.134976;

  partner = 0;
  double bestMassDiff = 20.0;

  for (int p = 0; p < _pi0List->length(); ++p){
    BtaCandidate * pi0 = (*_pi0List)[p];
    
    // Find its daughters:
    BtaCandidate * tempPartner = 0;
    BtaCandidate * dtr = 0;
    bool isMother = false;
    HepAListIterator<BtaCandidate> dtrIter = pi0->daughterIterator();        
    while (dtr = dtrIter()){	
      if ((0 != exclude1 && dtr->p4() == exclude1->p4()) || 
	  (0 != exclude2 && dtr->p4() == exclude2->p4())){
	isMother = false; // ensure rejection below
	break; // don't use a pi0 which decays into the original partner
      }
      
      if (dtr->p4() != include->p4()){
	tempPartner = dtr; // got sister of include
      }
      else { // pi0 is a mother of include
	isMother = true;
      }
    }
    
    if (false == isMother) {
      continue; // this pi0 is not a mother of include
    }
    
    double mass = dtrMass(pi0);
    double massDiff = mass - pi0Mass;
    
    if (fabs(massDiff) < fabs(bestMassDiff)) {
      bestMassDiff = massDiff;
      partner = tempPartner;
      bestMass = mass;
    }
  }
  
  if (0 != partner) {
    // fill the necessary quantities:
    HepLorentzVector partnerP4CM = partner->p4();
    partnerP4CM.boost(_boostVector);

    HepLorentzVector includeP4CM = include->p4();
    includeP4CM.boost(_boostVector);

    HepLorentzVector pairP4CM = includeP4CM + partnerP4CM;
    pairMom = pairP4CM.vect().mag();
    
    helic = calcHelicityCM(includeP4CM, pairP4CM);

    asym = (includeP4CM.e() - partnerP4CM.e()) / 
      (includeP4CM.e() + partnerP4CM.e());
    
    good = 0;
    static const PdtEntry * pi0Entry = Pdt::lookup("pi0");
    if (0 != good2Body(include, partner, pi0Entry)) {
      good = 1;
    }
  }
}   
 
//------------------------------------------------------------------------
void Pi0Vars::calc() {
  HepLorentzVector pi0P4 = _pi0->p4();
  pi0P4.boost(_boostVector);
  
  // Get pi0 momentum:
  _mom = pi0P4.vect().mag();
  
  // get daughters and their momentum:
  BtaCandidate * dtrH = 0;
  BtaCandidate * dtrS = 0;
  
  BtaCandidate * dtr = 0;    
  int d = 0;
  HepAListIterator<BtaCandidate> dtrIter = _pi0->daughterIterator();        
  while (dtr = dtrIter()){
    if (0 == d){
      dtrH = dtr;
    }
    else {
      dtrS = dtr;
    }
    ++d;
  }
  
  HepLorentzVector pH = dtrH->p4();
  HepLorentzVector pS = dtrS->p4();
  
  pH.boost(_boostVector);
  pS.boost(_boostVector);
  
  // determine who is more energetic:
  if (pH.e() > pS.e()) {
    dtr = dtrS;
    dtrS = dtrH;
    dtrH = dtr;
    
    HepLorentzVector p = pS;
    pS = pH;
    pH = p;
  }
  
  // get truth:
  _good = 0;
  static const PdtEntry * pi0Entry = Pdt::lookup("pi0");
  if (0 != good2Body(dtrH, dtrS, pi0Entry)) {
    _good = 1;
  }
  
  _mass = (pH + pS).mag();
  
  // the helicity angle.  Note that all vectors are in the CM frame:
  _helic = calcHelicityCM(pH, pi0P4);
  
  // find the best candidate pairs:
  BtaCandidate * partnerH1 = 0;

  bestPair(dtrH, dtrS, 0, 
	   partnerH1, 
	   _massH[0], 
	   _momH[0], 
	   _helicH[0], 
	   _asymH[0],
	   _goodH[0]);
  
  bestPair(partnerH1, dtrH, 0, 
	   _tMassH[0], 
	   _tHelicH[0], 
	   _tGoodH[0]);
  
  
  BtaCandidate * partnerS1 = 0;
  
  bestPair(dtrS, dtrH, 0, 
	   partnerS1,
	   _massS[0], 
	   _momS[0], 
	   _helicS[0], 
	   _asymS[0],
	   _goodS[0]);
  
  bestPair(partnerS1, dtrS, 0, 
	     _tMassS[0], 
	   _tHelicS[0], 
	   _tGoodS[0]);
  
  // find the 2nd best candidate pairs:
  BtaCandidate * partnerH2 = 0;
  
  bestPair(dtrH, dtrS, partnerH1, 
	   partnerH2, 
	   _massH[1], 
	   _momH[1], 
	   _helicH[1], 
	   _asymH[1],
	   _goodH[1]);    
  
  bestPair(partnerH2, dtrH, 0, 
	   _tMassH[1], 
	   _tHelicH[1], 
	   _tGoodH[1]);
  
  
  BtaCandidate * partnerS2 = 0;
  
  bestPair(dtrS, dtrH, partnerS1, 
	   partnerS2, 
	   _massS[1], 
	   _momS[1], 
	   _helicS[1], 
	   _asymS[1],
	   _goodS[1]);
  
  bestPair(partnerS2, dtrS, 0, 
	   _tMassS[1], 
	   _tHelicS[1], 
	   _tGoodS[1]);
}
 
