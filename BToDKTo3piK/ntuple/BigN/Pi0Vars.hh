//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: Pi0Vars.hh,v 1.1 2005/11/08 17:29:10 zhangjl Exp $
//
// Description:
//	Class Pi0Vars
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Abi Soffer                   Original Author
//
// Copyright Information:
//	Copyright (C) 2003		
//
//------------------------------------------------------------------------

#ifndef PI0VARS_HH
#define PI0VARS_HH


//----------------------
// Base Class Headers --
//----------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class HepTuple;
class BtaMcAssoc;
class EventInfo;
class HepLorentzVector;

#include <assert.h>
#include "Beta/BtaCandidate.hh"
#include "CLHEP/Alist/AList.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"


//		---------------------
// 		-- Class Interface --
//		---------------------
 
class Pi0Vars {

//--------------------
// Instance Members --
//--------------------

public:

  enum {NPAIRS = 2}; // # of pairs we search for

  // Constructors
  Pi0Vars(const BtaCandidate & pi0, 
	  const HepAList<BtaCandidate> & pi0List,
	  const EventInfo & eventInfo,
	  const HepAList<BtaCandidate> * gammaList = 0,
	  const BtaMcAssoc * truthmap = 0);

  // Destructor
  virtual ~Pi0Vars( );


  // Accessors:
  float mom() const {return _mom;}
  float mass() const {return _mass;}
  float helic() const {return _helic;}
  int good() const {return _good;}

  float massH(int i) const {assert(i>-1 && i<2); return _massH[i];}
  float massS(int i) const {assert(i>-1 && i<2); return _massS[i];}
  float helicH(int i) const {assert(i>-1 && i<2); return _helicH[i];}
  float helicS(int i) const {assert(i>-1 && i<2); return _helicS[i];}
  float asymH(int i) const {assert(i>-1 && i<2); return _asymH[i];}
  float asymS(int i) const {assert(i>-1 && i<2); return _asymS[i];}
  float momH(int i) const {assert(i>-1 && i<2); return _momH[i];}
  float momS(int i) const {assert(i>-1 && i<2); return _momS[i];}
  int goodH(int i) const {assert(i>-1 && i<2); return _goodH[i];}
  int goodS(int i) const {assert(i>-1 && i<2); return _goodS[i];}

  float tMassH(int i) const {assert(i>-1 && i<2); return _tMassH[i];}
  float tMassS(int i) const {assert(i>-1 && i<2); return _tMassS[i];}
  float tHelicH(int i) const {assert(i>-1 && i<2); return _tHelicH[i];}
  float tHelicS(int i) const {assert(i>-1 && i<2); return _tHelicS[i];}
  int tGoodH(int i) const {assert(i>-1 && i<2); return _tGoodH[i];}
  int tGoodS(int i) const {assert(i>-1 && i<2); return _tGoodS[i];}


protected:
  // Helpers:
  void calc();
  void check(int i) const;

  const BtaCandidate * mcFromReco(const BtaCandidate * reco) const;
  const BtaCandidate * gammaRecoFromMc(const BtaCandidate * mc) const;

  double calcHelicityCM(HepLorentzVector daughter,    // by value
			const HepLorentzVector & mother) const;

  void bestPair(const BtaCandidate * include,
		const BtaCandidate * exclude1,
		const BtaCandidate * exclude2,
		BtaCandidate *& partner,
		float & bestMass,
		float & pairMom,
		float & helic,
		float & asym,
		int & good) const;

  void bestPair(const BtaCandidate * include,
		const BtaCandidate * exclude1,
		const BtaCandidate * exclude2,		
		float & bestMass,
		float & helic,
		int & good) const;

  const BtaCandidate * good2Body(const BtaCandidate * g1, 
				 const BtaCandidate * g2, 
				 const PdtEntry * entry) const;
  
  const BtaCandidate * good2Body(const BtaCandidate * mother, 
				 const PdtEntry * entry) const;
  
  double dtrMass(const BtaCandidate * mother) const;


private:
  // Data:  
  const HepAList<BtaCandidate> * _pi0List;
  const HepAList<BtaCandidate> * _gammaList;
  const BtaCandidate * _pi0;
  const BtaMcAssoc * _truthMap;
  const Hep3Vector _boostVector;
  const HepLorentzVector _cmFrame;

  float _mom;
  float _mass;
  float _helic;
  int _good;

  float _massH[NPAIRS];
  float _massS[NPAIRS];
  float _helicH[NPAIRS];
  float _helicS[NPAIRS];
  float _asymH[NPAIRS];
  float _asymS[NPAIRS];
  float _momH[NPAIRS];
  float _momS[NPAIRS];
  int _goodH[NPAIRS];
  int _goodS[NPAIRS];

  float _tMassH[NPAIRS];
  float _tMassS[NPAIRS];
  float _tHelicH[NPAIRS];
  float _tHelicS[NPAIRS];
  int _tGoodH[NPAIRS];
  int _tGoodS[NPAIRS];

};

#endif


