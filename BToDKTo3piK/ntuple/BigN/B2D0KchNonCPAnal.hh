//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: B2D0KchNonCPAnal.hh,v 1.3 2006/01/17 17:25:13 zhangjl Exp $
//
// Description:
//	Class B2D0KchNonCPAnal
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Qinglin Zeng                   Original Author
//   colorado state University
//     
// Copyright Information:
//	Copyright (C) 2003		
//
//------------------------------------------------------------------------

#ifndef B2D0KCHNONCPANAL_HH
#define B2D0KCHNONCPANAL_HH

#include <fstream>

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"
#include "CLHEP/Alist/AList.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "PDT/PdtLund.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"
#include "Framework/AbsParmBool.hh"
#include "Framework/AbsParmDouble.hh"
#include "Framework/AbsParmString.hh"
#include "Framework/AbsParmGeneral.hh"
#include "Beta/EventInfo.hh"
#include "BetaMiniUser/D0KchCand.hh"
#include "BetaCoreTools/BtaPrintTree.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class AbsBTagger;
class HepHistogram;
class HepTuple;
class BtaMcAssoc;
class BtaCandidate;
class Hep3Vector;

//		---------------------
// 		-- Class Interface --
//		---------------------
 
class B2D0KchNonCPAnal : public AppModule {

public:

  //----------------
  // Static members:
  //----------------

  static bool inRange(double val, double low, double high) {
    assert(low <= high);
    return (val > low && val < high) || high == low;
  }

  // check the daughers
  static int checkDaughters(const BtaCandidate * mother,
                             int num);

  // 2 body 
  static void findDaughters(const BtaCandidate * mother,
                            BtaCandidate *& dtr1, int type1,
                            BtaCandidate *& dtr2, int type2,
                            bool useAbs = true);

  // 3 body decay
  static void findDaughters(const BtaCandidate * mother,
			    BtaCandidate *& dtr1, int type1,
			    BtaCandidate *& dtr2, int type2,
			    BtaCandidate *& dtr3, int type2,
			    bool useAbs = true);

  static void findStableDescendents(const BtaCandidate * mother, 
                              HepAList<BtaCandidate> & descendents);
  
  static void findD0Descendents(const BtaCandidate * D0cand,
                                BtaCandidate *& descend1, int type1,
                                BtaCandidate *& descend2, int type2,
                                BtaCandidate *& descend3, int type3);

  static void findD0Descendents(const BtaCandidate * D0cand,
                                BtaCandidate *& descend1, int type1,
                                BtaCandidate *& descend2, int type2);

  static bool  findD0Descendents(const BtaCandidate * D0cand,
                  int type1,int type2,int type3, bool useAbs = false);

  static bool  findD0Descendents(const BtaCandidate * D0cand,
                  int type1,int type2, bool useAbs = false);

  static void findDecayBchToD0X(const BtaCandidate * Bcand,
				 BtaCandidate *& D0,
 				 BtaCandidate *& BDau,
                                 int type,bool useAbs = true );

  static void findDecayBchToDstar0X(const BtaCandidate * Bcand,
		  		BtaCandidate *& Dstar0,
		  		BtaCandidate *& BDau,
				int type,bool useAbs = true );
						   
  static void findDecayB0ToDstarcX(const BtaCandidate * Bcand,
		                BtaCandidate *& Dstarc,
				BtaCandidate *& BDau,			
 				int type,bool useAbs = true );

  static void getExclParams(const BtaCandidate * exclB,
			    double & exclMES, double & exclDeltaE,
			    double & exclD0Mass );	

  static double calcHelicity(HepLorentzVector daughter,
                      const HepLorentzVector & mother,
                      HepLorentzVector grandmother);
  // Constructors
  B2D0KchNonCPAnal( const char* const theName, const char* const theDescription );

  // Destructor
  virtual ~B2D0KchNonCPAnal( );

  // Operations
  virtual AppResult beginJob( AbsEvent* anEvent );
  virtual AppResult event( AbsEvent* anEvent );
  virtual AppResult endJob  ( AbsEvent* anEvent );
    
protected:

  // Helper functions
  bool isMC() const;
  bool isCategory(HepString cat, int bmod1, 
		  int bmod2, int dflg, int drec);

  BtaCandidate * ancestor(BtaCandidate * descent);
  
  void ntuple(const vector<D0KchCand> & d0KchCandidates);  

  void dumptrack(BtaCandidate* track);
  
  void dumpGamma(BtaCandidate* gamma);
  
  void CandVtxFit(const HepAList<BtaCandidate>& dauList, bool useBeamSpotConstrain=false, bool useMassConstraint=false, float candMass=0.0);
  
  int calNPhot(double sigLv1,double nExp);
   
  bool printCand(const char * particle) const;
  
  BtaPrintTree _printTree; 

private:

  //static data:
  static Hep3Vector _boostVector;

  // Exclusive analysis tcl flags and cuts:
  AbsParmIfdStrKey       _exclBchD0KchList;
  AbsParmIfdStrKey	 _btaGammaList;
  AbsParmIfdStrKey       _btaKsList;
  AbsParmIfdStrKey	 _btaPi0List;
  AbsParmIfdStrKey       _btaKsLsList;
  AbsParmIfdStrKey	 _btaTrkList;
  AbsParmIfdStrKey       _BtaTruthList;

  // Beam energy (= CM energy / 2) used in off-4s runs. Should be the
  // average on-4s beam energy for this run period:
  AbsParmDouble _offResEBeam;

  // The minimum difference between the CM energy and the Upsilon(4S)
  // mass that indicates that this is an off-resonance run (must be
  // positive):
  AbsParmDouble _offResEBeamDiff;

  AbsParmDouble _exclMESMin;      // low  cut on mES
  AbsParmDouble _exclMESMax;      // high cut on mES
  AbsParmDouble _exclDeltaEMin;   // low  cut on delta E - 0
  AbsParmDouble _exclDeltaEMax;   // high cut on delta E - 0 
  AbsParmDouble _exclD0MassMin;   // low  cut on D mass - nominal D0 mass
  AbsParmDouble _exclD0MassMax;   // high cut on D mass - nominal D0 mass
  AbsParmDouble _exclpi0MassMin;  // low  cut on pi0 mass
  AbsParmDouble _exclpi0MassMax;  // high cut on pi0 mass
  AbsParmDouble _exclKsMassMin;   // low  cut on Ks mass
  AbsParmDouble _exclKsMassMax;   // high cut on Ks mass
  AbsParmDouble _ListSelector; //List selector for D0 list

  AbsParmBool _fillNtp;
  AbsParmBool _dumpNtuple;
  AbsParmBool _debug;
  AbsParmString _printCands;
  AbsParmBool _treatAsMC;
  AbsParmString _evtCat;
  AbsParmString _outSeq;

  //ntuple
  int _numEventsRead;
  HepTuple* _tuple;

  //files
  std::ofstream _mES_Sig;
  std::ofstream _mES_Sideband;


  // Masses:
  double _massUpsilon;
  double _massBch;
  double _massD0;
  double _massPi0;
  double _massKs;
  double _massKch;
  double _massPi;
  double _massRho;
  double _massMuon;
  double _massElec;
 
  //bit map for particle Id
  double _muonId;
  double _elecId;
  double _kaonId;
  double _pionId;
  double _protId;
 
 // trk info
  double _pocad0;
  double _pocaz0;
  double _thC;
  double _thCErr;
  double _phi;
  double _sigoffprot;
  double _KaonNN;
  double _NPhot;
  double _NBkPhot;
  double _drcKaonCon;
  double _dchKaonCon;
  double _svtKaonCon;
  double _nSvtHits;
  double _nDchHits;
  double _pidHypo;

  // gamma info
  double _GamEn;
  double _latMom;
  double _bumRawE;

  //vertexed a candidate
  BtaCandidate _fittedCand;

  // Beam energy of the event (can be modified for off-4s runs):
  double _eBeam;
  BbrPointErr _beamSpot;

  // The CM B momentum given _eBeam:
  double _pB;

  //data share between functions
  BtaMcAssoc             * _truthMap;
  HepAList<BtaCandidate> * _mcList;

  AbsBTagger* _softPiTagger;
  AbsBTagger* _electronTagger;
  AbsBTagger* _muonTagger;
  AbsBTagger* _kinematicLeptonTagger; 

  AbsBTagger* _ElbTagger; //test purpose

  D0KchCand _d0KchCand;
};

#endif
