//-------------------------------------------------------------------------
// File and Version Information:
// 	$Id: B2D0KchNonCPAnal.cc,v 1.3 2006/01/17 17:25:12 zhangjl Exp $
//
// Description:
//	Class B2D0KchNonCPAnal
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Qinglin Zeng                 Original Author
//  colorado state university
//
// Copyright Information:
//	Copyright (C) 2003
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------
#include "BetaMiniUser/B2D0KchNonCPAnal.hh"

//-------------
// C Headers --
//-------------
#include <assert.h>

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cfloat>
//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/getTmpAList.hh"
#include "AbsEnv/AbsEnv.hh"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEventTag/AbsEventTag.hh"
#include "EventTagTools/TagBitModule.hh"
#include "AbsEvent/AbsEventID.hh"
//#include "AssocTools/AstKeyMap.hh"
#include "Beta/BtaCandidate.hh"
#include "Beta/BtaAbsVertex.hh"
#include "BetaCoreTools/BtaMcAssoc.hh"
#include "BetaEvent/BtaParam.hh"
#include "CLHEP/Alist/AList.h"
#include "CLHEP/Alist/ConstAList.h"
#include "CLHEP/Alist/ConstAIterator.h"
#include "CLHEP/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "GenEnv/GenEnv.hh"
#include "ErrLogger/ErrLog.hh"
#include "Framework/AbsParmBool.hh"
#include "Framework/AbsParmGeneral.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"
#include "HepTuple/TupleManager.h" 
#include "HepTuple/Tuple.h"
#include "HepTuple/Histogram.h"
#include "HepTuple/HTValOrderedVector.h"
#include "HepTuple/HTValVector.h"
#include "BetaMicroAdapter/BtaMicroAdapter.hh"
#include "BetaMicroAdapter/BtaCalQual.hh"
#include "BetaMicroAdapter/BtaPidQual.hh"
#include "BetaMicroAdapter/BtaTrkQual.hh"
#include "BetaMicroAdapter/BtaPidInfo.hh"
#include "ProxyDict/IfdStrKey.hh"
#include "ProxyDict/IfdKey.hh"
#include "PDT/PdtLund.hh"
#include "PDT/PdtEntry.hh"
#include "PDT/Pdt.hh"
#include "PepEnvData/PepBeams.hh"
#include "PepEnv/PepEnv.hh"
#include "AbsEnv/AbsEnv.hh"
#include "BaBar/Constants.hh"
#include "BetaMiniUser/D0KchCand.hh" 
#include "BetaMiniUser/Pi0Vars.hh"

#include "BetaCoreTools/BtaThrust.hh"
#include "BetaCoreTools/BtaFoxWolfMom.hh"
#include "BetaCoreTools/BtaBooster.hh"
#include "BetaCoreTools/BtaOpMakeTree.hh"
#include "BetaCoreTools/BtaBVariables.hh"
#include "BetaCoreTools/BtaHelicity.hh"

#include "BetaPid/PidPionLHSelector.hh"
#include "BetaPid/PidKaonMicroSelector.hh"
#include "BetaPid/PidKaonSMSSelector.hh"
//#include "BetaPid/PidLHElectronSelector.hh"
#include "BetaPid/PidElectronMicroSelector.hh"
#include "BetaPid/PidMuonLikeSelector.hh"
#include "BetaPid/PidMuonMicroSelector.hh"

#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkFit.hh"
#include "TrkBase/TrkExchangePar.hh"

#include "ProbTools/NumRecipes.hh"
#include "ProbTools/probab.hh"

#include "Q2BUser/Q2BCHBEventShape.hh"

#include "VtxFitter/VtxFitterOper.hh"
#include "VtxFitter/VtxTagBtaSelFit.hh"
#include "BbrGeom/BbrDoubleErr.hh"

#include "BTaggingTools/BtgStaticFactory.hh"
#include "BTaggingTools/BtgVariables.hh"
#include "AbsBTagging/AbsBTagger.hh"
#include "AbsBTagging/AbsTaggingFactory.hh"
//#include "AbsBTagging/AbsBTag.hh"

#include "TrajGeom/TrkLineTraj.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkAbsFit.hh"
#include "BbrGeom/Trajectory.hh"

#include "CompositionUtils/CompBaseAnalUtil.hh"
#include "EidData/EidEventTriplet.hh"

using namespace std;
//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------
//----------------------------------------------
// --Static Data & Function Member Definitions--
// ---------------------------------------------
//----------------------------------------------------------------------
Hep3Vector B2D0KchNonCPAnal::_boostVector;

//----------------------------------------------------------------------
//check the daughter whether it is a gamma that is due to radiation
int B2D0KchNonCPAnal::checkDaughters(const BtaCandidate * mother,int num) {
  if(0==mother) return 0;
  HepAListIterator<BtaCandidate> daugIter(mother->daughterIterator());
  daugIter.rewind();
  BtaCandidate * tmpDaug = 0;
  int totdau = 0;
  while(tmpDaug = daugIter()) {
    int tmpId = 0;
    tmpId = tmpDaug->pdtEntry()->lundId();
    if( tmpId == PdtLund::gamma ) {
      if(tmpDaug->energy() <= 0.02 ) continue;
    }
    totdau += 1;
  }
  if(totdau != num ) return 0;
  return 1;
}

//----------------------------------------------------------------------    
// Finds the 2 daughters of types type1 & type2 (PdtLund types) of mother,
// returns 0 for anything not found.
void B2D0KchNonCPAnal::findDaughters(const BtaCandidate * mother,
				     BtaCandidate *& dtr1, int type1,
				     BtaCandidate *& dtr2, int type2,
				     bool useAbs) {

  dtr1 = 0;
  dtr2 = 0;

  if (0 == mother){
    return;
  }
  
  if( checkDaughters(mother, 2) == 0 ) return; 
  HepAListIterator<BtaCandidate> daugIter(mother->daughterIterator());
  daugIter.rewind();
  BtaCandidate * tmpDaug = 0;
  int tmpID = 0;
  while(0 != (tmpDaug = daugIter())){
    tmpID = tmpDaug->pdtEntry()->lundId();
    if (useAbs) {
      tmpID = abs(tmpID);
    }

    if (tmpID == type1) {
      dtr1 = tmpDaug;
      tmpID = 0;
    }
    else if (tmpID == type2) {
      dtr2 = tmpDaug;
      tmpID = 0;
    }
  }
}

//-------------------------------------------------------------------
// Findse if (tmpID == typthe 3 daughters of types type1 & type2 
// (PdtLund types) of mother, returns 0 for anything not found.
void B2D0KchNonCPAnal::findDaughters(const BtaCandidate * mother,
				     BtaCandidate *& dtr1, int type1,
				     BtaCandidate *& dtr2, int type2,
				     BtaCandidate *& dtr3, int type3,
				     bool useAbs) {
  
  dtr1 = 0;
  dtr2 = 0;
  dtr3 = 0;

  if (0 == mother){
    return;
  }
  
  if(checkDaughters(mother,3)==0) return; 
  
  HepAListIterator<BtaCandidate> daugIter(mother->daughterIterator());
  daugIter.rewind();
  BtaCandidate * tmpDaug = 0;
  int tmpID = 0;
  while(tmpDaug = daugIter()){
    tmpID = tmpDaug->pdtEntry()->lundId();
    if (useAbs) {
      tmpID = abs(tmpID);
    }
    if (tmpID == type1) {
      dtr1 = tmpDaug;
    }
    else if (tmpID == type2) {
      dtr2 = tmpDaug;
    }
    else if (tmpID == type3) {
      dtr3 = tmpDaug;
    }	
  }
}

       
//------------------------------------------------------------------
// used to find the "stable" daughters of a compound particle as B meson, 
// D0, this does not work for D*0 that has a photon
void B2D0KchNonCPAnal::findStableDescendents(const BtaCandidate * mother,
				    HepAList<BtaCandidate> & descendents){
  static bool recursive = false;
  if(false == recursive) {
    while(descendents.length() > 0) {
      static const unsigned int zero =0;
      descendents.remove(zero);
    }
  }		  

  if( (0 == mother) ) {	return;  }
  
  HepAListIterator<BtaCandidate> daugIter(mother->daughterIterator());
  daugIter.rewind();
  BtaCandidate * tmpDaug = 0;
  while (tmpDaug = daugIter()) {
    int tempID = 0;
    tempID = tmpDaug->pdtEntry()->lundId();
    // if(useAbs) {
    //	 tempID = abs(tempID);
    // }

    if( abs(tempID) == PdtLund::gamma ) {
      if(abs(mother->pdtEntry()->lundId()) == PdtLund::D_star0
	 && tmpDaug->energy()>0.015 ) {
	descendents.append(tmpDaug);
      }
    }

    if((abs(tempID) == PdtLund::pi_plus) || 
       (abs(tempID) == PdtLund::K_plus)) {
      //found a stable chargeed kaon or pion
      descendents.append(tmpDaug);
    } else if( tmpDaug->nDaughters() > 0 ){
      if((tempID == PdtLund::K_S0) || (tempID == PdtLund::pi0)) {
	// found a intermid state particle Pi0 or K_S0;
	descendents.append(tmpDaug);
      } else {
	recursive = true;
	findStableDescendents(tmpDaug,descendents);
	recursive = false;
      }
    }     
  }
}  

//-----------------------------------------------------------------------
//used for finding D0 Descendents: the number is 3
void B2D0KchNonCPAnal::findD0Descendents( const BtaCandidate * D0,
			BtaCandidate *& descend1, int type1,
                        BtaCandidate *& descend2, int type2,
                        BtaCandidate *& descend3, int type3) {
  if(0 == D0 ) { return; }
  // only look for D0->K*K, RhoPi and KKPi0,KPiKs,PiPiPi0
  HepAList<BtaCandidate>  D0decayList;
  findStableDescendents( D0, D0decayList);
  //shoud set false
  if( D0decayList.length()>3 || D0decayList.length() == 0) {
    return; 
    //false value
  }
  
  HepAListIterator<BtaCandidate> D0descendIter(D0decayList);
  D0descendIter.rewind();	
  BtaCandidate *ListItem = 0;
  int dauID = 0;
  while( ListItem = D0descendIter() ) {
    dauID = ListItem->pdtEntry()->lundId();
    if(dauID == type1) { descend1 = ListItem; }
    if(dauID == type2) { descend2 = ListItem; }
    if(dauID == type3) { descend3 = ListItem; }
  }
}

//---------------------------------------------------------------
//used for finding D0 Descendents: the number is 2
void B2D0KchNonCPAnal::findD0Descendents( const BtaCandidate * D0,
                        BtaCandidate *& descend1, int type1,
                        BtaCandidate *& descend2, int type2) {
  if(0 == D0 ) { return; }
  
  HepAList<BtaCandidate>  D0decayList;
  findStableDescendents( D0, D0decayList);
  //shoud set false
  if( D0decayList.length()>2 || D0decayList.length() == 0) {
    return; 
    //false value
  }
  
  HepAListIterator<BtaCandidate> D0descendIter(D0decayList);
  D0descendIter.rewind();	
  BtaCandidate *ListItem = 0;
  int dauID = 0;
  while( ListItem = D0descendIter() ) {
    dauID = ListItem->pdtEntry()->lundId();
    if(dauID == type1) { descend1 = ListItem; }
    if(dauID == type2) { descend2 = ListItem; }
  }
}

//-----------------------------------------------------------------------
bool B2D0KchNonCPAnal::findD0Descendents( const BtaCandidate * D0,
                        int type1,int type2,int type3, bool useAbs) {
  if(0 == D0 ) { return false; }
  HepAList<BtaCandidate>  D0decayList;
  findStableDescendents( D0, D0decayList);
  //shoud set false
  if( D0decayList.length()>3 || D0decayList.length() == 0) {
    return false; 
    //false value
  }
  
  HepAListIterator<BtaCandidate> D0descendIter(D0decayList);
  D0descendIter.rewind();	
  BtaCandidate* ListItem(0);
  int dauID = 0;
  BtaCandidate* descend1(0);
  BtaCandidate* descend2(0);
  BtaCandidate* descend3(0);
  while( ListItem = D0descendIter() ) {
    dauID = ListItem->pdtEntry()->lundId();
    if(useAbs) dauID = abs(dauID);
    if(dauID == type1) { descend1 = ListItem; }
    if(dauID == type2) { descend2 = ListItem; }
    if(dauID == type3) { descend3 = ListItem; } 
  }
  //if(descend1 != 0) cout<<" 1. "<<descend1->pdtEntry()->name();
  //if(descend2 != 0) cout<<" 2. "<<descend2->pdtEntry()->name();
  //if(descend3 != 0) cout<<" 3. "<<descend3->pdtEntry()->name();
  if(0 != descend1 && 0 != descend2 && 0 != descend3 ){
    return true;
  } else {
    return false;
  }
}

//-------------------------------------------------------------------
bool B2D0KchNonCPAnal::findD0Descendents( const BtaCandidate * D0,
                        int type1,int type2, bool useAbs) {
  if(0 == D0 ) { return false; }
  HepAList<BtaCandidate>  D0decayList;
  findStableDescendents( D0, D0decayList);
  //shoud set false
  if( D0decayList.length()>2 || D0decayList.length() == 0) {
    return false; 
    //false value
  }
  
  HepAListIterator<BtaCandidate> D0descendIter(D0decayList);
  D0descendIter.rewind();	
  BtaCandidate* ListItem(0);
  int dauID = 0;
  BtaCandidate* descend1(0);
  BtaCandidate* descend2(0);
  while( ListItem = D0descendIter() ) {
    dauID = ListItem->pdtEntry()->lundId();
    if(useAbs) dauID = abs(dauID);
    if(dauID == type1) { descend1 = ListItem; }
    if(dauID == type2) { descend2 = ListItem; }
  }
  //if(descend1 != 0) cout<<" 1. "<<descend1->pdtEntry()->name();
  //if(descend2 != 0) cout<<" 2. "<<descend2->pdtEntry()->name();
  if(0 != descend1 && 0 != descend2){
    return true;
  } else {
    return false;
  }
}

//------------------------------------------------------------------------
// B+ -> D0 X
void B2D0KchNonCPAnal::findDecayBchToD0X(const BtaCandidate * exclB,
                         BtaCandidate *& D0,
                         BtaCandidate *& BDau,
                         int type, bool useAbs) {
  // The functions called by this function should set the appropriate
  // values even if exclB is 0.
  
  // Get Bch daughers D0 and Kch:
  findDaughters(exclB,D0,PdtLund::D0,BDau,type, useAbs);
}

//------------------------------------------------------------------------
// B+ -> D*0 X
void B2D0KchNonCPAnal::findDecayBchToDstar0X(const BtaCandidate * exclB,
                         BtaCandidate *& Dstar0,
			 BtaCandidate *& BDau,
			 int type, bool useAbs) {
  // The functions called by this function should set the appropriate
  // values even if exclB is 0.
  // Get Bch daughers D0 and Kch:
  findDaughters(exclB,Dstar0,PdtLund::D_star0,BDau,type, useAbs);
}   
 
//------------------------------------------------------------------------
// B0 -> D*c X
void B2D0KchNonCPAnal::findDecayB0ToDstarcX(const BtaCandidate * exclB,
                         BtaCandidate *& Dstarc,
			 BtaCandidate *& BDau,
			 int type, bool useAbs) {
  // The functions called by this function should set the appropriate
  // values even if exclB is 0.
  // Get Bch daughers D0 and Kch:
  findDaughters(exclB,Dstarc,PdtLund::D_star_plus,BDau,type, useAbs);
}   

//-------------------------------------------------------------------------
double B2D0KchNonCPAnal::calcHelicity(HepLorentzVector daughter, 
                               const HepLorentzVector & mother,
                               HepLorentzVector grandmother) {

  // Transform the grandmother and the daughter to the mother's frame:
  const Hep3Vector boostVector = -(mother.boostVector());
  daughter.boost(boostVector);
  grandmother.boost(boostVector);
  return cos(grandmother.angle(daughter));
}


//-------------------------------------------------------------------
// These ntupleVar*** functions cannot be member functions, because then
// the linux compiler can't find a match for them. The reason is probably 
// that the compiler doesn't know how to handle pointers to class member 
// functions in the argument list of a member function of another class.
// But it does know how to handle it for a global function:
//-------------------------------------------------------------------

static const char * EVENT_BLOCK_NAME = "EVENT";    
static const char * CANDS_BLOCK_NAME = "D0Kch";
static const char * CANDS_INDEX_NAME = "Length";   

//-------------------------------------------------------------------
static void ntupleVarInt(const vector<D0KchCand> & 
			 d0KchCandidates,
			 const char * name,
			 int (D0KchCand::*function)() const,
			 HepTuple * tuple) {
  
  const int nCands = d0KchCandidates.size();
  HTValOrderedVector<int> theVector(nCands);
  
  // Fill the vector:
  for (int i = 0; i < nCands; i++){
    theVector.append((d0KchCandidates[i].*function)());
  }
  
  // Dump to the ntuple:
  static const int zero = 0;
  tuple->column(name, theVector, CANDS_INDEX_NAME, zero, CANDS_BLOCK_NAME);
}

//-------------------------------------------------------------------
static void ntupleVar(const vector<D0KchCand> & 
		      d0KchCandidates,
		      const char * name,
		      float (D0KchCand::*function)() const,
		      HepTuple * tuple) {
  
  const int nCands = d0KchCandidates.size();
  HTValOrderedVector<float> theVector(nCands);
  
  // Fill the vector:
  for (int i = 0; i < nCands; i++){
    theVector.append((d0KchCandidates[i].*function)());
  }
  
  // Dump to the ntuple:
  static const float zero = 0;
  tuple->column(name, theVector, CANDS_INDEX_NAME, zero, CANDS_BLOCK_NAME);
}

//--------------------------------------------------------------------
static void ntupleVarIndexInt(const vector<D0KchCand> &
			      d0KchCandidates,
			      const char * name,
			      int index,
			      int (D0KchCand::*function)(int) const,
			      HepTuple * tuple) {

  HepString nameIndex = HepString(name) + HepString(index);
  
  const int nCands = d0KchCandidates.size();
  HTValOrderedVector<int> theVector(nCands);

  // Fill the vector:
  for (int i = 0; i < nCands; i++){
    theVector.append((d0KchCandidates[i].*function)(index));
  }

  // Dump to the ntuple":
  static const int zero = 0;
 tuple->column(nameIndex, theVector, CANDS_INDEX_NAME, zero, 
	       CANDS_BLOCK_NAME);
}

//--------------------------------------------------------------------------
static void ntupleVarIndex(const vector<D0KchCand> &
			   d0KchCandidates,
			   const char * name,
			   int index,
			   float (D0KchCand::*function)(int) const,
			   HepTuple * tuple) {

  HepString nameIndex = HepString(name) + HepString(index);

  const int nCands = d0KchCandidates.size();
  HTValOrderedVector<float> theVector(nCands);

  // Fill the vector:
  for (int i = 0; i < nCands; i++){
    theVector.append((d0KchCandidates[i].*function)(index));
  }

  // Dump to the ntuple":
  static const float zero = 0;
  tuple->column(nameIndex, theVector, CANDS_INDEX_NAME, zero, 
		CANDS_BLOCK_NAME);
}

//----------------------------------------
//-- Public Function Member Definitions --
//----------------------------------------

//----------------
// Constructors --
//----------------

// in general, a module constructor should not do much.  The begin(job) or
// begin(run) members are better places to put initialization
B2D0KchNonCPAnal::B2D0KchNonCPAnal( const char* const theName, 
		  const char* const theDescription )
  : AppModule( theName, theDescription )
  , _exclBchD0KchList("exclBchD0KchList", this, "BchToD0KchNonCPHard") 
  , _btaGammaList("gammalist",this,"GoodPhotonLoose")
  , _btaKsList("KsList", this, "KsDefault")
  , _btaPi0List("pi0List", this, "pi0AllLoose")
  , _btaKsLsList("KsLsList", this, "KsDefault")	
  , _btaTrkList("btaTrkList",this,"GoodTracksVeryLoose")
  , _BtaTruthList("btaTruthList",this,"MCTruth") 
  , _offResEBeam("offResEBeam", this, 5.2891 )
  , _offResEBeamDiff("offResEBeamDiff", this, 0.030)
  , _exclMESMin("exclMESMin", this, 5.2)
  , _exclMESMax("exclMESMax", this, 5.3)
  , _exclDeltaEMin("exclDeltaEMin", this, -0.3)
  , _exclDeltaEMax("exclDeltaEMax", this,  0.3)
  , _exclD0MassMin("exclD0MassMin", this, -0.090)
  , _exclD0MassMax("exclD0MassMax", this,  0.090)
  , _exclpi0MassMin("exclpi0MassMin", this, -0.025)
  , _exclpi0MassMax("exclpi0MassMax", this,  0.025)
  , _exclKsMassMin("exclKsMassMin", this, -0.02)
  , _exclKsMassMax("exclKsMassMax", this,  0.02)
  , _ListSelector("listSelector", this, 1 )
  , _fillNtp("fillNtple",this,false)
  , _dumpNtuple("dumpNtuple", this,true) 
  , _debug("debug", this, false)
  , _printCands("printCands",this,"") 
  , _treatAsMC("treatAsMC", this, false)
  , _evtCat("evtCat", this, "All")
  , _outSeq("outSeq", this, "Nothing")
  , _printTree(BtaPrintTree::printName, BtaPrintTree::printP4)
{
  commands()->append( &_exclBchD0KchList);
  commands()->append( &_btaGammaList);
  commands()->append( &_btaKsList);
  commands()->append( &_btaPi0List);
  commands()->append( &_btaKsLsList);
  commands()->append( &_btaTrkList);
  commands()->append( &_BtaTruthList);
  commands()->append( &_offResEBeam);
  commands()->append( &_offResEBeamDiff);
  commands()->append( &_exclMESMin);
  commands()->append( &_exclMESMax);
  commands()->append( &_exclDeltaEMin);
  commands()->append( &_exclDeltaEMax);
  commands()->append( &_exclD0MassMin);
  commands()->append( &_exclD0MassMax);
  commands()->append( &_exclpi0MassMin);
  commands()->append( &_exclpi0MassMax);
  commands()->append( &_exclKsMassMin);
  commands()->append( &_exclKsMassMax);
  commands()->append( &_ListSelector);
  commands()->append( &_fillNtp);
  commands()->append( &_dumpNtuple);
  commands()->append( &_debug);
  commands()->append( &_printCands);
  commands()->append( &_treatAsMC);
  commands()->append( &_evtCat);
  commands()->append( &_outSeq);
}

//--------------
// Destructor --
//--------------

// The destructor should be limited to undoing the work of the constructor
B2D0KchNonCPAnal::~B2D0KchNonCPAnal() { } ; 

//--------------
// Operations --
//--------------

//member function

// The begin(AppJob*) member function is run before any events are
// processed.  In sample code, it opens the output histogram file
// and then books a histogram and optionally an ntuple.
AppResult
B2D0KchNonCPAnal::beginJob( AbsEvent* anEvent )
{
  ErrMsg(routine)<<"Begin job"<<endmsg; 

  //set masses
  _massUpsilon = Pdt::lookup(PdtLund::Upsilon_4S)->mass();
  _massBch = Pdt::lookup(PdtLund::B_plus)->mass();
  _massD0  = Pdt::lookup(PdtLund::D0)->mass();
  _massPi0 = Pdt::lookup(PdtLund::pi0)->mass();
  _massKs = Pdt::lookup(PdtLund::K_S0)->mass();  
  _massKch = Pdt::lookup(PdtLund::K_plus)->mass(); 
  _massPi = Pdt::lookup(PdtLund::pi_plus)->mass();
  _massMuon = Pdt::lookup(PdtLund::mu_minus)->mass();
  _massElec = Pdt::lookup(PdtLund::e_minus)->mass();
	
  _softPiTagger = BtgStaticFactory::factoryInstance()->
          generateTagger("BtgSlowPionTag");
  _softPiTagger->initialize();

  _electronTagger = BtgStaticFactory::factoryInstance()->
          generateTagger("BtgElectronTag");
  _electronTagger->initialize();

  _muonTagger = BtgStaticFactory::factoryInstance()->
          generateTagger("BtgMuonTag");
  _muonTagger->initialize();

  _kinematicLeptonTagger = BtgStaticFactory::factoryInstance()->
          generateTagger("BtgKinematicLeptonTag");
  _kinematicLeptonTagger->initialize();

  HepTupleManager* manager = gblEnv->getGen()->ntupleManager();
  assert(manager != 0);

  // book an n-tuple if enabled via tcl parameter
  if( _fillNtp.value() )  _tuple = manager->ntuple(name());

  //output files for truth trees
  if (_evtCat.value() != "All") {
    HepString sigFileName="/nfs/farm/babar/AWG17/BCK/Log/"+
      (HepString)(_outSeq.value())+"_"+
      (HepString)(_evtCat.value())+"_mES_Sig.log";
    HepString sbFileName="/nfs/farm/babar/AWG17/BCK/Log/"+
      (HepString)(_outSeq.value())+"_"+
      (HepString)(_evtCat.value())+"_mES_Sideband.log";
    _mES_Sig.open((const char*)sigFileName,
		  ios_base::out | ios_base::app);
    _mES_Sideband.open((const char*)sbFileName , 
		       ios_base::out | ios_base::app);
  }

  return AppResult::OK;
}

// end(AppJob*) function is called after all events have been processed.
AppResult
B2D0KchNonCPAnal::endJob( AbsEvent* anEvent ) 
{

  if (_evtCat.value() != "All") {
    _mES_Sig.close();
    _mES_Sideband.close();
  }

  ErrMsg(routine) << "End job" << endmsg;
  return AppResult::OK;
}



// event function is called once per event
AppResult B2D0KchNonCPAnal::event( AbsEvent* anEvent )
{ 
   
  // Count events read:
  ++_numEventsRead;

  //Get Event Info
  IfdStrKey key("Default");
  HepAList< EventInfo >* infoList=NULL;
  getTmpAList(anEvent, infoList, key);
  EventInfo* eventInfo = infoList->first();

  // store event ID and time stamp
  const AbsEventID * eventID = Ifd<AbsEventID>::get(anEvent, "AbsEventID") ;
  EidEventTriplet theTriplet = eventID->eventIdTriplet() ;

  // gamma list  
  HepAList<BtaCandidate> * gammaList;
  getTmpAList(anEvent, gammaList, _btaGammaList.value());

  // pi0 List
  HepAList<BtaCandidate> * pi0List;
  getTmpAList(anEvent, pi0List, _btaPi0List.value());

  //trk list
  HepAList<BtaCandidate> * btaTrkList;
  getTmpAList(anEvent, btaTrkList, _btaTrkList.value());

  IfdStrKey key1(" MircoTaggingList ");
  HepAList<BtaCandidate>  * MicroTaggingList;
  getTmpAList(anEvent, MicroTaggingList, key1);
  
  // Get exclusive Bch ->D0 Kch reconstructed lists:
  HepAList<BtaCandidate> * BchD0KchList;
  getTmpAList (anEvent, BchD0KchList, _exclBchD0KchList.value());  

  // mc List  
  _truthMap = Ifd<BtaMcAssoc>::get(anEvent,"GHit"); 
  getTmpAList(anEvent,_mcList,_BtaTruthList.value() );

  // Get the tag info for the event:
  AbsEventTag * tag = Ifd<AbsEventTag>::get( anEvent );
  if (0 == tag){
    ErrMsg(routine) << name() << ": tag not found. Returning." << endmsg;
  }

  // get beam spot and boost vector:
  const PepBeams * pepBeams = gblEnv->getPep()->pepBeams();
  _boostVector = -(eventInfo->cmFrame().boostVector());
  _beamSpot = eventInfo->beamSpotBFlight();

  // Set the expected beam energy and B momentum for this event:
  _eBeam = pepBeams->energyCM() / 2.0;
  bool usingOffResEBeam = false;
  if (_massUpsilon - pepBeams->energyCM() > _offResEBeamDiff.value()){
    // Then this is an off-resonance run. Set the effective eBeam to what
    // the on-resonance runs have been for the run period:
    _eBeam = _offResEBeam.value();
    usingOffResEBeam = true;
  }
	  
  //make total list
  HepAList<BtaCandidate> allList(*btaTrkList);
  allList += *gammaList;
	  
  //event shape variables R2 cut
  BtaFoxWolfMom myFoxWolNeu ( &allList, eventInfo, 2);
  float foxWolfR2 = (float)myFoxWolNeu.FWnorm(2);  
  
  //event thrust
  BtaThrust thrusterNeu(allList,*eventInfo,BtaThrust::BTAllParticles);
  float thrust = thrusterNeu.thrust();
  
  BtaCandidate* AllListTmpCand(0);
  HepAListIterator<BtaCandidate> AllListTmpIter(allList);
  HepLorentzVector AllP4Hem1(0,0,0,0);	
  HepLorentzVector AllP4Hem2(0,0,0,0);	
  while( 0 != (AllListTmpCand = AllListTmpIter()) ) {
    HepLorentzVector tmp4Vect = AllListTmpCand->p4();
    tmp4Vect.boost(_boostVector);
    double costmp = cos(thrusterNeu.thrust_axis().angle(tmp4Vect.vect()));
    if( costmp >0.0 ) AllP4Hem1 += tmp4Vect;
    if( costmp <0.0 ) AllP4Hem2 += tmp4Vect;
  }
  double hemMass1 = AllP4Hem1.mag();
  double hemMass2 = AllP4Hem2.mag(); 
  if(hemMass2 > hemMass1 ) {
    double tmpVal = hemMass1;
    hemMass1 = hemMass2;
    hemMass2 = tmpVal;
  }	
	 
  /* Monte Carlo Part: used to find the decay mode from monte carlo truth
     and at the end it will be used to check the reconstructed particle 
     mc truth and then can be used to find the background. */

  int BdecModMC[] = { 0, 0};
  int D0CADecMode[] = { 0, 0 };
  double bx[] = {0.,0.};
  double by[] = {0.0,0.0};
  double bz[] = {0., 0.};
  double mcbFitChisq[] = {-100.,-100.};	  
  double mcBchrge(-10.);	
  int numBFound = 0;

  //B+ -> D0 K+  10 signal
  // all others are background
  int d0XType = 0;
  int D0CADecay = 0; 
  int bBaseType = -1;

  double Kch2MCMass = 0;
  double Kch2MCHelic = -999;
  double KpPi0MCMass = 0;
  double KpPi0MCHelic = -999;
  double KmPi0MCMass = 0;
  double KmPi0MCHelic = -999;
  double KpPi2MCMass = 0;
  double KpPi2MCHelic = -999;
  double KpKsMCMass = 0;
  double KpKsMCHelic =-999;
  double PimKsMCMass = 0;
  double PimKsMCHelic =-999;
  double KmPi2MCMass =  0;
  double KmPi2MCHelic =-999;
  double KmKsMCMass =  0;
  double KmKsMCHelic =-999;
  double PipKsMCMass = 0;
  double PipKsMCHelic =-999;
  double PiPiMCMass =  0;
  double PiPiMCHelic =-999;
  double PipPi0MCMass = 0;
  double PipPi0MCHelic =-999;
  double PimPi0MCMass =  0;
  double PimPi0MCHelic = -999;

  HepAListIterator<BtaCandidate> bIter(*_mcList);
  bIter.rewind();
  BtaCandidate * bCand = 0;
  while (bCand = bIter()) {
    d0XType = 0;
    bBaseType = -1;
    D0CADecay = 0;		
    const PdtLund::LundType bID = bCand->pdtEntry()->lundId();
    if(PdtLund::B0 == bID || PdtLund::anti_B0 == bID ||
       PdtLund::B_plus == bID || PdtLund::B_minus == bID) {
      ++numBFound;
      // get the Vtx info from the B
      if(numBFound==1) mcBchrge = bCand->charge();	
      mcbFitChisq[numBFound-1] = bCand->decayVtx()->chiSquared();
      bx[numBFound-1] = bCand->decayVtx()->point().x();
      by[numBFound-1] = bCand->decayVtx()->point().y();
      bz[numBFound-1] = bCand->decayVtx()->point().z();
      
      // MC truth of the bachelor kaon
      HepLorentzVector KchMCP4;
      // detect the B decay base type on the decay information:
      BtaCandidate * dCand = 0; 
      // B+ -> D0 X( K+ Rho+ pi+ K*+)
      // B+ -> D0 K+ signal 
      if( bBaseType == -1 ) {
	dCand = 0;
	BtaCandidate * KchcCand = 0;
	findDaughters(bCand,dCand,PdtLund::D0, 
		      KchcCand,PdtLund::K_plus,true );
	if((dCand != 0)&&(KchcCand !=0)) {//B+ -> D0 K+
	  bBaseType = 10;
	  // used for compond mass
	  KchMCP4 = KchcCand->p4();
	}
      }

      // B+ -> D0 pi+
      if(bBaseType == -1 ) { 
	dCand = 0; 	
	BtaCandidate * picCand = 0;
	findDaughters(bCand,dCand,PdtLund::D0,
		      picCand,PdtLund::pi_plus,true );
	if((dCand != 0)&&(picCand !=0)) { //B+ -> D0 pi+
	  bBaseType = 150;
	}      
      }
      
      // B+ -> D0 Rho+
      if( bBaseType == -1 ) {	
	dCand = 0;
	BtaCandidate * RhocCand = 0;
	findDaughters(bCand,dCand,PdtLund::D0,
		      RhocCand,PdtLund::rho_plus,true );
	if((dCand != 0)&&(RhocCand !=0)) { //B+ -> D0 Rho+
	  bBaseType = 160;
	}
      }
      
      // B+ -> D0 K*+
      if(bBaseType == -1 ) {
	dCand = 0;
	BtaCandidate * KstarcCand = 0;
	findDaughters(bCand,dCand,PdtLund::D0,
		      KstarcCand,PdtLund::K_star_plus,true); 		
	if((dCand != 0)&&(KstarcCand !=0)) { //B+ -> D0 K*+ (K+ pi0 only)
	  BtaCandidate * KstarDauCand1 = 0;
	  BtaCandidate * KstarDauCand2 = 0;
	  BtaCandidate * KstarDauCand3 = 0;
	  BtaCandidate * KstarDauCand4 = 0;
	  findDaughters(KstarcCand,KstarDauCand1,
			PdtLund::K_plus,KstarDauCand2,PdtLund::pi0,true);
	  findDaughters(KstarcCand,KstarDauCand3,
			PdtLund::K_S0,KstarDauCand4,PdtLund::pi_plus,true);
	  if((KstarDauCand1!=0&&KstarDauCand2!=0) || 
	     (KstarDauCand3!=0&&KstarDauCand4!=0)) { 
	    bBaseType = 110;
	  }
	}
      }
		
      // B+ -> D*0 K+
      if( bBaseType == -1 ) {
	BtaCandidate * dstar0Cand = 0;
	BtaCandidate * KchcCand = 0;
	findDaughters(bCand,dstar0Cand,PdtLund::D_star0,
		      KchcCand,PdtLund::K_plus, true);	
	if((dstar0Cand!=0)&&(KchcCand!=0)) { //B+ ->D*0 K+
	  dCand = 0;	
	  BtaCandidate * dstarLDauCand0 = 0;
	  BtaCandidate * dstarLDauCand1 = 0;
	  findDaughters(dstar0Cand,dCand,PdtLund::D0,
			dstarLDauCand0,PdtLund::pi0,true);
	  findDaughters(dstar0Cand,dCand,PdtLund::D0,
			dstarLDauCand1,PdtLund::gamma,true);
	  if(dCand!=0&&(dstarLDauCand0 != 0||dstarLDauCand1 != 0)) {
	    bBaseType = 130;
	  }
	}
      }
	     
      //B+ -> D*0 pi+
      if( bBaseType == -1 ) {
	BtaCandidate * dstar0Cand = 0;
	BtaCandidate * picCand = 0;  	    
	findDaughters(bCand,dstar0Cand,PdtLund::D_star0,
		      picCand,PdtLund::pi_plus,true);
	if((dstar0Cand!=0)&&(picCand!=0)) { //B+ ->D*0 pi+
	  dCand = 0;
	  BtaCandidate * dstarLDauCand0 = 0;
	  BtaCandidate * dstarLDauCand1 = 0;
	  findDaughters(dstar0Cand,dCand,PdtLund::D0,
			dstarLDauCand0,PdtLund::pi0,true);
	  findDaughters(dstar0Cand,dCand,PdtLund::D0,
			dstarLDauCand1,PdtLund::gamma,true);
	  if(dCand!=0&&(dstarLDauCand0 !=0 || dstarLDauCand1 != 0)) {
	    bBaseType = 170;
	  }
	}
      }

      // B+ -> D*0 K*+
      if( bBaseType == -1 ) {
	BtaCandidate * dstar0Cand = 0;
	BtaCandidate * KstarcCand = 0;
	findDaughters(bCand,dstar0Cand,PdtLund::D0,
		      KstarcCand,PdtLund::K_star_plus,true);
	if((dstar0Cand!=0)&&(KstarcCand!=0)) { //B+ ->D*0 K*+
	  dCand = 0;
	  BtaCandidate * dstarLDauCand = 0;	
	  BtaCandidate * dstarLDauCand1= 0;	
	  //D*0 -> D0 pi0	
	  findDaughters(dstar0Cand,dCand,PdtLund::D0,
			dstarLDauCand,PdtLund::pi0,true);
	  findDaughters(dstar0Cand,dCand,PdtLund::D0,
			dstarLDauCand1,PdtLund::gamma,true);
	  if(dCand!=0&&(dstarLDauCand != 0 ||dstarLDauCand1 != 0)) {
	    //K*+ -> K+ pi0/K*+ -> K_S0 pi+
	    BtaCandidate * KstarDauCand1 = 0;
	    BtaCandidate * KstarDauCand2 = 0;
	    BtaCandidate * KstarDauCand3 = 0;
	    BtaCandidate * KstarDauCand4 = 0;		
	    findDaughters(KstarcCand,KstarDauCand1,
			  PdtLund::K_plus,KstarDauCand2,PdtLund::pi0,true);
	    findDaughters(KstarcCand,KstarDauCand3,PdtLund::K_S0,
			  KstarDauCand4,PdtLund::pi_plus,true);
	    if((KstarDauCand1!=0&&KstarDauCand2!=0) || 
	       (KstarDauCand3!=0&&KstarDauCand4!=0)) bBaseType=140;
	  }
	}	
      }
      
      // B0 -> D*- K+
      if( bBaseType == -1 ) {
	BtaCandidate * dstarcCand = 0;
	BtaCandidate * KchcCand = 0;
	findDaughters(bCand,dstarcCand,PdtLund::D_star_plus,
		      KchcCand,PdtLund::K_plus,true);	 
	if((dstarcCand!=0)&&(KchcCand!=0)) { //B0 -> D*-(D0 pi-) K+
	  dCand = 0;
	  BtaCandidate * dstarLDauCand = 0;		
	  findDaughters(dstarcCand,dCand,PdtLund::D0,
			dstarLDauCand,PdtLund::pi_plus,true);
	  if(dCand!=0&&dstarLDauCand!=0) bBaseType = 120;
	}
      }

      // B0 -> D*- pi+
      if( bBaseType == -1 ) {
	BtaCandidate * dstarcCand = 0;
	BtaCandidate * picCand = 0;
	findDaughters(bCand,dstarcCand,PdtLund::D_star_plus,
		      picCand,PdtLund::pi_plus,true);
	if((dstarcCand!=0)&&(picCand!=0)) { //B0 -> D*-(D0 pi-) pi+
	  dCand = 0;
	  BtaCandidate * dstarLDauCand = 0;
	  findDaughters(dstarcCand,dCand,PdtLund::D0,
			dstarLDauCand,PdtLund::pi_plus,true);
	  if(dCand!=0&&dstarLDauCand!=0) bBaseType = 100;
	}
      }

      if( bBaseType == -1 ) {
	BtaCandidate * dstarcCand = 0;
	BtaCandidate * RhocCand = 0;
	findDaughters(bCand,dstarcCand,PdtLund::D_star_plus,
		      RhocCand,PdtLund::rho_plus,true);
	if((dstarcCand!=0)&&(RhocCand!=0)) { //B0 -> D*-(D0 pi-) rho+
	  dCand = 0;
	  BtaCandidate * dstarLDauCand = 0;
	  findDaughters(dstarcCand,dCand,PdtLund::D0,
			dstarLDauCand,PdtLund::pi_plus,true);
	  if(dCand!=0&&dstarLDauCand!=0) bBaseType = 90;
	}
      }

      //B0 -> D- rho+
      if( bBaseType == -1 ) {
	BtaCandidate * dcCand = 0;
	BtaCandidate * RhocCand = 0;
	findDaughters(bCand,dcCand,PdtLund::D_plus,
		      RhocCand,PdtLund::rho_plus,true);
	if(dcCand!=0&&RhocCand!=0) bBaseType = 80;
      }
      
      //B0 -> D- K*+
      if( bBaseType == -1 ) {
	BtaCandidate * dcCand = 0;
	BtaCandidate * KstCand = 0;
	findDaughters(bCand,dcCand,PdtLund::D_plus, 
		      KstCand, PdtLund::K_star_plus,true);
	if(dcCand!=0&&KstCand!=0) bBaseType = 70;
      }
 
      if( bBaseType == -1 ) {
	BtaCandidate * dcCand = 0;
	BtaCandidate * kCand = 0;
	findDaughters(bCand,dcCand,PdtLund::D_plus, kCand, 
		      PdtLund::K_plus,true);
	if(dcCand!=0&&kCand!=0) bBaseType = 60;
      }
      
      if( bBaseType == -1 ) {
	BtaCandidate * dcCand = 0;
	BtaCandidate * piCand = 0;
	findDaughters(bCand,dcCand,PdtLund::D_plus, piCand, 
		      PdtLund::pi_plus,true);
	if(dcCand!=0&&piCand!=0) bBaseType = 50;
      }
 
      // ========for B+ -> D0(pi+ pi- pi0) K+ mode======  
      //Charmless decay:
      //B+->K*+(K+ pi0) pi+ pi-
      if( bBaseType == -1) {
	BtaCandidate * KstarCh = 0;
	BtaCandidate * PiChp = 0;
	BtaCandidate * PiChm = 0;    
	findDaughters(bCand,PiChp,PdtLund::pi_plus,
		      PiChm,PdtLund::pi_minus,KstarCh,
		      PdtLund::K_star_plus,false);
	if((KstarCh==0)&&(PiChp==0)||(PiChm==0)) {//B- -> K*- pi+ pi-
	  findDaughters(bCand,PiChp,PdtLund::pi_plus,
			PiChm,PdtLund::pi_minus,KstarCh,
			PdtLund::K_star_minus,false);
	}
	if((KstarCh!=0)&&(PiChp!=0)&&(PiChm!=0)) {
	  HepAList<BtaCandidate> KstarDescList;
	  findStableDescendents(KstarCh, KstarDescList);
	  if(KstarDescList.length()==2) {
	    BtaCandidate * KstarD1 = KstarDescList[1];
	    int KstarD1Id = abs(KstarD1->pdtEntry()->lundId());
	    if(KstarD1Id==PdtLund::K_plus||
	       KstarD1Id==PdtLund::pi0)  bBaseType = 201;
	  }
	}
      }

      //B+->K*0(K+ pi-) pi+ pi0
      if( bBaseType == -1) {
	BtaCandidate * Kstar0 = 0;
	BtaCandidate * PiChp = 0;
	BtaCandidate * Pi0 = 0;
	findDaughters(bCand, Kstar0, PdtLund::K_star0, 
		      PiChp, PdtLund::pi_plus, Pi0, PdtLund::pi0,true);
	if((Kstar0!=0)&&(PiChp!=0)&&(Pi0!=0)) {
	  HepAList<BtaCandidate> KstarDescList;
	  findStableDescendents(Kstar0, KstarDescList);
	  if(KstarDescList.length()==2) {
	    BtaCandidate * KstarD1 = KstarDescList[0];
	    int KstarD1Id = abs(KstarD1->pdtEntry()->lundId());
	    if(KstarD1Id==PdtLund::K_plus||
	       KstarD1Id==PdtLund::pi_plus)  bBaseType = 202;
	  }	
	}
      }

      //B+ -> K*+(K+ pi0) rho0(pi+ pi-) 
      if( bBaseType == -1 ) {
	BtaCandidate * KstarCh = 0;
	BtaCandidate * Rho0 = 0;
	findDaughters(bCand, KstarCh, PdtLund::K_star_plus, 
		      Rho0, PdtLund::rho0,true);
	if((KstarCh!=0)&&(Rho0!=0)) {
	  HepAList<BtaCandidate> KstarDescList;
	  HepAList<BtaCandidate> RhoDescList;
	  findStableDescendents(KstarCh, KstarDescList);
	  findStableDescendents(Rho0, RhoDescList);
	  if(RhoDescList.length()==2&&KstarDescList.length()==2) {
	    BtaCandidate * KstarD1 = KstarDescList[0];
	    BtaCandidate * RhoD1 = RhoDescList[0];
	    int KstarD1Id = abs(KstarD1->pdtEntry()->lundId());
	    int RhoD1Id = abs(RhoD1->pdtEntry()->lundId());
	    if((KstarD1Id==PdtLund::K_plus||
		KstarD1Id==PdtLund::pi0)&&(RhoD1Id==PdtLund::pi_plus)){
	      bBaseType = 203;
	    }
	  }
	}
      }

      //B+ -> K*0(K+ pi-) rho+(pi+ pi0)
      if( bBaseType == -1 ) {
	BtaCandidate * Kstar0 = 0;
	BtaCandidate * RhoC = 0;
	findDaughters(bCand, Kstar0, PdtLund::K_star0, 
		      RhoC, PdtLund::rho_plus,true);
	if((Kstar0!=0)&&(RhoC!=0)) {
	  HepAList<BtaCandidate> KstarDescList;
	  HepAList<BtaCandidate> RhoDescList;   
	  findStableDescendents(Kstar0, KstarDescList);
	  findStableDescendents(RhoC, RhoDescList); 
	  if(RhoDescList.length()==2&&KstarDescList.length()==2) {
	    BtaCandidate * KstarD1 = KstarDescList[0];
	    BtaCandidate * RhoD1 = RhoDescList[0];
	    int KstarD1Id = abs(KstarD1->pdtEntry()->lundId());
	    int RhoD1Id = abs(RhoD1->pdtEntry()->lundId());  
	    if((KstarD1Id==PdtLund::K_plus||
		KstarD1Id==PdtLund::pi_plus)&&
	       (RhoD1Id==PdtLund::pi_plus||RhoD1Id==PdtLund::pi0)){
	      bBaseType = 204;
	    }
	  }
	}
      }

      //B+ -> w(pi+ pi-) K*+(pi0 K+)
      if( bBaseType == -1 ) {
	BtaCandidate * w0 = 0;
	BtaCandidate * KstarCh = 0;
	findDaughters(bCand, KstarCh, PdtLund::K_star_plus, 
		      w0, PdtLund::omega,true);
	if((KstarCh!=0)&&(w0!=0)) {
	  HepAList<BtaCandidate> OmegaDescList;
	  HepAList<BtaCandidate> KstarDescList;
	  findStableDescendents(w0, OmegaDescList);	
	  findStableDescendents(KstarCh, KstarDescList);
	  BtaCandidate * OmgaD1 = OmegaDescList[0];
	  BtaCandidate * KstarD1 = KstarDescList[0];
	  int KstarD1Id = abs(KstarD1->pdtEntry()->lundId());
	  int OmgaD1Id = abs(OmgaD1->pdtEntry()->lundId());
	  if(KstarD1Id==PdtLund::pi_plus&&
	     (OmgaD1Id==PdtLund::K_plus||OmgaD1Id==PdtLund::pi0)) {
	    bBaseType = 205;
	  }     		 				 
	}
      }
	  
      //B+ -> pi- K+ rho+(pi0 pi+)
      if(bBaseType == -1) {
	BtaCandidate * PiCh=0;
	BtaCandidate * KCh = 0;
	BtaCandidate * RhoCh = 0;
	findDaughters(bCand,RhoCh,PdtLund::rho_plus, KCh, 
		      PdtLund::K_plus, PiCh, PdtLund::pi_plus,true);
	if(RhoCh!=0&&KCh!=0&&PiCh!=0) {
	  HepAList<BtaCandidate> RhoDescList;
	  findStableDescendents(RhoCh, RhoDescList);
	  if(RhoDescList.length()==2) {
	    BtaCandidate * RhoD1 = RhoDescList[0];
	    int RhoD1Id = abs(RhoD1->pdtEntry()->lundId());
	    if(RhoD1Id==PdtLund::pi_plus||
	       RhoD1Id==PdtLund::pi0) bBaseType = 206;
	  }		 
	}
      }
      
      //B+ -> pi0 K+ rho0(pi- pi+)
      if(bBaseType == -1) {
	BtaCandidate * Pi0=0;
	BtaCandidate * Kch = 0;
	BtaCandidate * Rho0 = 0;
	findDaughters(bCand,Rho0,PdtLund::rho0, Kch, 
		      PdtLund::K_plus, Pi0, PdtLund::pi0, true);
	if(Rho0!=0&&Kch!=0&&Pi0!=0) {
	  HepAList<BtaCandidate> RhoDescList;
	  findStableDescendents(Rho0, RhoDescList);
	  if(RhoDescList.length()==2) {
	    BtaCandidate * RhoD1 = RhoDescList[0];
	    int RhoD1Id = abs(RhoD1->pdtEntry()->lundId());
	    if(RhoD1Id==PdtLund::pi_plus) bBaseType = 207;
	  }	
	}
      }

      //B+ -> rho- K+ pi+
      if(bBaseType == -1) {
	BtaCandidate * Pich=0;
	BtaCandidate * Kch = 0;
	BtaCandidate * RhoC = 0;
	findDaughters(bCand,RhoC,PdtLund::rho_plus, 
		      Kch, PdtLund::K_plus, Pich, PdtLund::pi_plus, true);
	if(RhoC!=0&&Kch!=0&&Pich!=0) {
	  HepAList<BtaCandidate> RhoDescList;
	  findStableDescendents(RhoC, RhoDescList);
	  if(RhoDescList.length()==2) {
	    BtaCandidate * RhoD1 = RhoDescList[0];
	    int RhoD1Id = abs(RhoD1->pdtEntry()->lundId());
	    if(RhoD1Id==PdtLund::pi_plus) bBaseType = 208;
	  }
	}
      }

      //====for B+ -> D0(K+ K- pi0) K+ mode ========
      // B+ -> K*+(K+ pi0/K0 pi+) K+ K-
      if( bBaseType == -1 ) {
	BtaCandidate * KstarCh = 0;
	BtaCandidate * Kchp = 0;
	BtaCandidate * Kchm = 0;
	findDaughters(bCand, KstarCh, PdtLund::K_star_plus, 
		      Kchp, PdtLund::K_plus, Kchm, PdtLund::K_minus,false);
	if((KstarCh==0)||(Kchp==0)||(Kchm==0)) {
	  findDaughters(bCand, KstarCh, PdtLund::K_star_minus, Kchp, 
			PdtLund::K_plus, Kchm, PdtLund::K_minus,false);
	}	
	if((KstarCh!=0)&&(Kchp!=0)&&(Kchm!=0)) {
	  BtaCandidate * Kch = 0;
	  BtaCandidate * Pi0 = 0;
	  BtaCandidate * K0 = 0;
	  BtaCandidate * PiCh = 0;            
	  findDaughters(KstarCh,Kch, PdtLund::K_plus, Pi0, 
			PdtLund::pi0, true);
	  if((Kch!=0)&&(Pi0!=0)) bBaseType = 301;//B+ -> (K+ K- pi0) K+
	  
	  findDaughters(KstarCh,K0,PdtLund::K_S0, PiCh, 
			PdtLund::pi_plus, true);
	  if((K0!=0)&&(PiCh!=0)) bBaseType = 401;//B+->(K0 pi+ K-) K+
	}
      }

      // B+ -> K*+(K+ pi0/ K0 pi+) phi(K+ K-)
      if( bBaseType == -1 ) {
	BtaCandidate * KstarCh = 0;
	BtaCandidate * phi = 0;
	findDaughters(bCand, KstarCh, PdtLund::K_star_plus, 
		      phi, PdtLund::phi, true);
	if((KstarCh!=0)&&(phi!=0)) {
	  BtaCandidate * K1 = 0;
	  BtaCandidate * K2 = 0;
	  BtaCandidate * Kt = 0;
	  BtaCandidate * Pit = 0;
	  findDaughters(phi, K1, PdtLund::K_plus, K2, 
			PdtLund::K_minus, false);          
	  findDaughters(KstarCh, Kt, PdtLund::K_plus, 
			Pit, PdtLund::pi0);
	  if((K1!=0)&&(K2!=0)&&(Kt!=0)&&(Pit!=0)) 
	    bBaseType = 302;//B+ -> (K+ K- pi0) K+
	  if( bBaseType ==-1) {
	    findDaughters(KstarCh, Kt, PdtLund::K_S0, 
			  Pit, PdtLund::pi_plus, true);
	    if((K1!=0)&&(K2!=0)&&(Kt!=0)&&(Pit!=0)) 
	      bBaseType = 402;//B+ -> (K0 pi+ K-) K+
	  }
	}
      }    

      // B+ -> K*+1430(K+ pi0/ K0 pi+) phi(K+ K-)
      if( bBaseType == -1 ) {
	BtaCandidate * K2starCh = 0;
	BtaCandidate * phi = 0;
	findDaughters(bCand, K2starCh, PdtLund::K_2_star_plus, 
		      phi, PdtLund::phi, true);
	if((K2starCh!=0)&&(phi!=0)) {
	  BtaCandidate * K1 = 0;
	  BtaCandidate * K2 = 0;
	  BtaCandidate * Kt = 0;
	  BtaCandidate * Pit = 0;
	  findDaughters(phi, K1, PdtLund::K_plus, K2, 
			PdtLund::K_minus, false);
	  findDaughters(K2starCh, Kt, PdtLund::K_plus, 
			Pit, PdtLund::pi0);
	  if((K1!=0)&&(K2!=0)&&(Kt!=0)&&(Pit!=0)) 
	    bBaseType = 303;//B+ -> (K+ K- pi0) K+
	  if( bBaseType ==-1) {
	    findDaughters(K2starCh, Kt, PdtLund::K_S0, Pit, 
			  PdtLund::pi_plus, true);
	    if((K1!=0)&&(K2!=0)&&(Kt!=0)&&(Pit!=0)) 
	      bBaseType = 403;//B+ -> (K0 pi+ K-) K+
	  }
	}
      }

      //B+ -> phi K+ pi0
      if( bBaseType == -1 ) {
	BtaCandidate * Kch = 0;
	BtaCandidate * phi = 0;
	BtaCandidate * pi0 = 0;
	findDaughters(bCand, Kch, PdtLund::K_plus, phi, 
		      PdtLund::phi, pi0, PdtLund::pi0,true);
	if((Kch!=0)&&(phi!=0)&&(pi0!=0)) {
	  BtaCandidate * K1 = 0;
	  BtaCandidate * K2 = 0;
	  findDaughters(phi, K1, PdtLund::K_plus, K2, 
			PdtLund::K_minus, false);
	  if((K1!=0)&&(K2!=0)) bBaseType = 304;//B+ -> (K+ K- pi0) K+
	}
      }
		
      //B+ -> f_1(1285) K+
      if( bBaseType == -1) {
	BtaCandidate * f1neu = 0;
	BtaCandidate * Kch = 0;
	findDaughters(bCand, f1neu, PdtLund::f_1, Kch, 
		      PdtLund::K_plus, true);
	if(f1neu != 0 && Kch !=0 ) {
	  BtaCandidate * Pch1 = 0;
	  BtaCandidate * Pch2 = 0;
	  BtaCandidate * Pneu = 0; 
	  findDaughters(f1neu,Pch1, PdtLund::K_plus, Pch2, 
			PdtLund::K_minus, Pneu, PdtLund::pi0, false);
	  if(Pch1 != 0 && Pch2 != 0 && Pneu != 0 ) 
	    bBaseType = 305; //B+ -> (K+ K- pi0) K+
	  Pch1=0; Pch2=0; Pneu=0;
	  findDaughters(f1neu,Pch1, PdtLund::K_plus, Pch2, 
			PdtLund::pi_plus, Pneu, PdtLund::K_S0, true); 
	  if(Pch1 != 0 && Pch2 != 0 && Pneu != 0 ) 
	    bBaseType = 501; //B+ -> (K0 pi- K+) K+
	}
      }  

      //B+ -> Ds+(K+ Ks0) K*0(K+ pi-/ K- pi+)   
      if( bBaseType == -1 ) {
	BtaCandidate * DsubS = 0;
	BtaCandidate * Kstar0 = 0;
	findDaughters(bCand, DsubS, PdtLund::D_s_plus, Kstar0, 
		      PdtLund::K_star0,true);
	if((DsubS!=0)&&(Kstar0!=0)) {
	  //B+ -> (K+ pi- Ks) K+
	  BtaCandidate * Kchs1 = 0;
	  BtaCandidate * Kchs2 = 0;
	  findDaughters(DsubS, Kchs1, PdtLund::K_plus, Kchs2, 
			PdtLund::K_S0, false);
	  BtaCandidate * kCand = 0;
	  BtaCandidate * piCand = 0;
	  findDaughters(Kstar0, kCand, PdtLund::K_plus, piCand, 
			PdtLund::pi_minus, false);
	  if((Kchs1!=0)&&(Kchs2!=0)&&(kCand!=0)&&(piCand!=0)) 
	    bBaseType=502;
	  if( bBaseType== -1) {
	    //B+ -> (K- pi+ Ks) K+    
	    findDaughters(Kstar0, kCand, PdtLund::K_minus, piCand, 
			  PdtLund::pi_plus, false); 	        	
	    if((Kchs1!=0)&&(Kchs2!=0)&&(kCand!=0)&&(piCand!=0)) 
	      bBaseType=404;
	  }
	}
      }

      //    if(bBaseType>0)  cout<<"B decay Base Type "<<bBaseType<<endl;
      if(bBaseType < 0) continue;	
      
      BtaCandidate * mcPip = 0;
      BtaCandidate * mcPim = 0;
      BtaCandidate * mcPi0 = 0;
      BtaCandidate * mcKp = 0;
      BtaCandidate * mcKm = 0;
      BtaCandidate * mcKs = 0;
      
      HepLorentzVector KplusMCP4;
      HepLorentzVector KminusMCP4;
      HepLorentzVector KsMCP4;
      HepLorentzVector PiplusMCP4;
      HepLorentzVector PiminusMCP4;
      HepLorentzVector pi0MCP4;
      HepLorentzVector D0MCP4;

      d0XType = 0;
      
      if(bBaseType >=0 && dCand != 0  ) { // find the decay type of D0
	//try to find D0 -> K+ K- pi0
	if (0 == d0XType ) { 	        
	  findD0Descendents(dCand,mcKp,PdtLund::K_plus, 
			    mcKm,PdtLund::K_minus,
			    mcPi0,PdtLund::pi0);
	  if((mcKp != 0)&&(mcKm != 0)&&(mcPi0 != 0))  {
	    d0XType = 1;
	    
	    // compond mass
	    KplusMCP4 = mcKp->p4();
	    KminusMCP4 = mcKm->p4();
	    pi0MCP4 = mcPi0->p4();
	    D0MCP4 = dCand->p4();
	    
	    // K+ K-
	    HepLorentzVector KKcompMCP4;
	    KKcompMCP4 = KplusMCP4 + KminusMCP4;
	    Kch2MCMass = KKcompMCP4.mag();
	    Kch2MCHelic = calcHelicity(KplusMCP4, KKcompMCP4, D0MCP4);
	    //K+ Pi0
	    HepLorentzVector KpPi0compMCP4;
	    KpPi0compMCP4 = KplusMCP4 + pi0MCP4;
	    KpPi0MCMass = KpPi0compMCP4.mag();
	    KpPi0MCHelic = calcHelicity(KplusMCP4, KpPi0compMCP4, D0MCP4);
	    //K- pi0
	    HepLorentzVector KmPi0compMCP4;
	    KmPi0compMCP4 = KminusMCP4 + pi0MCP4;
	    KmPi0MCMass = KmPi0compMCP4.mag();
	    KmPi0MCHelic = calcHelicity(KminusMCP4, KmPi0compMCP4,D0MCP4);
	  }
	}

	//try to find D0 -> K+ pi- K_S0
	if (0 == d0XType ) { 
	  int D0mcId = 0; 
	  findD0Descendents(dCand,mcKp,PdtLund::K_plus, 
			    mcPim, PdtLund::pi_minus,
			    mcKs, PdtLund::K_S0);
	  if((mcKp != 0)&&(mcPim != 0)&&(mcKs != 0)) {
	    D0mcId = dCand->pdtEntry()->lundId();
	    if(D0mcId == PdtLund::D0) d0XType = 2;
	    if(D0mcId == PdtLund::anti_D0) d0XType = 3;
	    
	    // compond mass
	    KplusMCP4 = mcKp->p4();
	    PiminusMCP4 = mcPim->p4();
	    KsMCP4 = mcKs->p4();
	    D0MCP4 = dCand->p4();
	    
	    //K+ Pi- 
	    HepLorentzVector KpPicompMCP4;
	    KpPicompMCP4 = KplusMCP4 + PiminusMCP4;
	    KpPi2MCMass = KpPicompMCP4.mag();
	    KpPi2MCHelic = calcHelicity(KplusMCP4,KpPicompMCP4, D0MCP4);
	    
	    //K+ K_S0
	    HepLorentzVector KpKscompMCP4;
	    KpKscompMCP4 = KplusMCP4 + KsMCP4;
	    KpKsMCMass = KpKscompMCP4.mag();
	    KpKsMCHelic = calcHelicity(KplusMCP4,KpKscompMCP4,D0MCP4);
	    //pi- K_S0
	    HepLorentzVector PimKscompMCP4;
	    PimKscompMCP4 = PiminusMCP4 + KsMCP4;
	    PimKsMCMass = PimKscompMCP4.mag();
	    PimKsMCHelic = calcHelicity(PiminusMCP4,PimKscompMCP4,D0MCP4);
	  }
	}

	//try to find D0 -> K- pi+ K_S0
	if (0 == d0XType ) {
	  int D0mcId = 0;
	  findD0Descendents(dCand,mcKm,PdtLund::K_minus, 
			    mcPip, PdtLund::pi_plus,
			    mcKs, PdtLund::K_S0);
	  if((mcKm != 0)&&(mcPip != 0)&&(mcKs != 0)) {
	    D0mcId = dCand->pdtEntry()->lundId();				
	    if(D0mcId == PdtLund::D0) d0XType = 4;
	    if(D0mcId == PdtLund::anti_D0) d0XType = 5;
	    
	    // compond mass
	    KminusMCP4 = mcKm->p4();
	    PiplusMCP4 = mcPip->p4();
	    KsMCP4 = mcKs->p4();
	    D0MCP4 = dCand->p4();
	    
	    //K- Pi+ 
	    HepLorentzVector KmPicompMCP4;
	    KmPicompMCP4 = KminusMCP4 + PiplusMCP4;
	    KmPi2MCMass = KmPicompMCP4.mag();
	    KmPi2MCHelic = calcHelicity(KminusMCP4,KmPicompMCP4,D0MCP4);
	    //K- K_S0
	    HepLorentzVector KmKscompMCP4;
	    KmKscompMCP4 = KminusMCP4 + KsMCP4;
	    KmKsMCMass = KmKscompMCP4.mag();
	    KmKsMCHelic = calcHelicity(KminusMCP4,KmKscompMCP4,D0MCP4);
	    //pi+ K_S0
	    HepLorentzVector PipKscompMCP4;
	    PipKscompMCP4 = PiplusMCP4 + KsMCP4;
	    PipKsMCMass = PipKscompMCP4.mag();
	    PipKsMCHelic = calcHelicity(PiplusMCP4,PipKscompMCP4,D0MCP4);
	  }
	}

	//try to find D0 -> pi+ pi- pi0
	if (0 == d0XType ) {
	  findD0Descendents(dCand,mcPim,PdtLund::pi_minus, 
			    mcPip, PdtLund::pi_plus,
			    mcPi0, PdtLund::pi0);
	  if((mcPim != 0)&&(mcPip != 0)&&(mcPi0 != 0)) {
	    d0XType = 6;
	         
	    // compond mass
	    PiminusMCP4 = mcPim->p4();
	    PiplusMCP4 = mcPip->p4();
	    pi0MCP4 = mcPi0->p4();
	    D0MCP4 = dCand->p4();
	    
	    //pi+ pi-
	    HepLorentzVector PiPiMCP4;
	    PiPiMCP4 = PiminusMCP4 + PiplusMCP4;
	    PiPiMCMass = PiPiMCP4.mag();
	    PiPiMCHelic = calcHelicity(PiplusMCP4,PiPiMCP4,D0MCP4);
	    //pi+ pi0
	    HepLorentzVector PipPi0MCP4;
	    PipPi0MCP4 = PiplusMCP4 + pi0MCP4;
	    PipPi0MCMass = PipPi0MCP4.mag();
	    PipPi0MCHelic = calcHelicity(PiplusMCP4, PipPi0MCP4, D0MCP4);
	    //pi- pi0
	    HepLorentzVector PimPi0MCP4;
	    PimPi0MCP4 = PiminusMCP4 + pi0MCP4;
	    PimPi0MCMass = PimPi0MCP4.mag();
	    PimPi0MCHelic = calcHelicity(PiminusMCP4,PimPi0MCP4, D0MCP4);
	  }
	}

        //try to find D0 -> Ks(pi+ pi-) pi0
        if (0 == d0XType ) {
          findD0Descendents(dCand, mcKs, PdtLund::K_S0,
                            mcPi0, PdtLund::pi0);
          if((mcKs != 0)&&(mcPi0 != 0)) {
            d0XType = 7;
          }
        }
		 
	// to see D0 decay to cabibbo allowed mode
	// 1. D0 -> Ks K+ K-, 2. D0 -> Ks pi+ pi-, 4. D0 -> pi0 K- pi+, 
	// 3. cabibbo supressed D0-> pi0 K+ pi- 
	if(0 == d0XType ) {
	  BtaCandidate* Ddau1(0);
	  BtaCandidate* Ddau2(0);
	  BtaCandidate* Ddau3(0);
	  
	  findD0Descendents(dCand,Ddau1,PdtLund::K_plus,Ddau2, 
			    PdtLund::K_minus,Ddau3,PdtLund::K_S0);
	  if(0 != Ddau1 && 0 != Ddau2 && 0 != Ddau3 ) D0CADecay = 1;
	  findD0Descendents(dCand,Ddau1,PdtLund::pi_plus,Ddau2, 
			    PdtLund::pi_minus,Ddau3,PdtLund::K_S0);
	  if(0 != Ddau1 && 0 != Ddau2 && 0 != Ddau3 ) D0CADecay = 2;
	  findD0Descendents(dCand,Ddau1,PdtLund::K_plus, Ddau2, 
			    PdtLund::pi_minus,Ddau3,PdtLund::pi0);
	  if(0 != Ddau1 && 0 != Ddau2 && 0 != Ddau3 ) D0CADecay = 3;
	  findD0Descendents(dCand,Ddau1,PdtLund::pi_plus,Ddau2, 
			    PdtLund::K_minus,Ddau3,PdtLund::pi0);
	  if(0 != Ddau1 && 0 != Ddau2 && 0 != Ddau3 ) D0CADecay = 4;	
	}	
      } // end of bBaseType >= 0
      
      // here save these above mc decay types in histogram monte carlo
      if( bBaseType >= 0) BdecModMC[numBFound-1] = d0XType+bBaseType;
      if(D0CADecay > 0 ) D0CADecMode[numBFound-1] = D0CADecay;
      
      bCand = 0;
      bBaseType=-1;
      
      if( 2 == numBFound) {
	break; //Found the 2 B meson. Do not check further
      }    
    } // end of ==BID
    
  } // end of while
	  
  double mcbDistZ(-100.0);	  
  if(bz[0] != 0. || bz[1] != 0. ) {
    mcbDistZ    = bz[0] - bz[1];
  }

  //
  // B+ -> D0 K+  exclusive analysis:
  //
  if( BchD0KchList->length() == 0 ){
    if(_verbose.value()){
      ErrMsg(routine) << name()
		      << " reconstructed BchD0KchList is empty."
		      << endmsg;
    }
    return AppResult::OK ; //quit the program
  }

  // construct another list,be used for putting more than one singal list 
  // A list pointing to all the exclusive candidates:
  HepAList<BtaCandidate> exclList;
  for (int i = 0; i < BchD0KchList->length(); i++){
    exclList.append((*BchD0KchList)[i]);
  }
  
  // The vector of candidates for this event:
  vector<D0KchCand> d0KchCandidates;
  
  const BtaCandidate * exclBch = 0;
  static const HepLorentzVector labP4(0,0,0,_massUpsilon);
  
  for (int i = 0; i < exclList.length(); i++){
    exclBch = exclList[i];
    if (0 == exclBch){
      continue;
    }

    // Initialize the B+ -> D0 K+ candidate:
    _d0KchCand = D0KchCand();

    // add MC dalitz 
    _d0KchCand._d0kkMCMass  = Kch2MCMass;
    _d0KchCand._d0pKpMCMass = KpPi0MCMass;
    _d0KchCand._d0pKmMCMass = KmPi0MCMass;
    _d0KchCand._d0kkMCHelic  = Kch2MCHelic;
    _d0KchCand._d0pKpMCHelic = KpPi0MCHelic;
    _d0KchCand._d0pKmMCHelic = KmPi0MCHelic;  
    
    _d0KchCand._d0pKpMCMass = KpPi2MCMass;
    _d0KchCand._d0kkMCMass  = KpKsMCMass; 
    _d0KchCand._d0pKmMCMass = PimKsMCMass;
    _d0KchCand._d0pKpMCHelic = KpPi2MCHelic;
    _d0KchCand._d0kkMCHelic  = KpKsMCHelic;
    _d0KchCand._d0pKmMCHelic = PimKsMCHelic;

    _d0KchCand._d0pKmMCMass = KmPi2MCMass;
    _d0KchCand._d0kkMCMass  = KmKsMCMass;
    _d0KchCand._d0pKpMCMass = PipKsMCMass;
    _d0KchCand._d0pKmMCHelic = KmPi2MCHelic;
    _d0KchCand._d0kkMCHelic  = KmKsMCHelic;
    _d0KchCand._d0pKpMCHelic  = PipKsMCHelic;
    
    _d0KchCand._d0ppMCMass  = PiPiMCMass;
    _d0KchCand._d0pPpMCMass = PipPi0MCMass;
    _d0KchCand._d0pPmMCMass = PimPi0MCMass;
    _d0KchCand._d0ppMCHelic  = PiPiMCHelic;
    _d0KchCand._d0pPpMCHelic = PipPi0MCHelic;
    _d0KchCand._d0pPmMCHelic = PimPi0MCHelic;
		       
    //find B+ ->D0 K+ and get decay mode of D0 
    BtaCandidate * exclD0 = 0;
    BtaCandidate * exclKch = 0;
    
    findDaughters( exclBch, exclD0,PdtLund::D0,
		   exclKch,PdtLund::K_plus,true);
    if( (0 == exclKch) || (0 == exclD0) ) continue;

    //kaon pid declaration here:
    static PidKaonSMSSelector KaonSMSVt = PidKaonSMSSelector();
    static PidKaonSMSSelector KaonSMSTi = PidKaonSMSSelector();
    static PidKaonSMSSelector KaonSMSLs = PidKaonSMSSelector();
    static PidKaonSMSSelector KaonSMSVl = PidKaonSMSSelector();
    static PidKaonSMSSelector KaonSMSNp = PidKaonSMSSelector();
    static bool KaonIntializer = false;
    if(false == KaonIntializer) {
      KaonSMSVt.setParmValue("criteria","veryTight");
      KaonSMSTi.setParmValue("criteria","tight");
      KaonSMSLs.setParmValue("criteria","loose");
      KaonSMSVl.setParmValue("criteria","veryLoose");
      KaonSMSNp.setParmValue("criteria","notApion");
      KaonIntializer = true;
    }
	
    //set up the pion LH for all following pions
    static PidPionLHSelector pionLHSVl = PidPionLHSelector();
    static PidPionLHSelector pionLHSLs = PidPionLHSelector();
    static PidPionLHSelector pionLHSTi = PidPionLHSelector();
    static PidPionLHSelector pionLHSVt = PidPionLHSelector();
    static bool pionIntializer = false;
    if(false == pionIntializer) {
      pionLHSVl.setParmValue("criteria","veryLoose");
      pionLHSLs.setParmValue("criteria","loose");
      pionLHSTi.setParmValue("criteria","tight");
      pionLHSVt.setParmValue("criteria","veryTight");
      pionIntializer = true;
    }
       
    //get D0 daughter list
    HepAList<BtaCandidate> d0DauList;
    HepAListIterator<BtaCandidate> d0DauIter = exclD0->daughterIterator();
    d0DauIter.rewind();
    BtaCandidate * d0Dau(0);
    while(d0Dau = d0DauIter() ) {
      d0DauList += d0Dau;
    }

    //vertexing D0 candidate       
    BtaCandidate fittedD0(0);
    CandVtxFit(d0DauList, false, false, 0.0);
    double d0vtFitChisq(-100.0);
    if(_fittedCand.decayVtx()->status() == BtaAbsVertex::Success) { 
      fittedD0 = _fittedCand;
      d0vtFitChisq =  fittedD0.decayVtx()->chiSquared();
    }
    _d0KchCand._d0VtxFitChisq = d0vtFitChisq;

    // D0 cadidate with mass constraint --JZ
    BtaCandidate fittedMSD0(0);
    CandVtxFit(d0DauList,false, true,Pdt::mass("D0"));
    double d0vtFitMSChisq(-100.0);
    if(_fittedCand.decayVtx()->status() == BtaAbsVertex::Success) { 
      fittedMSD0 = _fittedCand;
      d0vtFitMSChisq =  fittedMSD0.decayVtx()->chiSquared();
    }
    
    _d0KchCand._d0VtxFitMSChisq = d0vtFitMSChisq;
    if( (fittedD0==0) || (fittedMSD0==0) ) continue;
    
    // get D0 mass
    HepLorentzVector D0RawP4 = fittedD0.p4();
    double exclD0Mass = 0.0;
    exclD0Mass = D0RawP4.mag();
    //D0 mass cut here
    if( false == inRange((exclD0Mass-_massD0),_exclD0MassMin.value(),
			 _exclD0MassMax.value())){
      continue;
    }	   
	
    //D0 CM momentum  
    HepLorentzVector D0CMP4 = D0RawP4;
    D0CMP4.boost(_boostVector);
    const double D0CMMomtum = D0CMP4.vect().mag();
    
    _d0KchCand._d0Mass = exclD0Mass;
    _d0KchCand._d0CMMomtum = D0CMMomtum;
    
    // Get this event's deltaE & mES:
    double exclMES = 0.0;
    double exclDeltaE = -1.0;
    
    HepLorentzVector D0FitP4 = fittedMSD0.p4(); 
    HepLorentzVector exclBRawP4 = D0FitP4+exclKch->p4();
    //Boost the B candidate to the CM frame:
    HepLorentzVector p4CM = exclBRawP4;
    p4CM.boost(_boostVector);
    
    // Get its deltaE & mES:
    exclDeltaE = p4CM.t() - _eBeam;
    exclMES = sqrt(_eBeam * _eBeam - p4CM.vect().mag2());
    //Delta E cut here		
    if( false ==inRange(exclDeltaE,_exclDeltaEMin.value(),
			_exclDeltaEMax.value())) {
      continue;
    }
    //Mes cut here
    if(false == inRange(exclMES, _exclMESMin.value(),
			_exclMESMax.value()) ) {
      continue;
    } 
	    
    // Momentum angle in CM
    double BCMcosTheta = p4CM.cosTheta();

    _d0KchCand._exclDE = exclDeltaE;
    _d0KchCand._Mes = exclMES;
    _d0KchCand._BCMcosTheta = BCMcosTheta;

    //fisher discriminant
    static Q2BCHBEventShape D0KchFisher = Q2BCHBEventShape();
    static bool setList = false;
    
    if(setList == false) {     
      D0KchFisher.setChargedList("GoodTracksVeryLoose");
      D0KchFisher.setNeutralsList("GoodPhotonLoose");
      setList = true;
    }

    D0KchFisher.setLegendreOn(9);
    D0KchFisher.compute(*exclBch, anEvent, true);
    int FisherSt = D0KchFisher.getStatus();
    // 0 success, 1 no B candidate, 2 only B candidate
    if( FisherSt == 0 ) {
      double Bfisher = D0KchFisher.getFisher();
      //cleo fisher
      double BfCleo = D0KchFisher.getFisherCLEO();
      // the thrust for B candidate only       
      double BcanThr = D0KchFisher.getBthrust();
      /* cosine of the angle between the thrust axis of the B 
	 candidate and the thrust axis of the rest of the event */
      double BcanCosThetaT = D0KchFisher.getCosThetaT();
      // cosine of the angle between the B direction and the beam axis
      double CosBmomAx = D0KchFisher.getCosBmom();
      /* cosine of the angle between the thrust axis of the B candidate
	 and the beam axis */
      double CosBThxAx = D0KchFisher.getCosBthr();
      // all above get member function is calculate from upslion(4s) COM
	       
      _d0KchCand._BFisher = Bfisher;
      _d0KchCand._BFCLeo  = BfCleo;
      _d0KchCand._BCanThrust = BcanThr;
      _d0KchCand._BCanCosThR = BcanCosThetaT;
      _d0KchCand._cosBmomAx = CosBmomAx;
      _d0KchCand._cosBThRAx = CosBThxAx;		 		  
      _d0KchCand._BFLegendre = D0KchFisher.getLegendreFisher();
    }

    //Hemisphere charge : use D0 Momentum vector as a reference
    double QDiff = -100.0;

    double Bchrge = exclBch->charge();
    HepLorentzVector MomRefFame = D0FitP4;
    MomRefFame.boost(_boostVector);//cm frame
    Hep3Vector D03Vect = MomRefFame.vect();
    
    //go through all the trks in the trklist
    HepAListIterator<BtaCandidate> IterTrkList(*btaTrkList);
    BtaCandidate * TrkCand(0);
    IterTrkList.rewind();
    double Hem1Chrge = 0.0;
    double Hem2Chrge = 0.0;
    int lpflg = 0;		
    while( 0 != (TrkCand = IterTrkList()) ) {
      if(TrkCand->overlaps(*exclKch)) continue;
      HepLorentzVector TrkP4cm = TrkCand->p4();
      double TrkChrge = TrkCand->charge();
      TrkP4cm.boost(_boostVector);
      Hep3Vector Trk3Vect = TrkP4cm.vect();
      double cosD0Trkcm = cos(D03Vect.angle(Trk3Vect));
      if(cosD0Trkcm>0) Hem1Chrge = Hem1Chrge + TrkChrge;
      if(cosD0Trkcm<0) Hem2Chrge = Hem2Chrge + TrkChrge;
      if(lpflg == 0) lpflg = 1;
    }
    if( lpflg == 1) QDiff = Bchrge*(Hem1Chrge - Hem2Chrge);
    _d0KchCand._QDiff = QDiff;
		
    //get Bch daughter list
    HepAList<BtaCandidate> bDauList;
    HepAListIterator<BtaCandidate> BDauIter = exclBch->daughterIterator();
    BDauIter.rewind();
    BtaCandidate * bDau(0);
    while(bDau = BDauIter() ) {
      bDauList += bDau;
    }

    //vertexing B candidate
    CandVtxFit(bDauList,true, false, 0.0);
    BtaCandidate fittedB(0);
    double BvtFitChisq(-100.0);
    if(_fittedCand.decayVtx()->status() == BtaAbsVertex::Success) {
      fittedB  = _fittedCand;
      BvtFitChisq =  fittedB.decayVtx()->chiSquared();
    }
    _d0KchCand._bVtxFitChisq = BvtFitChisq;
    //cout<<" B Vtx Fit freedom "<<fittedB.decayVtx()->nDof()<<endl;
    
    
    /*  get the angle between the thrust axis of D daughters and 
	the thrust axis/Mom of B in D cm */
    EventInfo D0eventInfo = EventInfo();
    D0eventInfo.setCmFrame(D0FitP4);
    BtaThrust D0Thr(d0DauList, D0eventInfo, BtaThrust::BTAllParticles);
    BtaThrust BThr(bDauList, D0eventInfo, BtaThrust::BTAllParticles);
    double cosBDThXDFm = -100.0;
    cosBDThXDFm = cos(D0Thr.thrust_axis().angle(BThr.thrust_axis()));
    _d0KchCand._cosBDThXDFm = cosBDThXDFm;
    
    HepLorentzVector BchD0Fm4Vector = exclBRawP4;
    BchD0Fm4Vector.boost(D0FitP4.boostVector());
    Hep3Vector BchD0Fm3Vect = BchD0Fm4Vector.vect();
    double cosBmomThrDFm = -100.0;
    cosBmomThrDFm = cos(D0Thr.thrust_axis().angle(BchD0Fm3Vect));
    _d0KchCand._cosBmomDThrDFm = cosBmomThrDFm;
			
    //distance and chisquare of D vertex from beam spot in y
    double d0VtxDY(-1000.0);
    double d0VtxDYChisq(-1000.0);
    double KbD0Doca(-100.0);
    double bdVtx3VectChisq(-100.0);
    double d0FlightDist(-100.0);
    double cosBDVtxDmom(-100.0);
    HepSymMatrix beamSpotMatxY = _beamSpot.covMatrix();
    double beamSpotCovY = beamSpotMatxY(2,2);
    if(fittedD0 != 0) {
      d0VtxDY = fittedD0.decayVtx()->point().y() - _beamSpot.y();
      //error matrix of D0 Vtx point
      HepSymMatrix dVtxCoV = fittedD0.decayVtx()->covariance();
      double d0VtxCovY = dVtxCoV(2,2);
      d0VtxDYChisq = (d0VtxDY*d0VtxDY)/(d0VtxCovY + beamSpotCovY);
      //D0 decay veterx point	  
      const HepPoint dVtx = fittedD0.decayVtx()->point(); 

      //Doca and chisquare for K and D "tracks"?????
      TrkLineTraj d0traj( dVtx, D0FitP4.vect(),-999., 999.);
      TrkPoca poca(*((Trajectory*) &d0traj), 0., 
		   *((Trajectory*)&(exclKch->trkAbsFit()->traj())), 0. );  
      KbD0Doca = poca.doca();
      
      if(fittedB != 0) {
	//B vertex piont
	const HepPoint bVtx = fittedB.decayVtx()->point();
	Hep3Vector bdVtx3Vect = bVtx - dVtx;
	d0FlightDist = bdVtx3Vect.mag();//d0 flight distance
	// to find the chisquar for the B, D0 decay vertex length
	//error matrix of B vtx point
	HepSymMatrix bVtxCoV = fittedB.decayVtx()->covariance();
	BbrVectorErr bdVtx3VectErr(bdVtx3Vect);
	bdVtx3VectErr.setCovMatrix(bVtxCoV + dVtxCoV);
	bdVtx3VectChisq = bdVtx3VectErr.determineChisq(Hep3Vector());
	cosBDVtxDmom = cos(bdVtx3Vect.angle(D0FitP4.vect()));
      }
    }
    _d0KchCand._d0VtxDY = d0VtxDY; 
    _d0KchCand._d0VtxDYChisq = d0VtxDYChisq;
    _d0KchCand._KbD0Doca = KbD0Doca;  
    _d0KchCand._cosBDVtxDmom = cosBDVtxDmom;
    _d0KchCand._bdVtxVctChisq = bdVtx3VectChisq;	
    _d0KchCand._d0FlightDist = d0FlightDist;

    // delta Z and its error for B vtx point and vtx point of 
    // rest of the event
    HepAList<BtaCandidate> restOfBEvent;
    HepAList<BtaCandidate> restOfBnoKEvent;
    HepAList<BtaCandidate> roeTrkList;
    HepAList<BtaCandidate> TgRoeNeutList;
    IterTrkList.rewind();//tracklist 
    BtaCandidate * RoeTrkCand(0);
    while(0 !=( RoeTrkCand = IterTrkList()) ) {
      if(!RoeTrkCand->overlaps(*exclBch)) {
	restOfBEvent.append(*RoeTrkCand);
	roeTrkList.append(*RoeTrkCand);
      }
    }

    HepAListIterator<BtaCandidate> neutralIter(*gammaList);
    BtaCandidate * RoeNeuCand(0);	
    neutralIter.rewind();
    while(0 != ( RoeNeuCand = neutralIter()) ) {
      if(!RoeNeuCand->overlaps(*exclBch)) 
	restOfBEvent.append(*RoeNeuCand);
      if(!RoeNeuCand->overlaps(*exclBch)) 
	restOfBnoKEvent.append(*RoeNeuCand);
      if(!RoeNeuCand->overlaps(*exclBch)) 
	TgRoeNeutList.append(*RoeNeuCand);
    }

    // Ks loose
    HepAList<BtaCandidate>* KsLsList;
    KsLsList = Ifd< HepAList<BtaCandidate> >::get(anEvent,
						  _btaKsLsList.value());
    HepAList<BtaCandidate>* KsList;
    KsList = Ifd< HepAList<BtaCandidate> >::get(anEvent,
						_btaKsList.value());
    HepAListIterator<BtaCandidate> KsIterator(*KsList);
    
    HepAList<BtaCandidate> KsRoEList;
    BtaCandidate * tmpCand(0);
    KsIterator.rewind();
    while(0 != (tmpCand = KsIterator()) ) {
      if( ! tmpCand->overlaps(*exclBch)) KsRoEList.append(*tmpCand);
    }
  
    //tag-use roe list    
    HepAList<BtaCandidate> TgRoeTrkList;
    HepAListIterator<BtaCandidate> TgRoEUseIter(*MicroTaggingList);
    BtaCandidate * TagRoeCand(0);
    TgRoEUseIter.rewind();
    while(0 != ( TagRoeCand = TgRoEUseIter() ) ) {
      if(!TagRoeCand->overlaps(*exclBch)) TgRoeTrkList.append(*TagRoeCand);
    }
    HepAList<BtaCandidate> TrkRoeNoKsList;
    IterTrkList.rewind();
    RoeTrkCand = 0;
    while( 0 != (RoeTrkCand = IterTrkList()) ) {
      BtaCandidate * tmpKsCand(0);
      KsIterator.rewind();
      while(0 != (tmpKsCand = KsIterator()) ) {
	if(! RoeTrkCand->overlaps(*tmpKsCand) )
	  TrkRoeNoKsList.append(*RoeTrkCand);	
      }
    } 
 
    //vertex rest of the Event
    double BandBRoeY(-1000.0);
    double BandBRoeZ(-1000.0);
    double BandBRoeYerr(-1000.0);
    double BandBRoeZerr(-1000.0);
    HepAListIterator<BtaCandidate> RoEListIter(restOfBEvent);
    RoEListIter.rewind();
    HepAListIterator<BtaCandidate> RoeTrkIter(roeTrkList);
    //RoeTrkIter.rewind();
    //here is an error from VtxGeoKin
    VtxTagBtaSelFit vtxFitRoE(RoEListIter,exclBch,_beamSpot,
			      eventInfo->electronBeam().p3WCov(),
			      eventInfo->positronBeam().p3WCov(), 
			      VtxTagBtaSelFit::OneVtx,KsList);
    vtxFitRoE.setRemoveBeamConstraints();
    vtxFitRoE.setBeamSpotConstraint();
    BtaAbsVertex *vtxRoE = vtxFitRoE.vertex();
    if((vtxRoE != 0) && (vtxRoE->status() == 
			 BtaAbsVertex::Success)&&(fittedB != 0)) {
      //BandBRoeZ = vtxFitRoE.getDz().value(); 
      //BandBRoeZerr = sqrt(vtxFitRoE.getDz().covariance());
      HepSymMatrix RoEVtxCov = vtxRoE->covariance();	
      //error matrix of B vtx point
      HepSymMatrix bVtxCoV = fittedB.decayVtx()->covariance();
      BandBRoeY = - vtxRoE->point().y() + fittedB.decayVtx()->point().y();
      BandBRoeZ = - vtxRoE->point().z() + fittedB.decayVtx()->point().z();
      BandBRoeYerr = sqrt(bVtxCoV(2,2)+RoEVtxCov(2,2));
      BandBRoeZerr = sqrt(bVtxCoV(3,3)+RoEVtxCov(3,3));
    }
    _d0KchCand._BandBRoeY = BandBRoeY;
    _d0KchCand._BandBRoeZ = BandBRoeZ;
    _d0KchCand._BandBRoeYerr = BandBRoeYerr;
    _d0KchCand._BandBRoeZerr = BandBRoeZerr;
    if(isMC()) {
      double dZdf(-100.0);
      if(mcBchrge != -10.0) {
	if(mcBchrge == 0.0) dZdf = BandBRoeZ - mcbDistZ;
	if( mcBchrge*Bchrge == 1 )  dZdf = BandBRoeZ - mcbDistZ;
	if( mcBchrge*Bchrge == -1 ) dZdf = BandBRoeZ + mcbDistZ;
      }
      _d0KchCand._dZdf = dZdf;
    }
			
    // take kaon out of RoE
    double BandBRoEnoKz(-1000.0);
    double BandBRoEnoKzErr(-1000.0);
    HepAList<BtaCandidate> RoEnoKaonTrkList;
    HepAList<BtaCandidate> RoEKaonTrkList;
    RoeTrkIter.rewind();
    tmpCand =  0;
    while(0 != (tmpCand = RoeTrkIter()) ) {
      bool isVTKaon = KaonSMSVt.accept(tmpCand);
      if(isVTKaon) RoEKaonTrkList.append(*tmpCand);	
      if(!isVTKaon) RoEnoKaonTrkList.append(*tmpCand);
      if(!isVTKaon) restOfBnoKEvent.append(*tmpCand);
    }

    //HepAListIterator<BtaCandidate> RoEnoKaonTrkIter(RoEnoKaonTrkList);
    //RoEnoKaonTrkIter.rewind();
    HepAListIterator<BtaCandidate> RoEnoKaonIter(restOfBnoKEvent);
    RoEnoKaonIter.rewind();
    //here is an error from VtxGeoKin
    VtxTagBtaSelFit vtxFitRoEnoK(RoEnoKaonIter,exclBch,
				 VtxTagBtaSelFit::OneVtx);
    vtxFitRoEnoK.setRemoveBeamConstraints();
    vtxFitRoEnoK.setBeamSpotConstraint();
    BtaAbsVertex *vtxRoEnoK = vtxFitRoEnoK.vertex();
    if( (vtxRoEnoK != 0) && (vtxRoEnoK->status() == 
			     BtaAbsVertex::Success)
	&& (fittedB != 0) ) {
      HepSymMatrix RoEnoKVtxCov = vtxRoEnoK->covariance();
      //error matrix of B vtx point
      HepSymMatrix bVtxCoV = fittedB.decayVtx()->covariance();
      BandBRoEnoKz = vtxRoEnoK->point().z()- 
	fittedB.decayVtx()->point().z();
      BandBRoEnoKzErr = sqrt(bVtxCoV(3,3)+RoEnoKVtxCov(3,3));	
    }
    _d0KchCand._BandBRoEnoKz = BandBRoEnoKz;
    _d0KchCand._BandBRoEnoKzErr = BandBRoEnoKzErr;
		
    // find the largest doca from a kaon in ROE wrt signal B
    // vtx and beam spot.
    double KroeBvtxDoca = -1000.0;
    HepAListIterator<BtaCandidate> RoEKaonTrkIter(RoEKaonTrkList);
    RoEKaonTrkIter.rewind();
    if(fittedB != 0 ) {
      tmpCand = 0;
      double tmpDoca = 0.0;
      int inlpflg = 0;			
      while( 0 != (tmpCand = RoEKaonTrkIter()) ) { 
	double tmpChrge = tmpCand->charge();
	if(tmpChrge*exclKch->charge()  == -1.0 ) {	
	  TrkLineTraj Btraj( fittedB.decayVtx()->point(), 
			     exclBRawP4.vect(),-999., 999.);
	  TrkPoca kpoca(*((Trajectory*) &Btraj), 0., 
			*((Trajectory*)&(tmpCand->trkAbsFit()->traj())),0.);
	  double BvtxDoca = kpoca.doca();
	  if(abs(BvtxDoca) > abs(tmpDoca)) {
	    tmpDoca = BvtxDoca;
	    inlpflg = 1;
	  }
	}
      }
      if(inlpflg == 1) KroeBvtxDoca = tmpDoca;
    }
    _d0KchCand._KroeBvtxDoca = KroeBvtxDoca;

    //any non-ks trk of ROE, find the doca and vtx(with beamspot) 
    //distance to B vtx
    double TrkRoeBvtxDoca = -1000.0;
    HepAListIterator<BtaCandidate> RoEnoKsTrkIter(TrkRoeNoKsList);
    RoEnoKsTrkIter.rewind();
    if(fittedB != 0 ) {
      tmpCand = 0;
      double tmpDoca = 0.0;
      int inlp1 = 0;
      HepAList<BtaCandidate>  tmpList;
      while( 0 != (tmpCand = RoEnoKsTrkIter()) ) {
	TrkLineTraj Btraj( fittedB.decayVtx()->point(), 
			   exclBRawP4.vect(),-999., 999.);
	TrkPoca nokspoca(*((Trajectory*) &Btraj), 0.,
			 *((Trajectory*)&(tmpCand->trkAbsFit()->traj())),0.);
	double BvtxDoca = nokspoca.doca();
	if(abs(BvtxDoca)>abs(tmpDoca)) {
	  tmpDoca = BvtxDoca;
	  inlp1 = 1;
	}
      }
      if(inlp1 == 1) TrkRoeBvtxDoca = tmpDoca;
    }
    _d0KchCand._TrkRoeBvtxDoca = TrkRoeBvtxDoca;
    
    // Tagging:
    BtgVariables myTaggingVariables;
    myTaggingVariables.analyze(anEvent, exclBch);
	       
    static const int NTAGGERS = 2;
    
    AbsBTagger * taggers[NTAGGERS] =
      { _electronTagger, _muonTagger };
    
    static bool firstCall = true;
    if (true == firstCall){
      firstCall = false;
      for (int t = 0; t < NTAGGERS; t++){
	taggers[t]->print(cout);
      }
    }
    
    double tagProbs[NTAGGERS] = {-100, -100};
    for (int t = 0; t < NTAGGERS; t++){
      bool tagged = taggers[t]->analyze(anEvent,exclBch);
      if (true == tagged) {
	tagProbs[t] = taggers[t]->probB0();
      }
    }
    
    _d0KchCand._probB0E      = Bchrge*(tagProbs[0]-0.5);
    _d0KchCand._probB0Mu     = Bchrge*(tagProbs[1]-0.5);

    /* Best vtx chisquare of bachelor K with any other trk 
       from the Rest of Event and its doca to beam spot in Y 
       direction      */
    double KbRoETrkChisq = 1000.0;
    double VtxKbRoETrkDistY = -1000.0;
    double VtxKbRoETrkDistChisqY = -1000.0;
    HepAList<BtaCandidate> KbTrkRoEList;
    KbTrkRoEList.append(*exclKch);
    BtaCandidate* RoETrkCand(0);
    HepAListIterator<BtaCandidate> RoETrkListIter(roeTrkList);
    RoETrkListIter.rewind();
    while(0 != (RoETrkCand = RoETrkListIter()) ) {
      KbTrkRoEList.append(RoETrkCand);
      CandVtxFit(KbTrkRoEList,false, false, 0.0);
      if(_fittedCand.decayVtx()->status() == BtaAbsVertex::Success) {
	double ChisqRoETrkKb =  _fittedCand.decayVtx()->chiSquared();
	if(ChisqRoETrkKb < KbRoETrkChisq ) {
	  KbRoETrkChisq = ChisqRoETrkKb;
	  VtxKbRoETrkDistY = _fittedCand.decayVtx()->point().y()
	    - _beamSpot.y();
	  HepSymMatrix KbRoETrkVtxCovMtr = 
	    _fittedCand.decayVtx()->covariance();
	  HepSymMatrix beamSpotCovMtr = _beamSpot.covMatrix();
	  double KbRoETrkVtxCovY = KbRoETrkVtxCovMtr(2,2);
	  double beamSpotCovY = beamSpotCovMtr(2,2);
	  VtxKbRoETrkDistChisqY = (VtxKbRoETrkDistY*VtxKbRoETrkDistY)
	    /(KbRoETrkVtxCovY+beamSpotCovY);
	}
      }
      KbTrkRoEList.remove(RoETrkCand);
    }

    _d0KchCand._KbRoETrkChisq = KbRoETrkChisq;
    _d0KchCand._VtxKbRoETrkDistY = VtxKbRoETrkDistY;
    _d0KchCand._VtxKbRoETrkDistChisqY = VtxKbRoETrkDistChisqY;


    //information for Kch track (charge and pid)
    int HdTrkChg =(int)exclKch->charge();
    _d0KchCand._HdTrkChg = HdTrkChg;	
    
    bool HdTrkPidVt = KaonSMSVt.accept(exclKch);
    bool HdTrkPidTi = KaonSMSTi.accept(exclKch);
    bool HdTrkPidLs = KaonSMSLs.accept(exclKch);
    bool HdTrkPidVl = KaonSMSVl.accept(exclKch);
    bool HdTrkPidNp = KaonSMSNp.accept(exclKch);
    
    int HdTrkPid = HdTrkPidNp+2*HdTrkPidVl+4*HdTrkPidLs+
      8*HdTrkPidTi+16*HdTrkPidVt;
    _d0KchCand._HdTrkPid = HdTrkPid;
    
    //Hard Kch particle
    HepLorentzVector exclKchRawP4;
    exclKchRawP4 = exclKch->p4();

    //hard trk pi+  instead of K+
    HepLorentzVector fakePicP4  = exclKchRawP4;
    // get the hard track momentum info:
    float HdTrkMom = exclKchRawP4.vect().mag();
    float HdTrkTheta = exclKchRawP4.theta();
    float HdTrkCosTheta = exclKchRawP4.cosTheta();
    float HdTrkPhi = exclKchRawP4.phi();
    
    dumptrack(exclKch); 
    //* close approach
    //unused variable
    //float HdTrkPocaZ0 = _pocaz0;
    //float HdTrkPocaD0 = _pocad0;
    //*pid info
    float HdTrkThetaC = _thC; // cheronkov angle thetac
    float HdTrkThcErr = _thCErr;
    double HdTrkKaonNN = _KaonNN;
    double HdTrkNPhot = _NPhot;
           
    //* lepton veto for Hd track:
    static PidElectronMicroSelector ectronMSTi = PidElectronMicroSelector();
    static PidElectronMicroSelector ectronMSVt = PidElectronMicroSelector();
    static bool electronMSIntial = false;
    if( false == electronMSIntial ) {	  
      ectronMSTi.setParmValue("criteria","tight");
      ectronMSVt.setParmValue("criteria","veryTight");
      electronMSIntial = true;
    } 
    bool isElecTight = ectronMSTi.accept(exclKch);
    bool isElecVTight = ectronMSVt.accept(exclKch);
    
    //use the muon pid micro selector
    static PidMuonMicroSelector muonMSTi = PidMuonMicroSelector();
    static PidMuonMicroSelector muonMSVt = PidMuonMicroSelector();
    static bool MuonMSIntializer = false;
    if(false == MuonMSIntializer) {
      muonMSTi.setParmValue("criteria","tight");
      muonMSVt.setParmValue("criteria","veryTight");
      MuonMSIntializer = true;
    }
    bool isMuonTight = muonMSTi.accept(exclKch);
    bool isMuonVTight = muonMSVt.accept(exclKch);
    _d0KchCand._HdTrkMom = HdTrkMom;  
    _d0KchCand._HdTrkTheta = HdTrkTheta;
    _d0KchCand._HdTrkCosTheta = HdTrkCosTheta;
    _d0KchCand._HdTrkPhi      = HdTrkPhi;
    _d0KchCand._HdTrkThetaC   = HdTrkThetaC;
    _d0KchCand._HdTrkThcErr   = HdTrkThcErr;
    _d0KchCand._HdTrkNPhot = HdTrkNPhot;
    _d0KchCand._HdTrkKaonNN = HdTrkKaonNN;

    int ElecPid = 2*isElecVTight + 1*isElecTight;
    _d0KchCand._HdTrkElecPid = ElecPid;
    int MuonPid = 1*isMuonTight+2*isMuonVTight;
    _d0KchCand._HdTrkMuonPid = MuonPid;
    
    //Hard trk CM mometum cut ( herited from skim selector) 
    HepLorentzVector exclKchRawP4CM = exclKchRawP4;
    exclKchRawP4CM.boost(_boostVector);
    const double exclKchCMMomtum = exclKchRawP4CM.vect().mag();
    float HdTrkCMTheta = exclKchRawP4CM.theta(); 
    float HdTrkCMPhi = exclKchRawP4CM.phi();
    _d0KchCand._HdTrkCMMomtum = exclKchCMMomtum;
    _d0KchCand._HdTrkCMTheta = HdTrkCMTheta;
    _d0KchCand._HdTrkCMPhi = HdTrkCMPhi;

    
    // go through each D0 decay mod:
    // D0 decay type: 1, D0 -> K  K  Pi0;   2, D0bar/D0 -> K+  pi- K_S0;
    //                3, D0bar/D0 -> K- pi+ K_S0; 4, D0 -> pi+ pi- pi0; 
    int d0kchDecMod = 0;

    BtaCandidate * Kplus = 0;
    BtaCandidate * Kminus = 0;
    BtaCandidate * ks = 0;
    BtaCandidate * Piplus = 0;
    BtaCandidate * Piminus = 0;
    BtaCandidate * pi0 = 0;
    
    //updated daughters --JZ
    BtaCandidate * KplusUP = 0;
    BtaCandidate * KminusUP = 0;
    BtaCandidate * ksUP = 0;
    BtaCandidate * PiplusUP = 0;
    BtaCandidate * PiminusUP = 0;
    BtaCandidate * pi0UP = 0;
    
    if(d0kchDecMod == 0 ) {	
      Kplus=0; Kminus=0; pi0=0;
      KplusUP=0; KminusUP=0; pi0UP=0;
      findD0Descendents(exclD0,Kplus, PdtLund::K_plus, Kminus, 
			PdtLund::K_minus,pi0,PdtLund::pi0);
      if (fittedMSD0 != 0){
	findD0Descendents(const_cast<const BtaCandidate*>(&fittedMSD0), 
			  KplusUP, PdtLund::K_plus, KminusUP, 
			  PdtLund::K_minus,pi0UP,PdtLund::pi0);
      }
      if((Kplus != 0)&&(Kminus !=0)&&(pi0 !=0)) {
	d0kchDecMod = 1;
      }
    }
			
    if(d0kchDecMod == 0) {	
      Kplus=0; Piminus=0; ks=0;
      KplusUP=0; PiminusUP=0; ksUP=0;
      findD0Descendents(exclD0,Kplus,PdtLund::K_plus,Piminus,
			PdtLund::pi_minus, ks, PdtLund::K_S0);
      if (fittedMSD0 !=0){
	findD0Descendents(const_cast<const BtaCandidate*>(&fittedMSD0),
			  KplusUP,PdtLund::K_plus,PiminusUP,
			  PdtLund::pi_minus, ksUP, PdtLund::K_S0);
      }
      if((Kplus != 0)&&(Piminus !=0)&&(ks !=0)) {
	int D0Id = exclD0->pdtEntry()->lundId();
	if(D0Id == PdtLund::D0) d0kchDecMod = 2;
	if(D0Id == PdtLund::anti_D0) d0kchDecMod = 3;
      }
    }
	 
    if(d0kchDecMod == 0) {
      Kminus=0; Piplus=0; ks=0;
      KminusUP=0; PiplusUP=0; ksUP=0;
      findD0Descendents(exclD0,Kminus,PdtLund::K_minus,Piplus,
			PdtLund::pi_plus,ks,PdtLund::K_S0);
      if (fittedMSD0 != 0){
	findD0Descendents(const_cast<const BtaCandidate*>(&fittedMSD0),
			  KminusUP,PdtLund::K_minus,PiplusUP,
			  PdtLund::pi_plus,ksUP,PdtLund::K_S0);
      }
      if((Kminus != 0)&&(Piplus != 0)&&(ks != 0)) {
	int D0Id = exclD0->pdtEntry()->lundId();
	if(D0Id == PdtLund::D0) d0kchDecMod = 4;
	if(D0Id == PdtLund::anti_D0) d0kchDecMod = 5;
      }
    }
			
    if(d0kchDecMod == 0) {
      Piplus=0; Piminus=0;pi0=0;
      PiplusUP=0; PiminusUP=0;pi0UP=0;
      findD0Descendents(exclD0,Piplus,PdtLund::pi_plus,Piminus,
			PdtLund::pi_minus,pi0,PdtLund::pi0);

      if (fittedMSD0 != 0){
	findD0Descendents(const_cast<const BtaCandidate*>(&fittedMSD0),
			  PiplusUP,PdtLund::pi_plus,PiminusUP,
			  PdtLund::pi_minus,pi0UP,PdtLund::pi0);
      }
      if((Piminus != 0)&&(Piplus != 0)&&(pi0 != 0)) {
	d0kchDecMod = 6;
      }
    }	

    if(d0kchDecMod == 0 ) continue;
    _d0KchCand._D0DecMode = d0kchDecMod;

    /*  maximum and minium invariant mass of bachelor Kaon 
	with other charge kaons or neutral kaons            */ 

    static PidElectronMicroSelector TrkElecMSVT = PidElectronMicroSelector();
    static PidMuonMicroSelector TrkMuonMSVT = PidMuonMicroSelector();
    TrkElecMSVT.setParmValue("criteria","veryTight");
    TrkMuonMSVT.setParmValue("criteria","veryTight");
    
    double KbLepMass = 0.0;
    int nKaonp = -1;
    int nKaonm = -1;  
    double KbKLowMass = 100.0; // which number??
    double KbKUpMass = 0.0;  // which number??

    //all chrge trks    
    BtaCandidate * KbTrkCand(0); 
    HepAListIterator<BtaCandidate> KbIterTrkList(*btaTrkList);        
    KbIterTrkList.rewind(); //track list 
    while( 0 != (KbTrkCand = KbIterTrkList()) ) {
      if(KbTrkCand->overlaps(*exclBch)) continue; 
      HepLorentzVector TrkP4 = KbTrkCand->p4();
      int trkchrge = (int)KbTrkCand->charge();
      bool isVTKaon = KaonSMSVt.accept(KbTrkCand);
      if(isVTKaon && trkchrge == -1 ) nKaonm += 1;
      if(isVTKaon && trkchrge == 1 ) nKaonp += 1;
      //kaon mass
      if( (isVTKaon) && (trkchrge*HdTrkChg == -1) ) {
	//mass
	TrkP4.setVectM(TrkP4.vect(),_massKch);
	HepLorentzVector KbKMix4Vect = TrkP4 + exclKchRawP4;
	double KbKMixMass = KbKMix4Vect.mag();
	if(KbKMixMass>KbKUpMass) {
	  KbKUpMass = KbKMixMass;
	}
	if(KbKMixMass<KbKLowMass) {
	  KbKLowMass = KbKMixMass;
	}	 
      }

      //lepton mass
      bool VTLepM = TrkMuonMSVT.accept(KbTrkCand);
      if(VTLepM &&(trkchrge*HdTrkChg == -1) ) {			
	TrkP4.setVectM(TrkP4.vect(),_massMuon); //add muon mass
	HepLorentzVector KbLepM4Vect = TrkP4 + exclKchRawP4;
	double KbLepMuMass = KbLepM4Vect.mag();
	if(KbLepMuMass>KbLepMass) {
	  KbLepMass = KbLepMuMass;
	}
      }
      bool VTLepE = TrkElecMSVT.accept(KbTrkCand);
      if(VTLepE &&(trkchrge*HdTrkChg == -1) ) {
	TrkP4.setVectM(TrkP4.vect(),_massElec); //add muon mass
	HepLorentzVector KbLepE4Vect = TrkP4 + exclKchRawP4;
	double KbLepEcMass = KbLepE4Vect.mag();
	if(KbLepEcMass>KbLepMass) {
	  KbLepMass = KbLepEcMass;
	}
      }  
    }
    _d0KchCand._KbKUpMass = KbKUpMass;	
    _d0KchCand._KbKLowMass = KbKLowMass;
    _d0KchCand._KbLepMass = KbLepMass;	
    _d0KchCand._nRoEKp = nKaonp;
    _d0KchCand._nRoEKm = nKaonm;	
    
    int nRoEKs = -1;
    double KbKsLowMass = 100.0;
    double KbKsUpMass = 0.0;
    if(KsRoEList.length() > 0 ) {
      nRoEKs = KsRoEList.length();
      HepAListIterator<BtaCandidate> KsRoEIter(KsRoEList);
      KsRoEIter.rewind();
      BtaCandidate * KsCand(0);
      while( 0 != (KsCand = KsRoEIter()) ) {
	HepLorentzVector Ks4Vect = KsCand->p4();
	HepLorentzVector KbKsMix4Vect = Ks4Vect + exclKchRawP4;
	double KbKsMixMass = KbKsMix4Vect.mag();
	if(KbKsMixMass>KbKsUpMass) {
	  KbKsUpMass = KbKsMixMass;
	}
	if(KbKsMixMass<KbKsLowMass) {
	  KbKsLowMass = KbKsMixMass;
	}
      }
    }		
    _d0KchCand._nRoEKs = nRoEKs;
    _d0KchCand._KbKsLowMass = KbKsLowMass;
    _d0KchCand._KbKsUpMass = KbKsUpMass;
    
    //for Pi0-> 2 gamma daughters:
    BtaCandidate * pi0Daughter1 = 0;
    BtaCandidate * pi0Daughter2 = 0;
    BtaCandidate * mcPi0Dau1 = 0;
    BtaCandidate * mcPi0Dau2 = 0;
    HepLorentzVector pi0UPRawP4;
    if( pi0UP != 0 ) {
      HepLorentzVector pi0P4;
      pi0Daughter1 = 0;
      pi0Daughter2 = 0;
      HepAListIterator<BtaCandidate> pi0DaughterIter =
	pi0UP->daughterIterator();
      pi0DaughterIter.rewind();
      // Total pi0 momentum before the mass constraint:
      BtaCandidate * pi0Daughter(0);
      while (pi0Daughter = pi0DaughterIter()){
	if (0 == pi0Daughter1) {
	  pi0Daughter1 = pi0Daughter;
	} else {
	  pi0Daughter2 = pi0Daughter;
	}
	// Accumulate pre-constraint p4 of pi0:
	pi0P4 += pi0Daughter->p4();
      }
      
      // check pi0 energy if mc
      if(isMC()) {
	double truPi0E1 = -10;
	double truPi0E2 = -10;
	mcPi0Dau1 = 0;
	mcPi0Dau2 = 0;
	mcPi0Dau1 = _truthMap->mcFromReco(pi0Daughter1);
	mcPi0Dau2 = _truthMap->mcFromReco(pi0Daughter2);
	if( 0 != mcPi0Dau1 ) {
	  if( 22 == mcPi0Dau1->pdtEntry()->lundId()) {
	    truPi0E1 = mcPi0Dau1->energy();
	  }
	}
	if( 0 != mcPi0Dau2 ) {
	  if( 22 == mcPi0Dau2->pdtEntry()->lundId() ) {
	    truPi0E2 = mcPi0Dau2->energy();
	  }
	}
	_d0KchCand._truPi0E1 = truPi0E1;
	_d0KchCand._truPi0E2 = truPi0E2;
      }
      
      pi0UPRawP4 = pi0UP->p4();
      HepLorentzVector pi0CMP4 = pi0UPRawP4;
      double pi0CMmom = pi0CMP4.vect().mag();
      pi0CMP4.boost(_boostVector);
      double pi0CMmom0 = pi0CMP4.vect().mag();
      
      dumpGamma(pi0Daughter1);
      double pi0gamE1 = _GamEn;
      double pi0gamLat1(-1.0);
      pi0gamLat1 =  _latMom;
      //photon 1 energy cut:
      if(pi0gamE1<0.03) continue;
      //get the cm energy for photon
      HepLorentzVector ph1P4 = pi0Daughter1->p4();
      ph1P4.boost(_boostVector);
      
      dumpGamma(pi0Daughter2);
      double pi0gamE2 = _GamEn;
      double pi0gamLat2(-1.0);
      pi0gamLat2 = _latMom;
      // photon 2 energy cut:
      if(pi0gamE2<0.03) continue;
      HepLorentzVector ph2P4 = pi0Daughter2->p4();
      ph2P4.boost(_boostVector);
      
      //pi0 helicity angle:
      //BtaHelicity pi0Hel = BtaHelicity(pi0Daughter1, pi0Daughter2);
      //double pi0Helicity = pi0Hel.helicity();
      //double cospi0Hel = cos(pi0Helicity);
      double cospi0Hel = calcHelicity(pi0Daughter1->p4(), 
				      pi0UPRawP4, exclBRawP4);
      double cospi0HelDfm = calcHelicity(pi0Daughter1->p4(), 
					 pi0UPRawP4, D0FitP4);
      double pi0EAsy = (ph1P4.t()-ph2P4.t())/(ph1P4.t()+ph2P4.t());
      
      //go through the gamma list to find the best mass(pi0)
      double pi0BestMass1 = -10.0;
      double pi0BestMass2 = -10.0;
      double Bestcospi0Hel1 = -10.0;
      double Bestcospi0Hel2 = -10.0;
      double massflg1=100.0;
      double massflg2=100.0;
      double Best2gPcm1 = -100;
      double Best2gPcm2 = -100;
      HepAListIterator<BtaCandidate> itergammaList(*gammaList);
      itergammaList.rewind();
      BtaCandidate * gamma(0);
      while( 0 != (gamma = itergammaList()) ) {
	if( !gamma->overlaps(*pi0UP) ) {
	  if( ph1P4.vect().mag()<ph2P4.vect().mag() ) {
	    BtaCandidate * pionSwp(0);
	    double gamEnSwp(0.0);
	    gamEnSwp = pi0gamE1;
	    pi0gamE1 = pi0gamE2;
	    pi0gamE2 = gamEnSwp;		
	    pionSwp = pi0Daughter1;
	    pi0Daughter1=pi0Daughter2;
	    pi0Daughter2=pionSwp;
	  }		 	
	  HepLorentzVector Bestpi01P4 = gamma->p4()+pi0Daughter1->p4();
	  HepLorentzVector Bestpi02P4 = gamma->p4()+pi0Daughter2->p4();
	  
	  double Bestpi01mass = Bestpi01P4.mag();
	  double Bestpi02mass = Bestpi02P4.mag();
	  if(abs(Bestpi01mass-0.135)<massflg1 ) {
	    massflg1 = abs(Bestpi01mass-0.135);
	    pi0BestMass1 = Bestpi01mass;	
	    //best pi0 helic
	    Bestcospi0Hel1 = calcHelicity(gamma->p4(),Bestpi01P4,
					  eventInfo->cmFrame());
	    HepLorentzVector Btpi01P4= Bestpi01P4;
	    Btpi01P4.boost(_boostVector);
	    Best2gPcm1 = Btpi01P4.vect().mag();
	  }
	  if(abs(Bestpi02mass-0.135)<massflg2 ) {
	    massflg2 = abs(Bestpi02mass-0.135);
	    pi0BestMass2 = Bestpi02mass;
	    //best pi0 helic
	    Bestcospi0Hel2 = calcHelicity(gamma->p4(),Bestpi02P4,
					  eventInfo->cmFrame());
	    HepLorentzVector Btpi02P4= Bestpi02P4;
	    Btpi02P4.boost(_boostVector);
	    Best2gPcm2 = Btpi02P4.vect().mag();	
	  }	 	
	}
      }
      _d0KchCand._pi0BestMass1 = pi0BestMass1;
      _d0KchCand._pi0BestMass2 = pi0BestMass2;
      _d0KchCand._Bestcospi0Hel1 = Bestcospi0Hel1;
      _d0KchCand._Bestcospi0Hel2 = Bestcospi0Hel2;
      _d0KchCand._Best2gPcm1 = Best2gPcm1;
      _d0KchCand._Best2gPcm2 = Best2gPcm2;
      
      const double pi0InvMass = pi0P4.mag();	            
      //pi0 mass cut added here
      if(false == inRange( (pi0InvMass-_massPi0), 
			   _exclpi0MassMin.value(), 
			   _exclpi0MassMax.value() )) {
	continue;
      }
              
      _d0KchCand._pi0Mass = pi0InvMass;
      _d0KchCand._gamEng1 = pi0gamE1;
      _d0KchCand._gamEng2 = pi0gamE2;
      _d0KchCand._gamLat1 = pi0gamLat1;
      _d0KchCand._gamLat2 = pi0gamLat2;
      _d0KchCand._pi0Hel = cospi0Hel;
      _d0KchCand._pi0HelDfm = cospi0HelDfm;
      _d0KchCand._pi0CMmom = pi0CMmom;
      _d0KchCand._pi0CMmom0 = pi0CMmom0;
      _d0KchCand._pi0EAsy = pi0EAsy;
      
      //pi0 variables required by Abi
      BtaCandidate * tpi0 = pi0UP;
      Pi0Vars vars(*tpi0, *pi0List, *eventInfo, gammaList, _truthMap);
      _d0KchCand._mom = vars.mom();
      _d0KchCand._mass = vars.mass();
      _d0KchCand._helic = vars.helic();
      _d0KchCand._good = vars.good();
      
      for (int i = 0; i < vars.NPAIRS; i++) {
	_d0KchCand._massH[i] = vars.massH(i);
	_d0KchCand._massS[i] = vars.massS(i);
	_d0KchCand._helicH[i] = vars.helicH(i);
	_d0KchCand._helicS[i] = vars.helicS(i);
	_d0KchCand._asymH[i] = vars.asymH(i);
	_d0KchCand._asymS[i] = vars.asymS(i);
	_d0KchCand._momH[i] = vars.momH(i);
	_d0KchCand._momS[i] = vars.momS(i);
	_d0KchCand._goodH[i] = vars.goodH(i);
	_d0KchCand._goodS[i] = vars.goodS(i);
      }
    }

    // without D0 mass constraint
    HepLorentzVector pi0RawP4;
    if( pi0 != 0 ) {
      pi0RawP4 = pi0->p4();
    }

    //information for K_S0 -> pi+ pi-
    HepLorentzVector KsUPRawP4;
    if( ksUP != 0 ) {
      // the mass calculated by the two ks daughter is not right
      // should use the vtx fit mass!!!!!!
      //HepLorentzVector KsP4(0,0,0,0);
      HepAListIterator<BtaCandidate> KsDaughterIter =
	ksUP->daughterIterator(); 
      KsDaughterIter.rewind();
      // Total Ks momentum before the mass constraint:
      BtaCandidate * KsDaughter1(0);
      BtaCandidate * KsDaughter(0);
      while (0 != (KsDaughter = KsDaughterIter())){
	// Accumulate pre-constraint p4 of Ks:
	if(0 == KsDaughter1 ) KsDaughter1 = KsDaughter;
      }
      double KsHel = calcHelicity(KsDaughter1->p4(),ksUP->p4(),D0FitP4); 
      
      //ks already vtx fitted in tcl file:   
      double KsVtxFitChisq =  ksUP->decayVtx()->chiSquared();
      _d0KchCand._KsVtxFitChisq = KsVtxFitChisq;
      
      BtaCandidate* tmpKs(0);
      double KsInvMass = -10.0;
      CompBaseAnalUtil tmpUtil;
      if(0 != ( tmpKs= const_cast<BtaCandidate*>
		(tmpUtil.getOriginalCandidate( ksUP,KsLsList)))) {
	KsInvMass = tmpKs->p4().mag();
      } 	
      
      //Ks decay veterx point:
      const HepPoint KsVtx = ksUP->decayVtx()->point();	
      //ks and the momentum of ks. 
      double KsTheta(-100);
      double pKs(-10000); 
      KsUPRawP4 = ksUP->p4();
      KsTheta = KsUPRawP4.theta();
      Hep3Vector Ks3Vect = KsUPRawP4.vect();
      pKs = Ks3Vect.mag();
      //error matrix of ks decay vtx point
      HepSymMatrix KsVtxCoV = ksUP->decayVtx()->covariance();

      double Kslth(-100.);
      double D0KsVtxVectChisq(-100.);
      double cosDKsVtx(-100.);
      if(fittedD0 != 0) {
	//D0 decay veterx point:
	const HepPoint D0Vtx = fittedD0.decayVtx()->point();	  
	Hep3Vector D0Ks3Vect = D0Vtx - KsVtx;
	Kslth  = D0Ks3Vect.mag();
	//Ks decay length 
	//cosine of the angle between the vector which links 
	//vertex of D0 and 
	cosDKsVtx = cos(D0Ks3Vect.angle(Ks3Vect));		
	
	// to find the chisquar for the decay length
	//D0 Vtx point error matrix
	HepSymMatrix D0VtxCoV = fittedD0.decayVtx()->covariance(); 

	BbrVectorErr D0KsVtxVectErr(D0Ks3Vect);
	D0KsVtxVectErr.setCovMatrix(KsVtxCoV + D0VtxCoV);
	D0KsVtxVectChisq = D0KsVtxVectErr.determineChisq(Hep3Vector());
      }
      if(false == inRange((KsInvMass-_massKs), 
			  _exclKsMassMin.value(),_exclKsMassMax.value() )) {
	continue;
      }
      _d0KchCand._cosDKsVtx = (float)cosDKsVtx;
      _d0KchCand._dksVtxVctChisq = D0KsVtxVectChisq;	
      _d0KchCand._KsDecLen = Kslth;
      _d0KchCand._KsTheta = KsTheta;
      _d0KchCand._pKs     = pKs;
      _d0KchCand._KsMass   = KsInvMass;
      _d0KchCand._KsHel    = KsHel;
    }

    // without D0 mass constraint	
    HepLorentzVector KsRawP4;
    if( ks != 0 ) {
      KsRawP4 = ks->p4();
    }

    //information on D0
    HepLorentzVector KplusRawP4;
    HepLorentzVector KminusRawP4;
    HepLorentzVector PiplusRawP4;
    HepLorentzVector PiminusRawP4;
    HepLorentzVector D0fakeP4;
    // with D0 mass constraint
    HepLorentzVector KplusUPRawP4;
    HepLorentzVector KminusUPRawP4;
    HepLorentzVector PiplusUPRawP4;
    HepLorentzVector PiminusUPRawP4;
    HepLorentzVector D0UPRawP4;

    int KpPidBit = -1;
    int KmPidBit = -1; 
    int PipPidBit = -1;
    int PimPidBit = -1;
   
    if(d0kchDecMod == 1) { //D0 -> K K pi0
      //kaon:
      KplusRawP4 = Kplus->p4();
      KminusRawP4 = Kminus->p4();
      
      //D0 4 momentum
      D0RawP4 = KplusRawP4 + KminusRawP4 + pi0RawP4;
      
      //Trk info
      //k+
      dumptrack(Kplus);
      double KpPocaz0 = _pocaz0;
      double KpPocad0	= _pocad0;
      dumptrack(Kminus);
      double KmPocaz0 = _pocaz0;
      double KmPocad0 = _pocad0;
      
      double KpTheta = KplusRawP4.theta();
      double KmTheta = KminusRawP4.theta();
      double KpPhi   = KplusRawP4.phi();
      double KmPhi   = KminusRawP4.theta();
      
      //volume cut here	
      _d0KchCand._D0DauTrkP1 = KplusRawP4.vect().mag();
      _d0KchCand._D0DauTrkP2 = KminusRawP4.vect().mag();
      _d0KchCand._D0DauTrkTheta1 = KpTheta;
      _d0KchCand._D0DauTrkTheta2 = KmTheta;
      _d0KchCand._D0DauTrkPhi1   = KpPhi;
      _d0KchCand._D0DauTrkPhi2   = KmPhi;
      
      if(abs(KpPocaz0)>10||abs(KpPocad0)>1.5) continue;
      if(abs(KmPocaz0)>10||abs(KmPocad0)>1.5) continue;
      
      //Pid info for Kaon comes from D0:
      //1.notApion
      bool KnotApiAcc1 = KaonSMSNp.accept(Kplus);
      bool KnotApiAcc2 = KaonSMSNp.accept(Kminus);
      
      //2.veryloose
      bool KvLAcc1 = KaonSMSVl.accept(Kplus);
      bool KvLAcc2 = KaonSMSVl.accept(Kminus);
      
      //3. Loose
      bool KlsAcc1 = KaonSMSLs.accept(Kplus);
      bool KlsAcc2 = KaonSMSLs.accept(Kminus);
      
      //4. Tight
      bool KtiAcc1 = KaonSMSTi.accept(Kplus);
      bool KtiAcc2 = KaonSMSTi.accept(Kminus);
      
      //5. veryTight
      bool KvTAcc1 = KaonSMSVt.accept(Kplus);
      bool KvTAcc2 = KaonSMSVt.accept(Kminus);
      
      //packed Kaon pid info:
      KpPidBit = 1*KvLAcc1+2*KlsAcc1+4*KtiAcc1+8*KvTAcc1+16*KnotApiAcc1;
      KmPidBit = 1*KvLAcc2+2*KlsAcc2+4*KtiAcc2+8*KvTAcc2+16*KnotApiAcc2;
      
      //saved the pid info:
      _d0KchCand._KpPidBit = KpPidBit;
      _d0KchCand._KmPidBit = KmPidBit;
      
      // two fake d0 trk k<->pi  
      HepLorentzVector fakeD0dTrkPipP4; //K+ replaced by pi+
      HepLorentzVector fakeD0dTrkPimP4; //K- replaced by pi-
      
      //=== two body compound mass ========       
      // K+K- (phi??) compound mass, helicity: no D0 mass constraint
      HepLorentzVector KKcompRawP4;
      
      KKcompRawP4 = KplusRawP4 + KminusRawP4;
      const double Kch2Mass = KKcompRawP4.mag();
      const double Kch2Helic = calcHelicity(KplusRawP4, 
					    KKcompRawP4,D0RawP4);	     
      
      //K+ pi0 (K*+) compound mass, helicity: no D0 mass constraint
      HepLorentzVector KpPi0compRawP4;
      
      KpPi0compRawP4 = KplusRawP4 + pi0RawP4;
      const double KpPi0Mass = KpPi0compRawP4.mag();
      const double KpPi0Helic = calcHelicity(KplusRawP4,
					     KpPi0compRawP4,D0RawP4);
      
      //K- pi0 (K*-) compound mass, helicity: no D0 mass constraint
      HepLorentzVector KmPi0compRawP4;
      KmPi0compRawP4 = KminusRawP4 + pi0RawP4;
      const double KmPi0Mass = KmPi0compRawP4.mag();
      const double KmPi0Helic = calcHelicity(KminusRawP4,
					     KmPi0compRawP4,D0RawP4);			    
      if ( fittedMSD0 != 0) {
	// with D0 mass constraint
	KplusUPRawP4 = KplusUP->p4();
	KminusUPRawP4 = KminusUP->p4();
	
	//D0 4 momentum
	D0UPRawP4 = KplusUPRawP4 + KminusUPRawP4 + pi0UPRawP4;
	
	// K+K- (phi??) compound mass, helicity: D0 mass constraint
	HepLorentzVector KKUPcompRawP4;		
	KKUPcompRawP4 = KplusUPRawP4 + KminusUPRawP4;
	const double Kch2UPMass = KKUPcompRawP4.mag();
	const double Kch2UPHelic = calcHelicity(KplusUPRawP4,
						KKUPcompRawP4,D0UPRawP4);

	//K+ pi0 (K*+) compound mass, helicity: D0 mass constraint
	HepLorentzVector KpPi0UPcompRawP4;
	KpPi0UPcompRawP4 = KplusUPRawP4 + pi0UPRawP4;
	const double KpPi0UPMass = KpPi0UPcompRawP4.mag();
	const double KpPi0UPHelic = calcHelicity(KplusUPRawP4,
						 KpPi0UPcompRawP4,D0UPRawP4);

	//K- pi0 (K*-) compound mass, helicity: D0 mass constraint
	HepLorentzVector KmPi0UPcompRawP4;
	KmPi0UPcompRawP4 = KminusUPRawP4 + pi0UPRawP4;
	const double KmPi0UPMass = KmPi0UPcompRawP4.mag();
	const double KmPi0UPHelic = calcHelicity(KminusUPRawP4,
						 KmPi0UPcompRawP4,D0UPRawP4);
	
	//more with D0 mass constraint
	_d0KchCand._d0kkUPMass  = Kch2UPMass;  // K+ K- mass
	_d0KchCand._d0pKpUPMass = KpPi0UPMass; // K+ pi0 mass
	_d0KchCand._d0pKmUPMass = KmPi0UPMass; // K- pi0 mass
	_d0KchCand._d0kkUPHelic  = Kch2UPHelic;  //K+K-, pi0
	_d0KchCand._d0pKpUPHelic = KpPi0UPHelic; //K+pi0, K-
	_d0KchCand._d0pKmUPHelic = KmPi0UPHelic; //K-pi0, K+
      }

      //===Prompt Track mix with neutral or opposite chrge track===      
      double MixNeuMass = 0.0;// mass( K+ K- ) high momentum!!!
      if(HdTrkChg == 1) {
	HepLorentzVector MixNeuRawP4 = exclKchRawP4 + KminusRawP4;
	MixNeuMass = MixNeuRawP4.mag();
      }
      if(HdTrkChg == -1 ) {
	HepLorentzVector MixNeuRawP4 = exclKchRawP4 + KplusRawP4;
	MixNeuMass = MixNeuRawP4.mag();	
      } 
      
      double MixChgMass = 0.0; //mass(K+/- pi0) K*+/- 
      HepLorentzVector MixChgRawP4 = exclKchRawP4 + pi0RawP4;
      MixChgMass = MixChgRawP4.mag();
      
      //======= fake mass ================== 
      //fake D0 mass Kaon is replaced by pion
      
      fakeD0dTrkPipP4 = KplusRawP4;
      fakeD0dTrkPimP4 = KminusRawP4;
      
      // set up the wrong kaon track with correct particle mass 		
      fakeD0dTrkPipP4.setVectM(KplusRawP4.vect(), _massPi);
      fakeD0dTrkPimP4.setVectM(KminusRawP4.vect(), _massPi);
      
      //K_plus replaced by pi_plus:	
      
      D0fakeP4 = fakeD0dTrkPipP4+KminusRawP4+pi0RawP4;
      const double D0fakeMass2Kppi0 = D0fakeP4.mag();
      //fake helicity angle 1
      
      //K_minus replaced by pi_minus:			
      D0fakeP4 = fakeD0dTrkPimP4+KminusRawP4+pi0RawP4;
      const double D0fakeMass2Kmpi0 = D0fakeP4.mag();
      
      //two kaon trks are swaped by pion
      HepLorentzVector D0swpP4 = fakeD0dTrkPimP4 + 
	fakeD0dTrkPipP4 + pi0RawP4;
      const double Md0Swp = D0swpP4.mag();
 
      //==== dump parameters=======
      //fake D0 mass:
      _d0KchCand._d0FkMpKp = D0fakeMass2Kppi0; // K+ -> pi+
      _d0KchCand._d0FkMpKm = D0fakeMass2Kmpi0; // K- -> pi-
      _d0KchCand._d0SwpMass = Md0Swp; //K <-> pi
      
      // tow body compound mass:
      _d0KchCand._d0kkMass  = Kch2Mass;  // K+ K- mass
      _d0KchCand._d0pKpMass = KpPi0Mass; // K+ pi0 mass
      _d0KchCand._d0pKmMass = KmPi0Mass; // K- pi0 mass
      
      // hard trk compound mass:
      _d0KchCand._MixNeuMass = MixNeuMass;
      _d0KchCand._MixChgMass = MixChgMass;
      
      // tow body compound Helicity:
      _d0KchCand._d0kkHelic  = Kch2Helic;  //K+K-, pi0
      _d0KchCand._d0pKpHelic = KpPi0Helic; //K+pi0, K-
      _d0KchCand._d0pKmHelic = KmPi0Helic; //K-pi0, K+
    }

    if(d0kchDecMod == 2 || d0kchDecMod == 3 ) { 
      //D0 ->K+ pi- K_S0 or D0bar -> K+ pi- K_S0
      KplusRawP4 = Kplus->p4(); //K_plus
      
      //trk info
      dumptrack(Kplus);
      double KpPocaz0 = _pocaz0;
      double KpPocad0 = _pocad0;
      double KpTheta = KplusRawP4.theta();
      double KpPhi   = KplusRawP4.phi();
      
      //kaon trk volume cut here
      
      _d0KchCand._D0DauTrkP1 = KplusRawP4.vect().mag();
      _d0KchCand._D0DauTrkTheta1 = KpTheta;
      _d0KchCand._D0DauTrkPhi1   = KpPhi;
      if(abs(KpPocaz0)>10||abs(KpPocad0)>1.5) continue;
      
      //Pid info for this kaon
      //1. not pion
      bool KnotApiAcc1 = KaonSMSNp.accept(Kplus);
      //2. veryLoose
      bool KvLAcc1 = KaonSMSVl.accept(Kplus);
      //3. Loose
      bool KlsAcc1 = KaonSMSLs.accept(Kplus);
      //4. tight
      bool KtiAcc1 = KaonSMSTi.accept(Kplus);
      //5. veryTight
      bool KvTAcc1 = KaonSMSVt.accept(Kplus);
      
      //packed Pid info:
      KpPidBit = 1*KvLAcc1+2*KlsAcc1+4*KtiAcc1+
	8*KvTAcc1+16*KnotApiAcc1;
      
      //save the pid info:
      _d0KchCand._KpPidBit = KpPidBit;
      
      // got the pi- momentum             
      PiminusRawP4 = Piminus->p4(); //pi_minus
      
      //D0 four momentum:
      D0RawP4 = KsRawP4 + PiminusRawP4 + KplusRawP4;
      
      //pion trk info
      dumptrack(Piminus);
      double PimPocaz0 = _pocaz0;
      double PimPocad0 = _pocad0;
      
      double PimTheta = PiminusRawP4.theta();
      double PimPhi   = PiminusRawP4.phi();
      //pion trk volume cut here
      
      _d0KchCand._D0DauTrkP2 = PiminusRawP4.vect().mag();
      _d0KchCand._D0DauTrkTheta2 = PimTheta;
      _d0KchCand._D0DauTrkPhi2    = PimPhi;
      
      if(abs(PimPocaz0)>10||abs(PimPocad0)>1.5) continue;
      
      //pion pid info:(pi-)
      //1. veryLoose;
      bool PvLAcc1 = pionLHSVl.accept(Piminus);
      //2. loose
      bool PlsAcc1 = pionLHSLs.accept(Piminus);
      //3. tight
      bool PtiAcc1 = pionLHSTi.accept(Piminus);
      //4. veryTight
      bool PvTAcc1 = pionLHSVt.accept(Piminus);
      
      //packed pion pid info:		
      PimPidBit = 1*PvLAcc1+2*PlsAcc1+4*PtiAcc1+8*PvTAcc1;
      _d0KchCand._PimPidBit = PimPidBit;
      
      //use kaon method to identify pion:
      //1. not pion
      bool KnotApiAcc2 = KaonSMSNp.accept(Piminus);
      //2. veryLoose
      bool KvLAcc2 = KaonSMSVl.accept(Piminus);
      //3. Loose
      bool KlsAcc2 = KaonSMSLs.accept(Piminus);
      //4. tight
      bool KtiAcc2 = KaonSMSTi.accept(Piminus);
      //5. veryTight
      bool KvTAcc2 = KaonSMSVt.accept(Piminus);

      //packed Pid info:
      KmPidBit = 1*KvLAcc2+2*KlsAcc2+4*KtiAcc2+
	8*KvTAcc2+16*KnotApiAcc2;
      
      //save the pid info for the pion(kaon):
      _d0KchCand._KmPidBit = KmPidBit;
      
      //Piminus is replaced by exclKch(minus), here used to 
      //rejected the background
      //B- -> D0(K+ K- ks) pi- 	
      if(HdTrkChg == -1 ) {
	HepLorentzVector D0FkPromtP4 = KplusRawP4 + exclKchRawP4 + KsRawP4;
	_d0KchCand._D0FkPromptMass = D0FkPromtP4.mag();
      } 		
		
      //======= tow body compound mass ===============
      //K+ Pi- (K*0) tow body compound mass  and helicity angle
      // no D0 mass constraint
      HepLorentzVector KpPicompRawP4;
      
      KpPicompRawP4 = KplusRawP4 + PiminusRawP4;
      const double KpPi2Mass = KpPicompRawP4.mag();
      const double KpPi2Helic = calcHelicity(KplusRawP4,
					     KpPicompRawP4,D0RawP4);

      //K+ K_S0 (????) compound mass and helicity
      // no D0 mass constraint
      HepLorentzVector KpKscompRawP4;
      
      KpKscompRawP4 = KplusRawP4 + KsRawP4;
      const double KpKsMass = KpKscompRawP4.mag();
      const double KpKsHelic = calcHelicity(KplusRawP4,
					    KpKscompRawP4,D0RawP4);

      //pi- K_S0 (K*-) compound mass and helicity
      // no D0 mass constraint
      HepLorentzVector PimKscompRawP4;
      PimKscompRawP4 = PiminusRawP4 + KsRawP4;
      const double PimKsMass = PimKscompRawP4.mag();
      const double PimKsHelic = calcHelicity(PiminusRawP4,
					     PimKscompRawP4,D0RawP4);  
      
      if ( fittedMSD0 != 0) {
	//with D0 mass constrain
	KplusUPRawP4 = KplusUP->p4(); //K_plus
	
	//with D0 mass constraint
	PiminusUPRawP4 = PiminusUP->p4(); //pi_minus
	//with D0 mass constraint
	D0UPRawP4 = KsUPRawP4 + PiminusUPRawP4 + KplusUPRawP4;
	
	//K+ Pi- (K*0) tow body compound mass  and helicity angle
	// D0 mass constraint
	HepLorentzVector KpPiUPcompRawP4;
	
	KpPiUPcompRawP4 = KplusUPRawP4 + PiminusUPRawP4;
	const double KpPi2UPMass = KpPiUPcompRawP4.mag();
	const double KpPi2UPHelic = calcHelicity(KplusUPRawP4,
						 KpPiUPcompRawP4,D0UPRawP4);

	//K+ K_S0 (????) compound mass and helicity
	// D0 mass constraint
	HepLorentzVector KpKsUPcompRawP4;
	
	KpKsUPcompRawP4 = KplusUPRawP4 + KsUPRawP4;
	const double KpKsUPMass = KpKsUPcompRawP4.mag();
	const double KpKsUPHelic = calcHelicity(KplusUPRawP4,
						KpKsUPcompRawP4,D0UPRawP4);

	//pi- K_S0 (K*-) compound mass and helicity
	// D0 mass constraint
	HepLorentzVector PimKsUPcompRawP4;
	PimKsUPcompRawP4 = PiminusUPRawP4 + KsUPRawP4;
	const double PimKsUPMass = PimKsUPcompRawP4.mag();
	const double PimKsUPHelic = calcHelicity(PiminusUPRawP4,
						 PimKsUPcompRawP4,D0UPRawP4);

	//more D0 mass constraint
	_d0KchCand._d0pKpUPMass = KpPi2UPMass; // K+ pi-
	_d0KchCand._d0kkUPMass  = KpKsUPMass;  // K+ Ks
	_d0KchCand._d0pKmUPMass = PimKsUPMass; // pi- Ks
	_d0KchCand._d0pKpUPHelic = KpPi2UPHelic; //K+pi-, K_S0
	_d0KchCand._d0kkUPHelic  = KpKsUPHelic; //K+K_S0, pi-
	_d0KchCand._d0pKmUPHelic = PimKsUPHelic; //pi-Ks, K+
      }
      
      //===Prompt Track mix with neutral or opposite chrge track===
      double MixNeuMass = 0.0;// mass( K+/pi- K-/+ ) high momentum!!!
      if(HdTrkChg == 1) {// K+ pi- (K*0)
	HepLorentzVector MixNeuRawP4 = exclKchRawP4 + PiminusRawP4;
	MixNeuMass = MixNeuRawP4.mag();
      }
      if(HdTrkChg == -1 ) {// K+ K- (phi???)
	HepLorentzVector MixNeuRawP4 = exclKchRawP4 + KplusRawP4;
	MixNeuMass = MixNeuRawP4.mag();
      }
      
      double MixChgMass = 0.0; //mass(K+/- Ks) Ds+??
      HepLorentzVector MixChgRawP4 = exclKchRawP4 + KsRawP4;
      MixChgMass = MixChgRawP4.mag();
      
      //=====fake two body mass ========
      
      // wrong identified D0 daugter trk, K+ -> pi+/ pi- ->K-
      HepLorentzVector fakeD0dTrkPipP4 = KplusRawP4;	
      HepLorentzVector fakeD0dTrkKmP4 = PiminusRawP4;
      
      // set up the correct mass for wrong pion(Kaon) trk:
      fakeD0dTrkPipP4.setVectM(fakeD0dTrkPipP4.vect(), _massPi);
      fakeD0dTrkKmP4.setVectM(fakeD0dTrkKmP4.vect(),_massKch);
      
      // K_plus replaced by pi_plus:
      D0fakeP4 = fakeD0dTrkPipP4 + PiminusRawP4 + KsRawP4;
      HepLorentzVector D0SwpP4 = fakeD0dTrkPipP4 + fakeD0dTrkKmP4 + KsRawP4;
      
      //fake D0 mass
      const double D0fakeMass2PipKs = D0fakeP4.mag();
      const double Md0Swp = D0SwpP4.mag();
      //fake helicity angle
      
      // dump parameters:
      //fake D0 mass:
      _d0KchCand._d0FkMpKp = D0fakeMass2PipKs; // K+ -> pi+
      _d0KchCand._d0SwpMass = Md0Swp;
      
      // tow body compound mass:
      _d0KchCand._d0pKpMass = KpPi2Mass; // K+ pi-
      _d0KchCand._d0kkMass  = KpKsMass;  // K+ Ks
      _d0KchCand._d0pKmMass = PimKsMass; // pi- Ks
      
      //HdTrk mix two body mass
      _d0KchCand._MixNeuMass = MixNeuMass;
      _d0KchCand._MixChgMass = MixChgMass;     
      
      // tow body compound Helicity:
      _d0KchCand._d0pKpHelic = KpPi2Helic; //K+pi-, K_S0
      _d0KchCand._d0kkHelic  = KpKsHelic; //K+K_S0, pi-
      _d0KchCand._d0pKmHelic = PimKsHelic; //pi-Ks, K+
      
    }
    
    if(d0kchDecMod == 4 || d0kchDecMod == 5) { //D0 -> K- pi+ K_S0
      KminusRawP4 = Kminus->p4(); //K_minus
            
      //trk info:
      dumptrack(Kminus);
      double KmPocaz0 = _pocaz0;
      double KmPocad0 = _pocad0;
      double KmTheta = KminusRawP4.theta();
      double KmPhi   = KminusRawP4.phi();
      //volume cut here
      _d0KchCand._D0DauTrkP1 = KminusRawP4.vect().mag();
      _d0KchCand._D0DauTrkTheta1 = KmTheta; 
      _d0KchCand._D0DauTrkPhi1   = KmPhi;
      if(abs(KmPocaz0)>10||abs(KmPocad0)>1.5) continue;
      
      //Pid info for this kaon
      //PidKaonSMSSelector();
      //1. not pion
      bool KnotApiAcc1 = KaonSMSNp.accept(Kminus);
      //2. veryLoose
      bool KvLAcc1 = KaonSMSVl.accept(Kminus);
      //3. Loose
      bool KlsAcc1 = KaonSMSLs.accept(Kminus);
      //4. tight
      bool KtiAcc1 = KaonSMSTi.accept(Kminus);
      //5. veryTight
      bool KvTAcc1 = KaonSMSVt.accept(Kminus);
      
      //packed Kminus pid info:
      KmPidBit = 1*KvLAcc1+2*KlsAcc1+4*KtiAcc1+
	8*KvTAcc1+16*KnotApiAcc1;
      _d0KchCand._KmPidBit = KmPidBit;
      
      // get the pi track momentum
      PiplusRawP4 = Piplus->p4(); //pi + 
      
      //D0 4 momentum:
      D0RawP4 = PiplusRawP4 + KminusRawP4 + KsRawP4;
      
      // pi track info:
      dumptrack(Piplus);
      double PipPocaz0 = _pocaz0;
      double PipPocad0 = _pocad0;
      double PipTheta = PiplusRawP4.theta();
      double PipPhi   = PiplusRawP4.phi();
      // volume cut here
      _d0KchCand._D0DauTrkP2 = PiplusRawP4.vect().mag();
      _d0KchCand._D0DauTrkTheta2 = PipTheta;
      _d0KchCand._D0DauTrkPhi2   = PipPhi;	
      if(abs(PipPocaz0)>10||abs(PipPocad0)>1.5) continue;				
      //pion pid info
      //pionLHS->setParmValue("veryLoose");		
      bool PvLAcc2 = pionLHSVl.accept(Piplus);		
      //2. loose
      bool PlsAcc2 = pionLHSLs.accept(Piplus);
      //3. tight
      bool PtiAcc2 = pionLHSTi.accept(Piplus);
      //4. veryTight
      bool PvTAcc2 = pionLHSVt.accept(Piplus);
      
      //packed pion pid info:
      PipPidBit = 1*PvLAcc2+2*PlsAcc2+4*PtiAcc2+8*PvTAcc2;
      _d0KchCand._PipPidBit = PipPidBit;			
      
      //use kaon method to identify this pion+
      //1. not pion
      bool KnotApiAcc2 = KaonSMSNp.accept(Piplus);
      //2. veryLoose
      bool KvLAcc2 = KaonSMSVl.accept(Piplus);
      //3. Loose
      bool KlsAcc2 = KaonSMSLs.accept(Piplus);
      //4. tight
      bool KtiAcc2 = KaonSMSTi.accept(Piplus);
      //5. veryTight
      bool KvTAcc2 = KaonSMSVt.accept(Piplus);

      //packed Pid info:
      KpPidBit = 1*KvLAcc2+2*KlsAcc2+4*KtiAcc2+
	8*KvTAcc2+16*KnotApiAcc2;

      //save the pid info:
      _d0KchCand._KpPidBit = KpPidBit;
      
      //Piplus is replaced by exclKch(plus), here used to rejected 
      //the background
      //B- -> D0(K+ K- ks) pi+
      if(HdTrkChg == 1 ) {
	HepLorentzVector D0FkPromtP4 = KminusRawP4 + exclKchRawP4 + KsRawP4;
	_d0KchCand._D0FkPromptMass =(float)D0FkPromtP4.mag();
      }
   
      // ========tow body compound mass===================
      
      //K- Pi+ (K*0) tow body compound mass  and helicity angle
      // no D0 mass constraint
      HepLorentzVector KmPicompRawP4;
      
      KmPicompRawP4 = KminusRawP4 + PiplusRawP4;
      const double KmPi2Mass = KmPicompRawP4.mag();
      const double KmPi2Helic = calcHelicity(KminusRawP4,
					     KmPicompRawP4,D0RawP4);
   
      //K- K_S0 (K*-) compound mass and helicity
      // no D0 mass constraint
      HepLorentzVector KmKscompRawP4;
      
      KmKscompRawP4 = KminusRawP4 + KsRawP4;	
      const double KmKsMass = KmKscompRawP4.mag();
      const double KmKsHelic = calcHelicity(KminusRawP4,
					    KmKscompRawP4,D0RawP4);

      //pi+ K_S0 (K*+) compound mass and helicity
      // no D0 mass constraint
      HepLorentzVector PipKscompRawP4;
      PipKscompRawP4 = PiplusRawP4 + KsRawP4; 
      const double PipKsMass = PipKscompRawP4.mag();
      const double PipKsHelic = calcHelicity(PiplusRawP4,
					     PipKscompRawP4,D0RawP4);
      
      if (fittedMSD0 != 0) {
	//Do mass constraint
	KminusUPRawP4 = KminusUP->p4();
	
	// D0 mass constraint
	PiplusUPRawP4 = PiplusUP->p4(); //pi +
	
	// D0 mass constraint
	D0UPRawP4 = PiplusUPRawP4 + KminusUPRawP4 + KsUPRawP4;
	
	//K- Pi+ (K*0) tow body compound mass  and helicity angle
	// D0 mass constraint
	HepLorentzVector KmPiUPcompRawP4;
	
	KmPiUPcompRawP4 = KminusUPRawP4 + PiplusUPRawP4;
	const double KmPi2UPMass = KmPiUPcompRawP4.mag();
	const double KmPi2UPHelic = calcHelicity(KminusUPRawP4,
						 KmPiUPcompRawP4,D0UPRawP4);

	//K- K_S0 (K*-) compound mass and helicity
	// D0 mass constraint
	HepLorentzVector KmKsUPcompRawP4;
	
	KmKsUPcompRawP4 = KminusUPRawP4 + KsUPRawP4;	
	const double KmKsUPMass = KmKsUPcompRawP4.mag();
	const double KmKsUPHelic = calcHelicity(KminusUPRawP4,
						KmKsUPcompRawP4,D0UPRawP4);

	//pi+ K_S0 (K*+) compound mass and helicity
	// D0 mass constraint
	HepLorentzVector PipKsUPcompRawP4;
	PipKsUPcompRawP4 = PiplusUPRawP4 + KsUPRawP4; 
	const double PipKsUPMass = PipKsUPcompRawP4.mag();
	const double PipKsUPHelic = calcHelicity(PiplusUPRawP4,
						 PipKsUPcompRawP4,D0UPRawP4);
	//more D0 mass constraint
	_d0KchCand._d0pKmUPMass = KmPi2UPMass; //K- pi+
	_d0KchCand._d0kkUPMass  = KmKsUPMass;  //K- Ks
	_d0KchCand._d0pKpUPMass = PipKsUPMass; //pi+ Ks
	_d0KchCand._d0pKmUPHelic = KmPi2UPHelic; //K-pi+, K_S0
	_d0KchCand._d0kkUPHelic  = KmKsUPHelic ; //K-K_S0, pi+
	_d0KchCand._d0pKpUPHelic  = PipKsUPHelic; //pi+K_S0,  K-	
      }

      //===Prompt Track mix with neutral or opposite chrge track===
      
      double MixNeuMass = 0.0;// mass( K-/pi+ K+/- ) high momentum!!!
      if(HdTrkChg == 1) {
	HepLorentzVector MixNeuRawP4 = exclKchRawP4 + KminusRawP4;
	MixNeuMass = MixNeuRawP4.mag();
      }
      if(HdTrkChg == -1 ) {
	HepLorentzVector MixNeuRawP4 = exclKchRawP4 + PiplusRawP4;
	MixNeuMass = MixNeuRawP4.mag();
      }
      
      double MixChgMass = 0.0; //mass(K+/- Ks) Ds+/-
      HepLorentzVector MixChgRawP4 = exclKchRawP4 + KsRawP4;
      MixChgMass = MixChgRawP4.mag();   
      
      //========== fake D0 mass kaon replaced by pion====
      
      //wrong identified d0 daughter trk K- -> pi-:
      HepLorentzVector fakeD0dTrkPimP4 = KminusRawP4; 
      HepLorentzVector fakeD0dTrkKpP4 = PiplusRawP4;
      
      fakeD0dTrkPimP4.setVectM(fakeD0dTrkPimP4,_massPi);
      fakeD0dTrkKpP4.setVectM(fakeD0dTrkKpP4,_massKch);
      
      // K_minus replaced by pi_minus:
      D0fakeP4 = fakeD0dTrkPimP4 + PiplusRawP4 + KsRawP4;
      HepLorentzVector D0SwpP4 = fakeD0dTrkPimP4 + fakeD0dTrkKpP4 + KsRawP4;
      
      //fake D0 mass
      const double D0fakeMass2PimKs = D0fakeP4.mag();
      const double Md0Swp = D0SwpP4.mag();
      
      //fake helicity angle
      
      //dump parameters:
      //fake D0 Mass
      _d0KchCand._d0FkMpKm = D0fakeMass2PimKs; // K- -> pi-
      _d0KchCand._d0SwpMass = Md0Swp;
      
      // tow body compound mass:
      _d0KchCand._d0pKmMass = KmPi2Mass; //K- pi+
      _d0KchCand._d0kkMass  = KmKsMass;  //K- Ks
      _d0KchCand._d0pKpMass = PipKsMass; //pi+ Ks

      //HdTrk mixed neutral and opposite charge trk
      _d0KchCand._MixNeuMass = MixNeuMass;
      _d0KchCand._MixChgMass = MixChgMass; 
                
      // two body compound helicity:
      _d0KchCand._d0pKmHelic = KmPi2Helic; //K-pi+, K_S0
      _d0KchCand._d0kkHelic  = KmKsHelic ; //K-K_S0, pi+
      _d0KchCand._d0pKpHelic  = PipKsHelic; //pi+K_S0,  K-
    }
    
    if(d0kchDecMod == 6) { //D0 -> (pi- pi+) pi0
      PiminusRawP4 = Piminus->p4(); //pi_plus
      PiplusRawP4 = Piplus->p4();   //pi_plus
      
      //D0 momentum:
      D0RawP4 = PiminusRawP4 + PiplusRawP4 + pi0RawP4;
      
      //trk info:
      dumptrack(Piminus);
      double PimPocaz0 = _pocaz0;
      double PimPocad0 = _pocad0;
      double PimTheta = PiminusRawP4.theta();
      double PimPhi   = PiminusRawP4.phi();
      
      dumptrack(Piplus);
      double PipPocaz0 = _pocaz0;
      double PipPocad0 = _pocad0;
      double PipTheta = PiplusRawP4.theta();
      double PipPhi   = PiplusRawP4.phi();
      
      //volume cut here:
      _d0KchCand._D0DauTrkP1 = PiminusRawP4.vect().mag();
      _d0KchCand._D0DauTrkP2 = PiplusRawP4.vect().mag();
      
      _d0KchCand._D0DauTrkTheta1 = PipTheta;
      _d0KchCand._D0DauTrkTheta2 = PimTheta; 
      
      _d0KchCand._D0DauTrkPhi1   = PipPhi;
      _d0KchCand._D0DauTrkPhi2   = PimPhi;
      
      if(abs(PimPocaz0)>10||abs(PimPocad0)>1.5) continue;
      if(abs(PipPocaz0)>10||abs(PipPocad0)>1.5) continue;
      
      // set up the pion pid info
      //1. veryLoose
      bool PvLAcc1 = pionLHSVl.accept(Piminus);
      bool PvLAcc2 = pionLHSVl.accept(Piplus);		
      //2. loose
      bool PlsAcc1 = pionLHSLs.accept(Piminus);
      bool PlsAcc2 = pionLHSLs.accept(Piplus);
      //3. tight
      bool PtiAcc1 = pionLHSTi.accept(Piminus);
      bool PtiAcc2 = pionLHSTi.accept(Piplus);
      //4. veryTight
      bool PvTAcc1 = pionLHSVt.accept(Piminus);
      bool PvTAcc2 = pionLHSVt.accept(Piplus);
      
      //packed pion pid info:
      PipPidBit = 1*PvLAcc2+2*PlsAcc2+4*PtiAcc2+8*PvTAcc2;
      PimPidBit = 1*PvLAcc1+2*PlsAcc1+4*PtiAcc1+8*PvTAcc1;
      _d0KchCand._PipPidBit = PipPidBit;
      _d0KchCand._PimPidBit = PimPidBit;
      
      //use the kaon info to identify pions:
      //1.notApion
      bool KnotApiAcc1 = KaonSMSNp.accept(Piplus);
      bool KnotApiAcc2 = KaonSMSNp.accept(Piminus);
      
      //2.veryloose
      bool KvLAcc1 = KaonSMSVl.accept(Piplus);
      bool KvLAcc2 = KaonSMSVl.accept(Piminus);

      //3. Loose
      bool KlsAcc1 = KaonSMSLs.accept(Piplus);
      bool KlsAcc2 = KaonSMSLs.accept(Piminus);
      
      //4. Tight
      bool KtiAcc1 = KaonSMSTi.accept(Piplus);
      bool KtiAcc2 = KaonSMSTi.accept(Piminus);
      
      //5. veryTight
      bool KvTAcc1 = KaonSMSVt.accept(Piplus);
      bool KvTAcc2 = KaonSMSVt.accept(Piminus);
      
      //packed Kaon pid info:
      KpPidBit = 1*KvLAcc1+2*KlsAcc1+4*KtiAcc1+
	8*KvTAcc1+16*KnotApiAcc1;
      KmPidBit = 1*KvLAcc2+2*KlsAcc2+4*KtiAcc2+
	8*KvTAcc2+16*KnotApiAcc2;
      
      //saved the pid info:
      _d0KchCand._KpPidBit = KpPidBit;
      _d0KchCand._KmPidBit = KmPidBit;
      
      
      //Piminus is replaced by exclKch(minus), here used to 
      //rejected the background
      //B- -> D0(pi+ K- pi0) pi-
      if(HdTrkChg == -1 ) {
	HepLorentzVector D0FkPromtP4 = PiplusRawP4 + exclKchRawP4 + pi0RawP4;
	_d0KchCand._D0FkPromptMass = D0FkPromtP4.mag();
      }
      
      //Piplus is replaced by exclKch(plus), here used to 
      //rejected the background
      //B- -> D0(K+ pi- pi0) pi+
      if(HdTrkChg == 1 ) {
	HepLorentzVector D0FkPromtP4 = PiminusRawP4 + 
	  exclKchRawP4 + pi0RawP4;
	_d0KchCand._D0FkPromptMass = D0FkPromtP4.mag();
      }

      double pipiDecLen = -100.0;
      double pipiMass = 0; 
      tmpCand = 0;
      int pipiKsTrue = 0; 
      KsIterator.rewind();
      while(0 != (tmpCand = KsIterator()) ) {
	if(tmpCand->overlaps(*Piminus) && tmpCand->overlaps(*Piplus)) {
	  pipiMass = tmpCand->p4().mag();
	  if( fittedB != 0 )  {
	    Hep3Vector pipiBVect = fittedB.decayVtx()->point() 
	      - tmpCand->decayVtx()->point();
	    pipiDecLen = pipiBVect.mag();
	  }
	  if( isMC() ) {
	    BtaCandidate * TrueTmpKs(0);
	    BtaCandidate * tmpPipls(0);
	    BtaCandidate * tmpPimns(0);
	    TrueTmpKs = _truthMap->mcFromReco(tmpCand );
	    tmpPipls = _truthMap->mcFromReco(Piplus);
	    tmpPimns = _truthMap->mcFromReco(Piminus);
	    if(TrueTmpKs !=0 && tmpPipls != 0 && tmpPimns != 0) {
	      if(TrueTmpKs->overlaps(*tmpPipls)&&
		 TrueTmpKs->overlaps(*tmpPimns)) 
		pipiKsTrue = 1;
	    }
	  }	
	}
      } 
      _d0KchCand._pipiDecLen = pipiDecLen;
      _d0KchCand._pipiMass = pipiMass; 
      _d0KchCand._pipiKsTrue = pipiKsTrue;

      // ======two body compound mass====================      
      //pi+ pi-(Ks,rho0) compound mass and helicity
      //no D0 mass constraint
      HepLorentzVector PiPiRawP4;
      PiPiRawP4 = PiminusRawP4 + PiplusRawP4;
      const double PiPiMass = PiPiRawP4.mag();
      const double PiPiHelic = calcHelicity(PiplusRawP4,
					    PiPiRawP4,D0RawP4);
      
      
      //pi+ pi0(rho+) compound mass and helicity
      //no D0 mass constraint
      HepLorentzVector PipPi0RawP4;
      
      PipPi0RawP4 = PiplusRawP4 + pi0RawP4;
      const double PipPi0Mass = PipPi0RawP4.mag(); 	
      const double PipPi0Helic = calcHelicity(PiplusRawP4,
					      PipPi0RawP4, D0RawP4);    
      
      //pi- pi0 (rho-) compound mass and helicity
      // no D0 mass constraint
      HepLorentzVector PimPi0RawP4;
      PimPi0RawP4 = PiminusRawP4 + pi0RawP4;
      const double PimPi0Mass = PimPi0RawP4.mag();
      const double PimPi0Helic = calcHelicity(PiminusRawP4,
					      PimPi0RawP4,D0RawP4);

      if ( fittedMSD0 != 0) {
	//D0 mass constraint
	PiminusUPRawP4 = PiminusUP->p4(); //pi_plus
	PiplusUPRawP4 = PiplusUP->p4();   //pi_plus
	// D0 mass constraint
	D0UPRawP4 = PiminusUPRawP4 + PiplusUPRawP4 + pi0UPRawP4;
	
	//pi+ pi-(Ks,rho0) compound mass and helicity
	//D0 mass constraint
	HepLorentzVector PiPiUPRawP4;
	PiPiUPRawP4 = PiminusUPRawP4 + PiplusUPRawP4;
	const double PiPiUPMass = PiPiUPRawP4.mag();
	const double PiPiUPHelic = calcHelicity(PiplusUPRawP4,
						PiPiUPRawP4,D0UPRawP4);
	//pi+ pi0(rho+) compound mass and helicity
	// D0 mass constraint
	HepLorentzVector PipPi0UPRawP4;
	
	PipPi0UPRawP4 = PiplusUPRawP4 + pi0UPRawP4;
	const double PipPi0UPMass = PipPi0UPRawP4.mag(); 	
	const double PipPi0UPHelic = calcHelicity(PiplusUPRawP4,
						  PipPi0UPRawP4, D0UPRawP4);
	//pi- pi0 (rho-) compound mass and helicity
	// D0 mass constraint
	HepLorentzVector PimPi0UPRawP4;
	PimPi0UPRawP4 = PiminusUPRawP4 + pi0UPRawP4;
	const double PimPi0UPMass = PimPi0UPRawP4.mag();
	const double PimPi0UPHelic = calcHelicity(PiminusUPRawP4,
						  PimPi0UPRawP4, D0UPRawP4);
	
	//more D0 mass constraint
	_d0KchCand._d0ppUPMass  = PiPiUPMass;  //pi+ pi-
	_d0KchCand._d0pPpUPMass = PipPi0UPMass; //pi+ pi0
	_d0KchCand._d0pPmUPMass = PimPi0UPMass; //pi- pi0
	_d0KchCand._d0ppUPHelic  = PiPiUPHelic;  // pi+pi-, pi0
	_d0KchCand._d0pPpUPHelic = PipPi0UPHelic;  // pi+pi0, pi-
	_d0KchCand._d0pPmUPHelic = PimPi0UPHelic;  // pi-pi0, pi+
      }

      //===Prompt Track mix with neutral or opposite chrge track===
      HepLorentzVector D0WrongP4;
     
      double MixNeuMass = 0.0;// mass( K+/- pi-/+ ) high momentum!!!
      if(HdTrkChg == 1) { //K+ pi- (K*0)
	HepLorentzVector MixNeuRawP4 = exclKchRawP4 + PiminusRawP4;
	MixNeuMass = MixNeuRawP4.mag();
	D0WrongP4 = exclKchRawP4 + PiminusRawP4 + pi0RawP4;
      }
      if(HdTrkChg == -1 ) {//K- pi+ (K*0)
	HepLorentzVector MixNeuRawP4 = exclKchRawP4 + PiplusRawP4;
	MixNeuMass = MixNeuRawP4.mag();
	D0WrongP4 = exclKchRawP4 + PiplusRawP4 + pi0RawP4;
      }
      
      double MixChgMass = 0.0; //mass(K+/- pi0) K*+/-
      HepLorentzVector MixChgRawP4 = exclKchRawP4 + pi0RawP4;
      MixChgMass = MixChgRawP4.mag();
      
      double D0WrongMass = D0WrongP4.mag();
      //====== fake D0 mass pion replaced by kaon====
      
      //wrong identified channel: pi <-> k
      HepLorentzVector fakeD0dTrkKpP4;
      HepLorentzVector fakeD0dTrkKmP4;
      HepLorentzVector D0SwpP4;
      
      fakeD0dTrkKpP4 = PiplusRawP4;
      fakeD0dTrkKmP4 = PiminusRawP4;
      fakeD0dTrkKpP4.setVectM(fakeD0dTrkKpP4.vect(), _massKch);
      fakeD0dTrkKmP4.setVectM(fakeD0dTrkKmP4.vect(), _massKch);
      
      //fake D0 mass 
      D0fakeP4 = fakeD0dTrkKmP4 + PiplusRawP4 + pi0RawP4;
      const double D0fakeMass2PiKm = D0fakeP4.mag();
      D0fakeP4 = PiminusRawP4 + fakeD0dTrkKpP4 + pi0RawP4;
      const double D0fakeMass2PiKp = D0fakeP4.mag();
      D0SwpP4 = fakeD0dTrkKmP4 + fakeD0dTrkKpP4 + pi0RawP4;
      const double Md0Swp = D0SwpP4.mag();
      
      // dump parameters: 	
      // fake D0 Mass:
      _d0KchCand._d0FkMpKp = D0fakeMass2PiKp; // pi+ -> K+
      _d0KchCand._d0FkMpKm = D0fakeMass2PiKm; // pi- -> K-
      _d0KchCand._d0SwpMass = Md0Swp;
      
      // two body compound mass:
      _d0KchCand._d0ppMass  = PiPiMass;  //pi+ pi-
      _d0KchCand._d0pPpMass = PipPi0Mass; //pi+ pi0
      _d0KchCand._d0pPmMass = PimPi0Mass; //pi- pi0

      //Hdtrk mixed;
      _d0KchCand._MixNeuMass = MixNeuMass;
      _d0KchCand._MixChgMass = MixChgMass; 
      _d0KchCand._D0WrongMass = D0WrongMass;

      // tow body compound helicity:
      _d0KchCand._d0ppHelic  = PiPiHelic;  // pi+pi-, pi0
      _d0KchCand._d0pPpHelic = PipPi0Helic;  // pi+pi0, pi-
      _d0KchCand._d0pPmHelic = PimPi0Helic;  // pi-pi0, pi+
    }

    //to tell the reconstruct particle mc truth
    if( isMC() ) {
      
      //to tell pi+ pi- mc truth
      int pipTru = -10;
      int pimTru = -10;
      int trueKsPP = -1;
      int trueDToKsX = -1;
      BtaCandidate * truePipTmp = 0;
      BtaCandidate * truePimTmp = 0;
      BtaCandidate* ppGGMother =0;
      if(0 != PiplusUP || 0 != PiminusUP ) {
	BtaCandidate * aD1moth(0);
	BtaCandidate * aD2moth(0);
	int aD1mothId = -10;
	int aD2mothId = -10;
	
	if(0 != PiplusUP ) {
	  truePipTmp = _truthMap->mcFromReco(PiplusUP);
	  int truePipTmpId = -10; 
	  if(0 != truePipTmp) truePipTmpId = 
				abs(truePipTmp->pdtEntry()->lundId());
	  if(PdtLund::pi_plus == truePipTmpId  ) {
	    aD1moth = truePipTmp->theMother();
	    if(aD1moth != 0) aD1mothId = abs(aD1moth->pdtEntry()->lundId());
	    if( PdtLund::D0 == aD1mothId ) {
	      if(0 != BdecModMC[0]) pipTru = 1;
	      if(0 != BdecModMC[1]) pipTru = 11;
	      if(0 != BdecModMC[1] && 0 != BdecModMC[0]) pipTru = 21;
	    } 
	  }
	}  				 
	  
	if(0 != PiminusUP) {
	  truePimTmp = _truthMap->mcFromReco(PiminusUP);
	  int truePimTmpId = -10;
	  if(0 != truePimTmp) truePimTmpId = 
				abs(truePimTmp->pdtEntry()->lundId());
	  if(PdtLund::pi_plus == truePimTmpId  ) {
	    aD2moth = truePimTmp->theMother();
	    if(aD2moth != 0) aD2mothId = abs(aD2moth->pdtEntry()->lundId());
	    if( PdtLund::D0 == aD2mothId ) {
	      if(0 != BdecModMC[0]) pimTru = 1;
	      if(0 != BdecModMC[1]) pimTru = 11;
	      if(0 != BdecModMC[1] && 0 != BdecModMC[0]) pimTru = 21;
	    }
	  }
	}
	
	// check if Ks-> pi+ pi-
	if(aD1moth != 0 && aD2moth!= 0 && abs(aD2mothId)==PdtLund::K_S0
	   && aD1moth->overlaps((const BtaCandidate&)(*aD2moth))) {
	  trueKsPP = 1;

	  int GMotherId;
	  int GGMotherId;
	  BtaCandidate* aGMother = aD1moth->theMother();
	  if( 0 != aGMother) {			
	    GMotherId = aGMother->pdtEntry()->lundId();
	    ppGGMother = aGMother->theMother();
	    if( 0 != ppGGMother) {
	      GGMotherId = ppGGMother->pdtEntry()->lundId();
	    }
	  }

	  if(abs(GMotherId) == PdtLund::K0 && 
	     abs(GGMotherId) == PdtLund::D0) {
	    trueDToKsX = 1;
	  }	 
	}
	_d0KchCand._pipTru = pipTru ;
	_d0KchCand._pimTru = pimTru ;
	_d0KchCand._d0trueKsPP = trueKsPP;
      }

      // tell the D0
      int OrigD0(0);
      int D0RecDec(0);
      if(exclD0 != 0 ) {
	BtaCandidate * OrigOfTD0(0);
	OrigOfTD0 = _truthMap->mcFromReco(exclD0);
	
	if(OrigOfTD0 != 0 ){
	  BtaCandidate * MothOfD0(0);
	  MothOfD0 = OrigOfTD0->theMother();
	  if(MothOfD0 != 0 ) {
	    int MD0Id=0;
	    MD0Id = MothOfD0->pdtEntry()->lundId();
	    OrigD0 = -20;
	    if(abs(MD0Id)%1000>410&&abs(MD0Id)%1000<416) {
	      OrigD0 = -4;//D* charge
	    }
	    if(abs(MD0Id)%1000>420&&abs(MD0Id)%1000<426){ 
	      OrigD0 = -6;//D*0 
	    }
	    if(abs(MD0Id) == PdtLund::B0||abs(MD0Id) == PdtLund::B_plus) {
	      OrigD0 = -8; //B+/B0
	    }
	  }	

	  int D0TruId(0);
	  
	  D0TruId = OrigOfTD0->pdtEntry()->lundId();
	  if(abs(D0TruId) == PdtLund::D0) { 
	    if(OrigD0 == 0) OrigD0 = 0;
	    //positve, good; .le.0, bad
	    if(OrigD0 != 0) OrigD0 = 10 + OrigD0 ; 

	    bool D0ToKKpi0=findD0Descendents(OrigOfTD0,PdtLund::K_plus,
					     PdtLund::K_minus, PdtLund::pi0);
	    if(D0ToKKpi0) D0RecDec = 1;
	    bool D0ToKKks=findD0Descendents(OrigOfTD0,PdtLund::K_plus,
					    PdtLund::K_minus, PdtLund::K_S0);
	    if(D0ToKKks) D0RecDec = 2;
	    bool D0ToPPpi0=findD0Descendents(OrigOfTD0, PdtLund::pi_plus,
					     PdtLund::pi_minus,PdtLund::pi0);
	    if(D0ToPPpi0) D0RecDec = 3;
	    bool D0ToPPKs=findD0Descendents(OrigOfTD0, PdtLund::pi_plus, 
					    PdtLund::pi_minus,PdtLund::K_S0);
	    if(D0ToPPKs) D0RecDec = 4;
	    bool D0ToKpPpi0=findD0Descendents(OrigOfTD0,PdtLund::K_plus, 
					      PdtLund::pi_minus,PdtLund::pi0);
	    if(D0ToKpPpi0) D0RecDec = 5;
	    bool D0ToKpPks=findD0Descendents(OrigOfTD0,PdtLund::K_plus, 
					     PdtLund::pi_minus,PdtLund::K_S0);
	    if(D0ToKpPks) D0RecDec = 6;
	    bool D0ToKmPpi0=findD0Descendents(OrigOfTD0,PdtLund::pi_plus,  
					      PdtLund::K_minus,PdtLund::pi0); 
	    if(D0ToKmPpi0) D0RecDec = 7;
	    bool D0ToKmPks=findD0Descendents(OrigOfTD0, PdtLund::pi_plus, 
					     PdtLund::K_minus, PdtLund::K_S0);
	    if(D0ToKmPks) D0RecDec = 8;
	    bool D0Tokspi0=findD0Descendents(OrigOfTD0, PdtLund::K_S0,
					     PdtLund::pi0);
	    if(D0Tokspi0) D0RecDec = 9;

	  }	
	}
      }
      _d0KchCand._OrigD0 = OrigD0;
      _d0KchCand._d0RecDec = D0RecDec;
    
      int Kchflg(0);
      if( exclKch != 0 ) {
	BtaCandidate * OrigOfKch(0);
	OrigOfKch = _truthMap->mcFromReco(exclKch);
	int tmpId(0);
	if(OrigOfKch != 0) {
	  //1.K+, 2. pi+ 3. e+ 4. mu+   
	  tmpId = OrigOfKch->pdtEntry()->lundId();
	  //cout<<"Hard Trk: "<<OrigOfKch->pdtEntry()->name()<<endl;
	  if(PdtLund::K_plus == tmpId ) Kchflg=1;
	  if(PdtLund::K_minus == tmpId ) Kchflg=2;
	  if(PdtLund::pi_plus == tmpId ) Kchflg = 3;
	  if(PdtLund::pi_minus == tmpId ) Kchflg = 4;
	  if(PdtLund::e_minus == tmpId ) Kchflg = 5;
	  if(PdtLund::e_plus == tmpId ) Kchflg = 6;
	  if(PdtLund::mu_minus == tmpId ) Kchflg = 7;
	  if(PdtLund::mu_plus == tmpId ) Kchflg = 8;
	  BtaCandidate * mothOfKch(0);
	  mothOfKch = OrigOfKch->theMother();
	  if(mothOfKch != 0) {
	    int MothOfKchId = mothOfKch->pdtEntry()->lundId();
	    if(abs(MothOfKchId)==PdtLund::K_star_plus) 
	      Kchflg = 10 + Kchflg;
	    if(abs(MothOfKchId)==PdtLund::K_star0) 
	      Kchflg = 20 + Kchflg;
	    if(abs(MothOfKchId)==PdtLund::D0||
	       abs(MothOfKchId)==PdtLund::D_plus) Kchflg = 30 + Kchflg;
	    if(abs(MothOfKchId)==PdtLund::D_star0||
	       abs(MothOfKchId)== PdtLund::D_star_plus) Kchflg = 40 + Kchflg;
	    if(abs(MothOfKchId)==PdtLund::B0||
	       abs(MothOfKchId)== PdtLund::B_plus) Kchflg = 50 + Kchflg;
	  }
	  
	  const BtaCandidate * KaonAncest(0);
	  KaonAncest = ancestor(OrigOfKch);
	  if(KaonAncest!=0) {
	    HepAListIterator<BtaCandidate> 
	      kAdaugIter(KaonAncest->daughterIterator());
	    int totChrm = 0;
	    int totNrm = 0;
	    int tmpNrm = 0;
	    int tmpChrm = 0;
	    // unused variable
	    // int totBdau = 0;
	    int totChrmDau[] = {0,0}; //array
	    int ChrmDauLund[2][4] = {{0,0,0,0},
				     {0,0,0,0} };
	    int BNrmDauLundId[] = {0,0,0,0,0,0}; //array
	    kAdaugIter.rewind();
	    BtaCandidate * tmpCand(0);
	    
	    while(tmpCand = kAdaugIter()) {
	      int tmpCandId = tmpCand->pdtEntry()->lundId();
	      if(tmpCandId ==PdtLund::gamma && 
		 tmpCand->energy()<0.015) continue;
	      totNrm += 1;
	      	      
	      if(abs(tmpCandId)%1000==411 || abs(tmpCandId)%1000==421) {
		totChrm += 1;
		if(tmpNrm <= 5) { // at most 6
		  HepAList<BtaCandidate> ChrmDescents;
		  findStableDescendents(tmpCand, ChrmDescents);
		  if(tmpChrm <=1) {
		    totChrmDau[tmpChrm] = ChrmDescents.length();
		    for (int i = 0; i < ChrmDescents.length(); i++){
		      BtaCandidate * tmpChrmDau(0);
		      tmpChrmDau = ChrmDescents[i];	
		      if(i<4) ChrmDauLund[tmpChrm][i]=
				tmpChrmDau->pdtEntry()->lundId();
		    }
		  }	
		  tmpChrm += 1; 	
		}              
	      }
	      if(tmpNrm <= 5) { // at most 6
		BNrmDauLundId[tmpNrm] = tmpCandId;
	      } 	
	      tmpNrm +=1;
	    } //end while line 3019
	    //fill in the ntuple of Kaon Ancestor information
	    _d0KchCand._TotBDau = totNrm;
	    _d0KchCand._TotChrmDau[0] = totChrmDau[0];
	    _d0KchCand._TotChrmDau[1] = totChrmDau[1]; 
	    for( int i = 0; i<2; i++) {
	      for( int j = 0; j<4; j++ ) {
		_d0KchCand._ChrmDau[i][j] = ChrmDauLund[i][j];
	      }
	    } 
	    for( int i=0; i<6; i++) {
	      _d0KchCand._BDau[i] = BNrmDauLundId[i];
	    }
	    totNrm = 0;
	  }// if for KaonAncest	
	}//if for OrigOfKch
      }
      _d0KchCand._HdTrkTruth = Kchflg;
      
      //tell ks or pi0 origin: 
      int OrigPi0(0);
      BtaCandidate * MothOfPi0(0);
      if( pi0UP != 0 ) {
	BtaCandidate * OrigOfPi0(0);
	OrigOfPi0 = _truthMap->mcFromReco(pi0UP);
	int OrigPi0Id(0);
	if(OrigOfPi0 != 0 ) OrigPi0Id = OrigOfPi0->pdtEntry()->lundId();
	if(OrigPi0Id == PdtLund::pi0) {
	  BtaCandidate * GMothOfPi0(0);
	  int MothOfPi0Id(0);
	  int GMothOfPi0Id(0);	
	  MothOfPi0 = OrigOfPi0->theMother();
	  if(0 != MothOfPi0) {
	    GMothOfPi0 = MothOfPi0->theMother();	
	    MothOfPi0Id = MothOfPi0->pdtEntry()->lundId();
	    if(GMothOfPi0 != 0 )
	      GMothOfPi0Id = GMothOfPi0->pdtEntry()->lundId();
	  }
	  
	  if(abs(MothOfPi0Id) == PdtLund::D0 || 
	     abs(GMothOfPi0Id) == PdtLund::D0) {
	    BtaCandidate * aD0Cand(0);
	    if(abs(MothOfPi0Id) == PdtLund::D0) {
	      aD0Cand = OrigOfPi0->theMother();
	      OrigPi0 = 10;
	    }
	    if(abs(GMothOfPi0Id) == PdtLund::D0) {
	      aD0Cand = MothOfPi0->theMother();
	      OrigPi0 = 20;
	    } 
	    bool Dto2Ppi0 = findD0Descendents(aD0Cand,PdtLund::pi_plus,
					      PdtLund::pi_minus, 
					      PdtLund::pi0, false);
	    if( Dto2Ppi0 ) OrigPi0 += 1;
	    bool Dto2Kpi0  = findD0Descendents(aD0Cand,PdtLund::K_plus,
					       PdtLund::K_minus, 
					       PdtLund::pi0, false);
	    if( Dto2Kpi0 ) OrigPi0 += 2;
	    bool DtoKpPpi0  = findD0Descendents(aD0Cand,PdtLund::K_plus,
						PdtLund::pi_minus, 
						PdtLund::pi0, false);
	    if( DtoKpPpi0 ) OrigPi0 += 3;
	    bool DtoKmPpi0  = findD0Descendents(aD0Cand,PdtLund::K_minus,
						PdtLund::pi_plus, 
						PdtLund::pi0, false);
	    if( DtoKmPpi0 ) OrigPi0 += 4;
	    
	    if( (!Dto2Ppi0) && (!Dto2Kpi0) && (!DtoKpPpi0) 
		&& (!DtoKmPpi0) ) OrigPi0 += 5;  
	  } else { // else of pi0 coming from a D0
	    OrigPi0 = 6;
	  }
	}//endif OrigOfPi0Id == PdtLund::pi0
	// 1-4 good, 5 not good but comes from d0, 
	// 6 is a pion but not comes from D0, 0 is bad. 
      }// end of if pi0 != 0
      _d0KchCand._OrigPi0 = OrigPi0; 
      //Ks		
      int OrigKs(0);
      if(ksUP != 0 ) {
	BtaCandidate * OrigOfKs(0);
	OrigOfKs = _truthMap->mcFromReco(ksUP);
	int OrigOfKsId(0);
	if(OrigOfKs != 0 ) OrigOfKsId = OrigOfKs->pdtEntry()->lundId();	       
	if(OrigOfKsId == PdtLund::K_S0) {  
	  BtaCandidate * MothOfKs(0);
	  BtaCandidate * GMothOfKs(0);
	  int MothOfKsId(0);
	  int GMothOfKsId(0);
	  MothOfKs = OrigOfKs->theMother();
	  if( 0 != MothOfKs) {			
	    MothOfKsId = MothOfKs->pdtEntry()->lundId();
	    GMothOfKs = MothOfKs->theMother();
	    if(GMothOfKs != 0) 
	      GMothOfKsId = GMothOfKs->pdtEntry()->lundId();
	  }
	  
	  if(abs(MothOfKsId) == PdtLund::D0 || 
	     abs(GMothOfKsId) == PdtLund::D0) {
	    BtaCandidate * aD0Cand(0);
	    if(abs(MothOfKsId) == PdtLund::D0)
	      aD0Cand = OrigOfKs->theMother();
	    if(abs(GMothOfKsId) == PdtLund::D0) 
	      aD0Cand = MothOfKs->theMother();
	    bool Dto2Pks = findD0Descendents(aD0Cand,PdtLund::pi_plus,
					     PdtLund::pi_minus, 
					     PdtLund::K_S0, false);
	    if( Dto2Pks ) OrigKs = 1;
	    bool Dto2Kks  = findD0Descendents(aD0Cand,PdtLund::K_plus,
					      PdtLund::K_minus, 
					      PdtLund::K_S0, false);
	    if( Dto2Kks ) OrigKs = 2;
	    bool DtoKpPks  = findD0Descendents(aD0Cand,PdtLund::K_plus,
					       PdtLund::pi_minus, 
					       PdtLund::K_S0, false);
	    if( DtoKpPks ) OrigKs = 3;
	    bool DtoKmPks  = findD0Descendents(aD0Cand,PdtLund::K_minus,
					       PdtLund::pi_plus, 
					       PdtLund::K_S0, false);
	    if( DtoKmPks ) OrigKs = 4;
	    if( (!Dto2Pks) && (!Dto2Kks) && (!DtoKpPks) 
		&& (!DtoKmPks) ) OrigKs = 5;  
	  } else { //else of Ks comes from D0:
	    OrigKs = 6;
	  }
	} //endif of Ks id 
	// 1-4, good, 5, comes from a d0, 6, is a pion,
	// but not from a d0, 0 bad pion.
      }//endif ks != 0
      _d0KchCand._OrigKs  = OrigKs;

      // tell the truth of the reconstructed B particle
      int truBflg = 0;	
      int DecInChrm = -1;  
      if(exclBch != 0 ) {
	const BtaCandidate * truBch(0);
	truBch = _truthMap->mcFromReco(exclBch);
	int BId(0);
	if(truBch !=0 ) BId = truBch->pdtEntry()->lundId();
	if( abs(BId) == PdtLund::B_plus || abs(BId) == PdtLund::B0 ) {
	  if( PdtLund::B_plus == BId ) truBflg = 1;
	  if( PdtLund::B_minus == BId ) truBflg = 2;
	  if( PdtLund::B0 == BId ) truBflg = 3;
	  if( PdtLund::anti_B0 == BId ) truBflg = 4;
	}
	if(_debug.value() && truBflg == 1 && OrigD0 > 0) {
	  BtaCandidate * tmpD = _truthMap->mcFromReco(exclD0);
	  Hep3Vector trueBDVtxVect = truBch->decayVtx()->point()
	    - tmpD->decayVtx()->point();
	  double TruCosBDVtxDmom = 
	    cos(trueBDVtxVect.angle(tmpD->p4().vect()));
	  _d0KchCand._TruCosBDVtxDmom = TruCosBDVtxDmom;
	}  			
	
	HepAList<BtaCandidate> RecBDecList; 
	findStableDescendents(exclBch, RecBDecList);
	if(RecBDecList.length() == 4 ) {
	  BtaCandidate * ListD1 = _truthMap->mcFromReco(RecBDecList[0]);
	  BtaCandidate * ListD2 = _truthMap->mcFromReco(RecBDecList[1]);
	  BtaCandidate * ListD3 = _truthMap->mcFromReco(RecBDecList[2]);
	  BtaCandidate * ListD4 = _truthMap->mcFromReco(RecBDecList[3]); 
	  
	  if(ListD1 != 0 && ListD2 != 0 && ListD3 != 0 && ListD4 != 0 ) {
	    int moth1Id(0), moth2Id(0),moth3Id(0),moth4Id(0);
	    int Gmoth1Id(0), Gmoth2Id(0),Gmoth3Id(0);	 	  
	    if(ListD1->theMother() != 0 ) {
	      BtaCandidate * moth1 = ListD1->theMother();
	      BtaCandidate * gmoth1 = moth1->theMother();
	      moth1Id = moth1->pdtEntry()->lundId();
	      if(gmoth1 != 0) { Gmoth1Id = gmoth1->pdtEntry()->lundId(); }
	    }	
	    if(ListD2->theMother() != 0 ) {
	      BtaCandidate * moth2 = ListD2->theMother();
	      BtaCandidate * gmoth2 = moth2->theMother();	
	      moth2Id = moth2->pdtEntry()->lundId();
	      if(gmoth2 != 0 ) { Gmoth2Id = gmoth2->pdtEntry()->lundId(); }
	    }	
	    if(ListD3->theMother() != 0 ) {
	      BtaCandidate * moth3 = ListD3->theMother();
	      BtaCandidate * gmoth3 = moth3->theMother();
	      moth3Id = moth3->pdtEntry()->lundId();
	      if(gmoth3 != 0 ) { Gmoth3Id = gmoth3->pdtEntry()->lundId(); }
	    }
	    if(ListD4->theMother() != 0 )
	      moth4Id = ListD4->theMother()->pdtEntry()->lundId();
	    BtaCandidate * Ancest1 = ancestor(ListD1);
	    BtaCandidate * Ancest2 = ancestor(ListD2);
	    BtaCandidate * Ancest3 = ancestor(ListD3);
	    BtaCandidate * Ancest4 = ancestor(ListD4);
	    if(Ancest1 != 0 && Ancest2 != 0 && 
	       Ancest3 != 0 && Ancest4 != 0 ) { 
	      if( Ancest1==Ancest2 && Ancest1== Ancest3
		  && Ancest3 == Ancest4 ){
		HepAList<BtaCandidate> TrubDecList; 
		findStableDescendents(Ancest1, TrubDecList);
		if(TrubDecList.length() == 4 ) {	
		  //Charm Decay:
		  if(abs(moth1Id) == PdtLund::D0 || abs(moth2Id) ==
		     PdtLund::D0 || abs(moth3Id) == PdtLund::D0 || 
		     abs(moth4Id) == PdtLund::D0 || abs(Gmoth1Id) ==
		     PdtLund::D0 || abs(Gmoth2Id) == PdtLund::D0 ||
		     abs(Gmoth3Id) == PdtLund::D0  ) {
		    DecInChrm = 1;
		  }
		  //Charmless decay
		  if(abs(moth1Id) != PdtLund::D0 && abs(moth2Id) != 
		     PdtLund::D0 && abs(moth3Id) != PdtLund::D0 &&
		     abs(moth4Id) != PdtLund::D0 && abs(Gmoth1Id) != 
		     PdtLund::D0 && abs(Gmoth2Id) != PdtLund::D0 &&
		     abs(Gmoth3Id) != PdtLund::D0 ) {
		    DecInChrm = 0;
		  }
		}	
	      }
	    }
	  }
	}	
      }  
      
      _d0KchCand._truBflg = truBflg;
      _d0KchCand._decInChrm = DecInChrm;
      // tell the truth of the reconstructed particles
      if(exclBch != 0 && exclD0 !=0 && exclKch != 0 ) {	
	const BtaCandidate * trueD0(0);
	const BtaCandidate * trueKch(0); 
        if (trueDToKsX>0 &&(OrigPi0>=10 && OrigPi0<20)) {
          if (ppGGMother->overlaps((const BtaCandidate&)(*MothOfPi0))){
	    if (findD0Descendents(ppGGMother, PdtLund::K_S0,PdtLund::pi0)){
	      trueD0 = ppGGMother;
	    }
	  }
	} else {          
	  trueD0 = _truthMap->mcFromReco(exclD0);
	}
	trueKch = _truthMap->mcFromReco(exclKch);
	if (isCategory((HepString)(_evtCat.value()), 
		       (int)BdecModMC[0], (int)BdecModMC[1],
		       OrigD0, D0RecDec)) {
	  //turn off during normal running
	  if(_evtCat.value() == "All" ||  
 	     (exclDeltaE < -0.07 || exclDeltaE >0.06) ||
 	     (exclD0Mass < 1.830 || exclD0Mass >1.895) ||
 	     HdTrkKaonNN <0.5 || 
 	     KpPidBit >= 19 || KmPidBit >=19 ) {
	    // do nothing 
	  } else {
	    //print decay tree
	    std::vector<const BtaCandidate*> matches;
	    matches.push_back(mcPi0Dau1);
	    matches.push_back(mcPi0Dau2);
	    matches.push_back(truePipTmp);
	    matches.push_back(truePimTmp);
	    matches.push_back(trueKch); 
	    
	    //HepAListIterator<BtaCandidate> mcIter(*_mcList);
	    bIter.rewind();
	    const BtaCandidate * upsilon= 0;
	    BtaCandidate * mcCand = 0;
	    while (mcCand = bIter()) {
	      PdtLund::LundType mcID = mcCand->pdtEntry()->lundId();
	      if (mcID == PdtLund::Upsilon_4S) {
		upsilon = mcCand;
	      }
	    }

	    if(exclMES>5.272){
	      _mES_Sig<<"DeltaE="<<exclDeltaE<<"   mES="<<exclMES
		      <<"    D0Mass="<<exclD0Mass<<endl;
	      _mES_Sig<<_printTree.print(*upsilon, matches)<<endl;
	    } else {
	      _mES_Sideband<<"DeltaE="<<exclDeltaE<<"   mES="<<exclMES
		      <<"    D0Mass="<<exclD0Mass<<endl;
	      _mES_Sideband<<_printTree.print(*upsilon, matches)<<endl;
	    }
	    matches.clear();
	  }
	}

	int mcdMomId = 0;
	int mcdGMomId = 0;
	int trueKchId = 0;
	int mcKchMomId = 0;
	int mcKchGMomId = 0;
	int nBDau = 0;
	int nDstarDau=0;
	int nKstarDau=0;
	
	if( (trueD0 != 0) && (trueKch != 0) ) {
	  // D0 Candidate 	     		    	      	
	  BtaCandidate * mcdMom(0);
	  const BtaCandidate * mcdAncest(0);
	  BtaCandidate * mcdGrandMom(0);
	  //unused variable
	  //BtaCandidate * mcdTemp(0);
	  //int isD0TempTrue = 0;
	  
	  mcdMom = const_cast<BtaCandidate *>(trueD0->theMother());
	  // find whether this D0 comes from B
	  mcdAncest = ancestor(const_cast<BtaCandidate *>(trueD0)); 
	  // cout<<" D Ancestor is: "<<mcdAncest->pdtEntry()->name()<<endl;
	  if(mcdAncest!=0) nBDau = mcdAncest->nDaughters();
	  if(mcdMom != 0  )  {
	    mcdMomId = mcdMom->pdtEntry()->lundId();
	    if(abs(mcdMomId)==PdtLund::D_star0||
	       abs(mcdMomId)==PdtLund::D_star_plus) 
	      nDstarDau=mcdMom->nDaughters();
	    mcdGrandMom = mcdMom->theMother();
	    if(mcdGrandMom != 0) {
	      mcdGMomId = mcdGrandMom->pdtEntry()->lundId();
	    }
	  }

	  // Kch Candidate		     
	  const BtaCandidate * mcKchMom(0);
	  const BtaCandidate * mcKchAncest(0);
	  BtaCandidate * mcKchGrandMom(0);
	  trueKchId = trueKch->pdtEntry()->lundId();
	  // find whether the ancestor is a B
	  mcKchAncest = ancestor(const_cast<BtaCandidate *>(trueKch));

	  mcKchMom = trueKch->theMother();
	  if(mcKchMom != 0 && mcKchAncest != 0 ) {
	    mcKchMomId = mcKchMom->pdtEntry()->lundId();
	    if(abs(mcKchMomId)==PdtLund::K_star_plus) 
	      nKstarDau = mcKchMom->nDaughters();
	    mcKchGrandMom = const_cast<BtaCandidate *>(mcKchMom->theMother());
	    if(mcKchGrandMom != 0) {
	      mcKchGMomId = mcKchGrandMom->pdtEntry()->lundId();
	    }		     
	  }

	  // combined all the information together:
	  int dFlg=0;
	  if( (mcdAncest != 0) && (mcKchAncest != 0) && 
	      (mcdAncest->overlaps(*mcKchAncest)) ) {
	    int  mcDTmpId = 0;
	    mcDTmpId = trueD0->pdtEntry()->lundId();
	    //cout<<" D0 id: mc "<<mcDTmpId<<" rec "
	    //<<exclD0->pdtEntry()->lundId()<<endl;
	    if(dFlg==0&&mcDTmpId == PdtLund::D0) dFlg = 1;
	    if(dFlg==0&&mcDTmpId == PdtLund::anti_D0) dFlg = -1;
	    int DecTruth = 0;
	    //reco Kch is identified as a kaon
	    //B+ -> D0 K+(signal), B+ -> D0 K*+
	    if(nBDau==2&&(abs(mcdMomId)==PdtLund::B_plus||
			  abs(mcdMomId)==PdtLund::B0)&&
	       abs(trueKchId)==PdtLund::K_plus) 
	      {
		// true signal for B+ -> D0 K+
		if(DecTruth==0&&abs(mcKchMomId)==PdtLund::B_plus
		   ||abs(mcKchMomId)==PdtLund::B0) DecTruth = 10;
		// cross feed by B+ ->D0 K*+ (pi0 K+), background
		if(DecTruth==0&&(abs(mcKchMomId)==PdtLund::K_star_plus)
		   &&nKstarDau==2&&
		   (abs(mcKchGMomId)==PdtLund::B_plus
		    ||abs(mcKchGMomId)==PdtLund::B0) ) DecTruth = 110;
	      }
	    //B0 -> D*-(D0 pi-) K+	
	    if(nDstarDau==2&&abs(mcdMomId)==PdtLund::D_star_plus
	       &&abs(mcdGMomId)==PdtLund::B0&&abs(trueKchId)==
	       PdtLund::K_plus){ 
	      if(nBDau==2&&DecTruth==0&&abs(mcKchMomId)==PdtLund::B0) 
		DecTruth = 120;
	    }
	    //B+ -> D*0 K+ or B+ -> D*0 K*+
	    if(nBDau==2&&abs(mcdMomId)==PdtLund::D_star0&&
	       (abs(mcdGMomId)==PdtLund::B_plus||
		abs(mcdGMomId)==PdtLund::B0)&&
	       abs(trueKchId)==PdtLund::K_plus) {	
	      // cross feed by B+ -> D0*(D0 pi0/gamma) K+, background
	      if(DecTruth==0&&abs(mcKchMomId)==PdtLund::B_plus) 
		DecTruth = 130;
	      // cross feed by B+ -> D0*(D0 pi0/gamma) K*+(pi0 K+), background
	      if(nDstarDau==2&&nKstarDau==2&&DecTruth==0&&
		 abs(mcKchMomId)==PdtLund::K_star_plus&&
		 abs(mcKchGMomId)==PdtLund::B_plus) DecTruth = 140;
	    }
		 
	    //reco Kch is misidentified, comes from a charge pion
	    //B+->D0 pi+ or B+ -> D0 rho+(pi+ pi0) this pi+ comes from Rho
	    if(nBDau==2&&(abs(mcdMomId)==PdtLund::B_plus||
			  abs(mcdMomId)==PdtLund::B0)&&
	       abs(trueKchId)==PdtLund::pi_plus)
	      {		
		// pi+ is from B+ -> D0 pi+
		if(DecTruth==0&&abs(mcKchMomId)==PdtLund::B_plus) 
		  DecTruth = 150;
		// pi+ comes from a charge Rho, B+ -> D0 rho+(pi0 pi+)
		if(DecTruth==0&&abs(mcKchMomId)==PdtLund::rho_plus&&
		   abs(mcKchGMomId)==PdtLund::B_plus) DecTruth = 160;
	      }
	    //B+ -> D*0 pi+
	    if(abs(mcdMomId)==PdtLund::D_star0&&
	       (abs(mcdGMomId)==PdtLund::B_plus||
		abs(mcdGMomId)==PdtLund::B0)&&
	       nDstarDau==2&&abs(trueKchId)==PdtLund::pi_plus) {
	      // pi+ from B+ ->D*0 pi+
	      if(nBDau==2&&DecTruth==0&&abs(mcKchMomId)==PdtLund::B_plus) 
		DecTruth = 170;
	      //cout<<"Dec code for B+->D* pi+: "<<DecTruth<<endl;
	    }
	    //B0 -> D*- pi+
	    if(nDstarDau==2&&abs(mcdMomId)==PdtLund::D_star_plus&&
	       (abs(mcdGMomId)==PdtLund::B_plus||abs(mcdGMomId)==PdtLund::B0)
	       &&abs(trueKchId)==PdtLund::pi_plus) {
	      // pi+ from B0 ->D*- pi+
	      if(nBDau==2&&DecTruth==0&&abs(mcKchMomId)==PdtLund::B0) 
		DecTruth = 180;
	    }				
	    //pi+ comes from K*+(Ks,pi+), B+ ->D0 K*+
	    if(nBDau==2&&(abs(mcdMomId)==PdtLund::B_plus||
			  abs(mcdMomId)==PdtLund::B0)&&
	       abs(trueKchId)==PdtLund::pi_plus){
	      if(nKstarDau==2&&DecTruth==0&&
		 abs(mcKchMomId)==PdtLund::K_star_plus&&
		 (abs(mcdGMomId)==PdtLund::B_plus
		  ||abs(mcKchGMomId)==PdtLund::B0)) DecTruth = 190;
	    }
	    //pi+ comes from K*+(Ks,pi+), B+ -> D*0 K*+	
	    if(nBDau==2&&(abs(mcdGMomId)==PdtLund::B_plus||
			  abs(mcdGMomId)==PdtLund::B0)&&
	       abs(trueKchId)==PdtLund::pi_plus) {
	      if(nKstarDau==2&&DecTruth==0&&
		 abs(mcKchMomId)==PdtLund::K_star_plus&&
		 (abs(mcdGMomId)==PdtLund::B_plus
		  ||abs(mcKchGMomId)==PdtLund::B0)) DecTruth = 200;	
	    }	

	    if(DecTruth >= 10 ) {
	      const BtaCandidate * trueKp = 0;
	      const BtaCandidate * trueKm = 0;
	      const BtaCandidate * trueKs = 0;
	      const BtaCandidate * truePip = 0;
	      const BtaCandidate * truePim = 0;
	      const BtaCandidate * truePi0 = 0;

	      int tmpValue = 0;
	      int tmpDec = 0;
	      if(KplusUP != 0 ) {
		trueKp = _truthMap->mcFromReco(KplusUP);
	      }
	      if(KminusUP != 0) {
		trueKm = _truthMap->mcFromReco(KminusUP);
	      }
	      if(ksUP != 0) {
		trueKs = _truthMap->mcFromReco(ksUP);
	      }
	      if(PiplusUP != 0) {
		truePip = _truthMap->mcFromReco(PiplusUP);
	      }
	      if(PiminusUP!= 0) {
		truePim = _truthMap->mcFromReco(PiminusUP);
	      }
	      if(pi0UP != 0) {
		truePi0 = _truthMap->mcFromReco(pi0UP);
	      }
	      
	      if( trueKp!=0 && trueKm!=0 && truePi0!=0 ) {
		int KpId(0); KpId = trueKp->pdtEntry()->lundId();
		int KmId(0); KmId = trueKm->pdtEntry()->lundId();
		int Pi0Id(0); Pi0Id = truePi0->pdtEntry()->lundId();
		if(KpId==PdtLund::K_plus&&KmId==PdtLund::K_minus&&
		   Pi0Id==PdtLund::pi0)  tmpValue = 1;
		if(KpId==PdtLund::pi_plus&&KmId==PdtLund::K_minus&&
		   Pi0Id==PdtLund::pi0) tmpDec = 1;
		if(KpId==PdtLund::K_plus&&KmId==PdtLund::pi_minus&&
		   Pi0Id==PdtLund::pi0) tmpDec = -2;
	      }
	      
	      if( trueKp!=0 && truePim!=0 && trueKs!=0 ) {
		int KpId(0); KpId = trueKp->pdtEntry()->lundId(); 
		int PimId(0); PimId = truePim->pdtEntry()->lundId();
		int KsId(0); KsId = trueKs->pdtEntry()->lundId();
		if(KpId==PdtLund::K_plus&&PimId==PdtLund::pi_minus&&
		   KsId==PdtLund::K_S0) {
		  if(dFlg==1) tmpValue=2;
		  if(dFlg==-1) tmpValue=3;
		}
		if(KpId==PdtLund::pi_plus&&PimId==PdtLund::pi_minus&&
		   KsId==PdtLund::K_S0) tmpDec = 3;
		if(KpId==PdtLund::K_plus&&PimId==PdtLund::K_minus&&
		   KsId==PdtLund::K_S0)   tmpDec = -4;
	      }
	      if( trueKm!=0 && truePip!=0 && trueKs!=0 ) {
		int KmId(0);  KmId = trueKm->pdtEntry()->lundId(); 
		int PipId(0); PipId = truePip->pdtEntry()->lundId();
		int KsId(0);  KsId = trueKs->pdtEntry()->lundId();
		if(KmId==PdtLund::K_minus&&PipId==PdtLund::pi_plus&&
		   KsId==PdtLund::K_S0) {		
		  if(dFlg==1) tmpValue=4;
		  if(dFlg==-1) tmpValue=5;
		}
		if(KmId==PdtLund::K_minus&&PipId==PdtLund::K_plus&&
		   KsId==PdtLund::K_S0)   tmpDec = 4;
		if(KmId==PdtLund::pi_minus&&PipId==PdtLund::pi_plus&&
		   KsId==PdtLund::K_S0) tmpDec = -3;
	      }
	      if( trueDToKsX<0 && truePip!=0 && truePim!=0 && truePi0!=0 ) {
		int PipId(0); PipId = truePip->pdtEntry()->lundId();
		int PimId(0); PimId = truePim->pdtEntry()->lundId();
		int Pi0Id(0); Pi0Id = truePi0->pdtEntry()->lundId();
		if(PipId==PdtLund::pi_plus&&PimId==PdtLund::pi_minus&&
		   Pi0Id==PdtLund::pi0) tmpValue = 6;
		if(PipId==PdtLund::K_plus&&PimId==PdtLund::pi_minus&&
		   Pi0Id==PdtLund::pi0) tmpDec = 2; 
		if(PipId==PdtLund::pi_plus&&PimId==PdtLund::K_minus&&
		   Pi0Id==PdtLund::pi0) tmpDec = -1;
	      }
	      
	      if( trueDToKsX>0 && truePip!=0 && truePim!=0 && truePi0!=0 ) {
		int PipId(0); PipId = truePip->pdtEntry()->lundId();
		int PimId(0); PimId = truePim->pdtEntry()->lundId();
		int Pi0Id(0); Pi0Id = truePi0->pdtEntry()->lundId();
		if(PipId==PdtLund::pi_plus&&PimId==PdtLund::pi_minus&&
		   Pi0Id==PdtLund::pi0) tmpValue = 7;
		if(PipId==PdtLund::K_plus&&PimId==PdtLund::pi_minus&&
		   Pi0Id==PdtLund::pi0) tmpDec = 2; 
		if(PipId==PdtLund::pi_plus&&PimId==PdtLund::K_minus&&
		   Pi0Id==PdtLund::pi0) tmpDec = -1;
	      }

	      
	      _d0KchCand._recD0CaDec = tmpDec;
	      // 1. D0 -> K- pi+ pi0 +/-1 + positive trk replaced
	      // 2. D0 -> K+ pi- pi0 +/-2 + positive trk replaced
	      // 3. D0 -> pi+ pi- ks +/-3 + positive trk replaced
	      // 3. D0 -> K+ K- Ks   +/-4 + positive trk replaced
	      
	      _d0KchCand._exclTruth = DecTruth + tmpValue;
	      // B+ -> D0 K+   10, true signal with tmpValue >1;
	      // B+ -> D0 Rho+(pi+ pi0) 160, background type 1, drop one pi0 
	      // B+ -> D0 pi+  150, background type 2
	      // B+ -> D0 K*+(K+ pi0)  K+, 110 background type 3, pi0 dropped
	      // B+ -> D*0(D0 pi0/gamma) K+, 130, background type 4, pi0/gamma is droped
	      // B+ -> D*0(D0 pi0/gamma) K*+(K+ pi0) 140, background type 5, 2 pi0/gamma is droped
	      // B0 -> D*-(D0 pi-) K+ 120;
	    } //end of DecTruth >= 10
	  } //end of D0 Kch comes from the same ancestor
	} // end of trueB, trueD0, trueKch
      } // end of exclBch, exclD0 and exclKch 
    }// end of isMC()

    //Append the _d0KchCand to the Vector
    d0KchCandidates.push_back(_d0KchCand);			       
  } // end of the loop of exclList ( d0KchCandidates )
    //deal with the vector stuff,  
  const int numCands = d0KchCandidates.size();
  static const float zF = 0;
  // Ntuple event level data:
  if(numCands > 0 ) {
    _tuple->column("eventNumber", (int)_numEventsRead, 0, EVENT_BLOCK_NAME);
    _tuple->column("runNumber",(int)eventID->run(), 0, EVENT_BLOCK_NAME);
    _tuple->column("lower",(int) theTriplet.timeStamp().binary().lower, 
		   0,EVENT_BLOCK_NAME);
    _tuple->column("upper",(int) theTriplet.timeStamp().binary().upper, 
		   0,EVENT_BLOCK_NAME);
    _tuple->column("R2", (float)foxWolfR2, zF, EVENT_BLOCK_NAME);
    if( isMC() ) {
      _tuple->column("B1decMode", (int)BdecModMC[0], 0, EVENT_BLOCK_NAME);
      _tuple->column("B2decMode", (int)BdecModMC[1], 0, EVENT_BLOCK_NAME);
      _tuple->column("D1CaDecMode", (int)D0CADecMode[0], 0, EVENT_BLOCK_NAME);
      _tuple->column("D2CaDecMode", (int)D0CADecMode[1], 0, EVENT_BLOCK_NAME);
      _tuple->column("bVtxDistZ",(float)mcbDistZ, zF, EVENT_BLOCK_NAME); 
    }
    _tuple->column("thrust",(float)thrust,zF,EVENT_BLOCK_NAME);
    _tuple->column("hem1Mass",(float)hemMass1,zF,EVENT_BLOCK_NAME);
    _tuple->column("hem2Mass",(float)hemMass2,zF,EVENT_BLOCK_NAME);
    
    if (false == usingOffResEBeam){
      _tuple->column("eBeam", (float)_eBeam, zF, EVENT_BLOCK_NAME);
    }
    else {
      _tuple->column("eBeam", -(float)_eBeam, zF, EVENT_BLOCK_NAME);
    }
    
    ntuple(d0KchCandidates);

    // Dump the ntuples:
    if (true == _dumpNtuple.value()){
      _tuple->dumpData();
    }
  }	   
  d0KchCandidates.clear();
  return AppResult::OK;
}    

// protected functions:
//-------------------------------------------------------------------

bool B2D0KchNonCPAnal::isMC() const {
  static const bool result =
    0 != _truthMap && 0 != _mcList && _mcList->length() > 0;
  
  if (true == _treatAsMC.value()) {
    return true;
  }
  else {
    return result;
  }
}

//------------------------------------------------------------------
bool B2D0KchNonCPAnal::isCategory(HepString cat, int bmod1, 
				  int bmod2, int dflg, int drec){
  bool result = false;
  
  if (cat == "All") {
    result = true;
  } else if (cat == "DKBadD") {
    result = (bmod1==16 || bmod2==16) && dflg<=0 && drec!=3;
  } else if (cat == "DKGoodD") {
    result = (bmod1==16 || bmod2==16) && dflg>0 && drec==3;
  } else if (cat == "DPiBadD") {
    result = (bmod1==156 || bmod2==156) && dflg<=0 && drec!=3;
  } else if (cat == "DPiGoodD") {
    result = (bmod1==156 || bmod2==156) && dflg>0 && drec==3;
  } else if (cat == "DPiX") {
    result = bmod1==50 || (bmod1>79 && bmod1<109) || 
             (bmod1>149 && bmod1<156) ||
             (bmod1>156 && bmod1<199) || 
             bmod2==50 || (bmod2>79 && bmod2<109) || 
             (bmod2>149 && bmod2<156) ||
             (bmod2>156 && bmod2<199);
  } else if (cat == "DKX") {
    result = (bmod1>9 && bmod1<16) || (bmod1>16 && bmod1<19) || 
             (bmod1>59 && bmod1<79) || 
             (bmod1>109 && bmod1<149) ||
             (bmod2>9 && bmod2<16) || (bmod2>16 && bmod2<19) || 
             (bmod2>59 && bmod2<79) || 
             (bmod2>109 && bmod2<149);
  } else if (cat == "BBBadD") {
    result = (bmod1>199 || bmod2>199) || 
             (bmod1<=0 && bmod2<=0 && 
              dflg<=0 && drec!=3);
  } else if (cat == "BBGoodD") {
    result = bmod1<=0 && bmod2<=0 && dflg>0 && drec==3;
  } else {
    result = false;
  }

  return result;
}

//------------------------------------------------------------------
BtaCandidate * B2D0KchNonCPAnal::ancestor(BtaCandidate * descen) {
  // No ancestor found if no descen, return 0
  if( 0 == descen ) return 0;
  // No ancestor found if no mother. This is the normal way to return 0:
  BtaCandidate * theMother = 0; 
  theMother = descen->theMother();
  if (0 == theMother) {
    return 0;
  }
  const PdtEntry * entry = theMother->pdtEntry();
  if (0 == entry){
    return 0;
  }
  //Ancestor found if mother = B0 or B+	
  if( abs(entry->lundId())==PdtLund::B0 || 
      abs(entry->lundId())==PdtLund::B_plus) {
    return theMother;
  }
  // Otherwise, call recursively on the mother:
  return ancestor(theMother);
}
 
//------------------------------------------------------------------
void B2D0KchNonCPAnal::CandVtxFit( const HepAList<BtaCandidate> & DauList,
				   bool useBeamSpotConstrain, 
				   bool useMassConstraint, float candMass) {
  //if the daulist is empty return 0;
  if(DauList.length() == 0 ) return;
  BtaOpMakeTree comb;
  BtaCandidate * candToFit = 0;
  
  if(DauList.length()==1) {
    candToFit = new BtaCandidate (*(DauList.first()));
  } else {
    HepAListIterator<BtaCandidate> IterDaughter(DauList);
    IterDaughter.rewind();		
    candToFit = comb.createFromList(IterDaughter);
  }
  
  //make sure candToFit is not zero
  if( 0 == candToFit ) return;
  
  setGeoConstraint(*candToFit);
  
  if (true == useBeamSpotConstrain) {
    setBeamConstraint(*candToFit, _beamSpot);
  }
  
  if ( true == useMassConstraint) {
    setMassConstraint(*candToFit, candMass);
  }
  
  //invalidate the possible existed fit at first
  candToFit->invalidatePresentFit();
  
  //candToFit should have more than one daughters
  VtxFitterOper::algType algorithm = VtxFitterOper::GeoKin;
  //if do not has daughters
  if( 0 == candToFit->nDaughters() ) {
    algorithm = VtxFitterOper::SingleTrackGeoKin;
  } 
  
  //vertex fit
  VtxFitterOper fitter(*candToFit, algorithm);
  fitter.fit();
  if (true == useMassConstraint) {
    _fittedCand = fitter.getFittedTree();
  } else{
    _fittedCand = fitter.getFitted(*candToFit);
  }
  
  delete candToFit;
}
	  
//------------------------------------------------------------------
bool  B2D0KchNonCPAnal::printCand(const char * particle) const {
  if (_printCands.value().length() > 0){
    std::istringstream candsStream(_printCands.value().c_str());
    while (1) {
      HepString cand; 
      candsStream >> cand;
      if (0 == cand.length()){
	return false;
      } 
      if (cand == particle){
	return true;
      }
    }
  }	
  return false;
}   

//-------------------------------------------------------------------
void B2D0KchNonCPAnal::dumptrack(BtaCandidate* track) {
  _pidHypo = -99.0;
  _nSvtHits = -99.0;
  _nDchHits = -99.0;
  _pocad0 = -99.0;
  _pocaz0 = -99.0;
  _KaonNN = -99.0;
  _NPhot = -99.0;
  _NBkPhot = -99.0;
  _sigoffprot = -99.0;
  _drcKaonCon = -99.0;
  _dchKaonCon = -99.0;
  _svtKaonCon = -99.0;
  _thC = -99.9;
  _thCErr = -99.9;
  
  if(track == 0 ) return;
  PidKaonMicroSelector* KaonSel = new PidKaonMicroSelector;
  const BtaTrkQual* TrkQual = track->getMicroAdapter()->getTrkQual();
  const BtaCalQual* CalQual = track->getMicroAdapter()->getCalQual();
  const BtaPidQual* PidQual = track->getMicroAdapter()->getPidQual();
  const BtaPidInfo* PidInfo = track->getMicroAdapter()->getPidInfo();
  //const BtaIfrQual* IfrQual = track->getMicroAdapter()->getIfrQual();
  
  //trk qual info
  if (TrkQual){
    _pidHypo  = TrkQual->pidHypo();
    _nSvtHits = TrkQual->nSvtHits();
    _nDchHits = TrkQual->nDchHits();
  }//end TrkQual
  
  // find the closest approach to the interaction point???
  
  if(track->charge() != 0 && track->trkAbsFit()){
    _pocad0=track->recoTrk()->fitResult()->helix(0.).d0();
    _pocaz0=track->recoTrk()->fitResult()->helix(0.).z0();
  }else if(CalQual){
    _pocad0 = sqrt((CalQual->centroid().x()*CalQual->centroid().x())
		   +(CalQual->centroid().y()*CalQual->centroid().y()));
    _pocaz0 = CalQual->centroid().z();
  }else {
    _pocad0 = 0;
    _pocaz0 = 0;
  }
  
  HTValVector<float> cele(6),cmu(6),cpi(6),cka(6),cpro(6);
  HTValVector<float> liele(6),limu(6),lipi(6),lika(6),lipro(6);
  //pid info
  if (PidQual){
    //_dEdXSvt=PidQual->dEdXSvt();
    //_nSamDeDxSvt=PidQual->nSamplesDeDxSvt();
    //_dEdXDch=PidQual->dEdXDch();
    //_nSamDeDxDch=PidQual->nSamplesDeDxDch();
    
    // other pid info added here::::
    
    _thC=PidQual->thetaC();
    _thCErr=PidQual->thetaCErr();
    _NPhot = PidQual->ringNPhot();
    _NBkPhot = PidQual->ringNBkgd();
    //test consistency:
    for (int j=0;j<=PidInfo->nDet();j++)
      {
	cele[j] = PidInfo->consistency(Pdt::lookup(PdtPdg::e_minus),
				       (PidSystem::System)j).consistency();
	cmu[j] = PidInfo->consistency(Pdt::lookup(PdtPdg::mu_minus),
				      (PidSystem::System)j).consistency();
	cpi[j] = PidInfo->consistency(Pdt::lookup(PdtPdg::pi_minus),
				      (PidSystem::System)j).consistency();
	cka[j] = PidInfo->consistency(Pdt::lookup(PdtPdg::K_minus),
				      (PidSystem::System)j).consistency();
	cpro[j] = PidInfo->consistency(Pdt::lookup(PdtPdg::p_plus),
				       (PidSystem::System)j).consistency();
	liele[j] = PidInfo->consistency(Pdt::lookup(PdtPdg::e_minus),
					(PidSystem::System)j).likelihood();
	limu[j] = PidInfo->consistency(Pdt::lookup(PdtPdg::mu_minus),
				       (PidSystem::System)j).likelihood();
	lipi[j] = PidInfo->consistency(Pdt::lookup(PdtPdg::pi_minus),
				       (PidSystem::System)j).likelihood();
	lika[j] = PidInfo->consistency(Pdt::lookup(PdtPdg::K_minus),
				       (PidSystem::System)j).likelihood();
	lipro[j] = PidInfo->consistency(Pdt::lookup(PdtPdg::p_plus),
					(PidSystem::System)j).likelihood();
      }

    //Svt consistency for e muon pion kaon proton
    //here just use kaon
    _svtKaonCon=lika[5];
    
    //Dch consistency
    //here just use kaon
    _dchKaonCon=lika[4];	  

    //Drc consistency
    //here just use kaon
    _drcKaonCon=lika[3];
    
    // calculate the sigma off from the expected particle
    HTValVector<int> nexPhotons(5);
    HTValVector<int> nexPhotFix(5);
    for (int k1=0;k1<5;k1++) {
      nexPhotons[k1] = PidQual->ringNExPhot(k1);
      nexPhotFix[k1] = nexPhotons[k1];
    }
    for (int k=PdtPid::muon; k<=PdtPid::proton; k++ ) {
      // unused variable
      //PdtPid::PidType part = (PdtPid::PidType)k;
      if ( nexPhotFix[k] > nexPhotFix[k-1]) {
	if (k<PdtPid::proton) {
	  nexPhotFix[k]=std::min((nexPhotFix[k-1]+
				  nexPhotFix[k+1])/2, nexPhotFix[k-1]);
	}else {
	  nexPhotFix[k]=nexPhotFix[k-1];
	}
      }
    }
    
    // _sigoffelec = calNPhot(cele[3],nexPhotFix(0));
    // _sigoffmuon = calNPhot(cmu[3],nexPhotFix(1));
    // _sigoffpion = calNPhot(cpi[3],nexPhotFix(2));
    // _sigoffkaon = calNPhot(cka[3],nexPhotFix(3));
    _sigoffprot = calNPhot(cpro[3],nexPhotFix(4));
    
    //kaon particle hypo signal Neural net output
    _KaonNN = KaonSel->nn_result(track);
  }//end PidQual   
  delete KaonSel;
  return; 
}//end dumptrack info  
     
//----------------------------------------------------------------------
void B2D0KchNonCPAnal::dumpGamma(BtaCandidate* gamma) {
  _GamEn = -1.0;
  _latMom = -1.0;
  _bumRawE = -1.0;
  if(gamma != 0) {
    _GamEn = gamma->energy();
    const BtaCalQual* CalQual = gamma->getMicroAdapter()->getCalQual();
    if (CalQual){
      _latMom = CalQual->lateralMoment();	
      _bumRawE = CalQual->rawEnergy();
    }
  }
}	
	
//-------------------------------------------------------------------
int B2D0KchNonCPAnal::calNPhot(double sigLv1,double nExp) {
  nExp=nExp+2;
  float sigLv1Tmp = 0;
  int i;
  if ( sigLv1 <= 0) {
    i=-1;
    while ( sigLv1Tmp < fabs(sigLv1/2) ) {
      i++;
      sigLv1Tmp +=exp(-nExp+i*log(nExp)-NumRecipes::gammln(i+1));
    }
  }  else {
    i=101;
    while ( sigLv1Tmp < fabs(sigLv1/2) ) {
      i--;
      sigLv1Tmp +=exp(-nExp+i*log(nExp)-NumRecipes::gammln(i+1));
    }
  }
  
  if ((i-2)>0) i-=2;
  else i=0;
  
  return i;
  
}

//-------------------------------------------------------------------
void
B2D0KchNonCPAnal::ntuple(const vector<D0KchCand> & cands) {
  
  const int nCands = cands.size();
  // Set the range for the vectors:
  _tuple->column(CANDS_INDEX_NAME, nCands, 0, CANDS_BLOCK_NAME,
                 HTRange<int>(0,150));
  // Ntuple tne vectors in _tuple
  
  //Fisher discriminant
  ntupleVar(cands,"BFisher", &D0KchCand::BFisher,_tuple);
  ntupleVar(cands,"BFCLEO", &D0KchCand::BFCLeo,_tuple);		
  ntupleVar(cands,"BFLegendre", &D0KchCand::BFLegendre,_tuple);		
  ntupleVar(cands,"BcanThR",&D0KchCand::BCanThrust,_tuple);
  ntupleVar(cands,"CosBThetaT", &D0KchCand::BCanCosThR,_tuple);
  ntupleVar(cands,"cosBmomAx", &D0KchCand::cosBmomAx,_tuple);
  ntupleVar(cands,"cosBThRAx", &D0KchCand::cosBThRAx,_tuple);

  //Hard trk quantities:
  ntupleVarInt(cands,"HdTrkChge",&D0KchCand::HdTrkChg,_tuple);
  ntupleVar(cands, "HdTrkTheta", &D0KchCand::HdTrkTheta, _tuple);
  ntupleVar(cands, "HdTrkCosTheta", &D0KchCand::HdTrkCosTheta, _tuple);
  ntupleVar(cands, "HdTrkPhi", &D0KchCand::HdTrkPhi, _tuple);
  ntupleVar(cands, "HdTrkMom", &D0KchCand::HdTrkMom, _tuple);
  ntupleVar(cands, "HdTrkThC", &D0KchCand::HdTrkThetaC, _tuple);
  ntupleVar(cands, "HdTrkThcErr", &D0KchCand::HdTrkThcErr, _tuple);	
  ntupleVar(cands, "HdTrkNPhot",&D0KchCand::HdTrkNPhot,_tuple);
  ntupleVar(cands, "HdTrkKaonNN", &D0KchCand::HdTrkKaonNN,_tuple);
  ntupleVarInt(cands, "HdTrkPid",&D0KchCand::HdTrkPid, _tuple);
  ntupleVarInt(cands, "HdTrkElecPid",&D0KchCand::HdTrkElecPid,_tuple);
  ntupleVarInt(cands, "HdTrkMuonPid",&D0KchCand::HdTrkMuonPid,_tuple);
  ntupleVar(cands, "HdTrkCMTheta",&D0KchCand::HdTrkCMTheta,_tuple);
  ntupleVar(cands, "HdTrkCMPhi", &D0KchCand::HdTrkCMPhi,_tuple);
  ntupleVar(cands, "HdTrkPcm", &D0KchCand::HdTrkCMMomtum, _tuple);
  // hard track truth:
  if( isMC() ) { 
    ntupleVarInt(cands, "HdTrkTruth", &D0KchCand::HdTrkTruth,_tuple);
    ntupleVarInt(cands, "TotBDau", &D0KchCand::TotBDau,_tuple);
    ntupleVarInt(cands, "TotChrmDau1", &D0KchCand::TotChrmDau1,_tuple);
    ntupleVarInt(cands, "TotChrmDau2", &D0KchCand::TotChrmDau2,_tuple);
    ntupleVarInt(cands, "ChrmDau11", &D0KchCand::ChrmDau11,_tuple);
    ntupleVarInt(cands, "ChrmDau12", &D0KchCand::ChrmDau12,_tuple);
    ntupleVarInt(cands, "ChrmDau13", &D0KchCand::ChrmDau13,_tuple);
    ntupleVarInt(cands, "ChrmDau14", &D0KchCand::ChrmDau14,_tuple);
    ntupleVarInt(cands, "ChrmDau21", &D0KchCand::ChrmDau21,_tuple);
    ntupleVarInt(cands, "ChrmDau22", &D0KchCand::ChrmDau22,_tuple);
    ntupleVarInt(cands, "ChrmDau23", &D0KchCand::ChrmDau23,_tuple);
    ntupleVarInt(cands, "ChrmDau24", &D0KchCand::ChrmDau24,_tuple);
    ntupleVarInt(cands, "BDau1", &D0KchCand::BDau1,_tuple); 	
    ntupleVarInt(cands, "BDau2", &D0KchCand::BDau2,_tuple); 	
    ntupleVarInt(cands, "BDau3", &D0KchCand::BDau3,_tuple); 	
    ntupleVarInt(cands, "BDau4", &D0KchCand::BDau4,_tuple); 	
    ntupleVarInt(cands, "BDau5", &D0KchCand::BDau5,_tuple); 	
    ntupleVarInt(cands, "BDau6", &D0KchCand::BDau6,_tuple); 	
  }
  //Hdtrk with another kaon trk, maximum and minium invariant mass and vtx
  ntupleVar(cands,"KbKUpMass", &D0KchCand::KbKUpMass, _tuple);
  ntupleVar(cands,"KbKLowMass", &D0KchCand::KbKLowMass, _tuple);
  ntupleVarInt(cands,"nRoEKp", &D0KchCand::nRoEKp, _tuple);
  ntupleVarInt(cands,"nRoEKm", &D0KchCand::nRoEKm, _tuple);		
  
  //Hd Trk with another neutral kaon(Rest of Event)
  ntupleVar(cands,"KbKsUpMass", &D0KchCand::KbKsUpMass, _tuple);
  ntupleVar(cands,"KbKsLowMass", &D0KchCand::KbKsLowMass, _tuple);	
  ntupleVarInt(cands,"nRoEKs", &D0KchCand::nRoEKs,_tuple);

  //HdTrk with lepton trk
  ntupleVar(cands,"KbLepMass", &D0KchCand::KbLepMass, _tuple);		

  //Hd Trk with a trk from ROE
  ntupleVar(cands,"KbRoETrkChisq", &D0KchCand::KbRoETrkChisq, _tuple);
  ntupleVar(cands,"VtxKbRoETrkDistY",&D0KchCand::VtxKbRoETrkDistY,_tuple);
  ntupleVar(cands,"VtxKbRoETrkDistChisqY",&D0KchCand::VtxKbRoETrkDistChisqY,_tuple);			
  //Hemisphere charge diff 
  ntupleVar(cands, "QHemiDiff", &D0KchCand::QDiff,_tuple);

  //a kaon from roe to B(signal) vtx doca
  ntupleVar(cands,"KroeBvtxDoca", &D0KchCand::KroeBvtxDoca,_tuple);	

  //a trk(no Ks daughter) from roe to B(signal) vtx doca
  ntupleVar(cands,"TrkRoeBvtxDoca", &D0KchCand::TrkRoeBvtxDoca,_tuple);

  //K_S0, one daughter of D0 in D0 -> K pi K_S0
  if( isMC() ) {
     ntupleVarInt(cands, "OrigKs", &D0KchCand::OrigKs,_tuple);
  }
  ntupleVar(cands, "KsMass", &D0KchCand::KsMass, _tuple);
  ntupleVar(cands, "KsHel", &D0KchCand::KsHel, _tuple);
  ntupleVar(cands, "KsVtxFitChisq", &D0KchCand::KsVtxFitChisq, _tuple);
  ntupleVar(cands, "KsDecLen", &D0KchCand::KsDecLen, _tuple);
  ntupleVar(cands, "KsTheta", &D0KchCand::KsTheta, _tuple);
  ntupleVar(cands, "pKs", &D0KchCand::pKs, _tuple);
  ntupleVar(cands, "cosDKsVtx", &D0KchCand::cosDKsVtx,_tuple);
  ntupleVar(cands, "dksVtxVctChisq", &D0KchCand::dksVtxVctChisq, _tuple);  

  //pi0, one daughter of D0 in D0 -> pi pi pi0, K K pi0
  if( isMC() ) {
     ntupleVarInt(cands,"OrigPi0",&D0KchCand::OrigPi0,_tuple);
  }
  ntupleVar(cands, "pi0Mass", &D0KchCand::pi0Mass, _tuple);
  ntupleVar(cands, "Gam1Energy", &D0KchCand::gamEng1, _tuple);
  ntupleVar(cands, "Gam2Energy", &D0KchCand::gamEng2, _tuple);
  //Mc true energy for photon
 if( isMC() ) {
   ntupleVar(cands, "Gam1McEgy", &D0KchCand::truPi0E1, _tuple);
   ntupleVar(cands, "Gam2McEgy", &D0KchCand::truPi0E2, _tuple);
 }
  ntupleVar(cands, "GamLat1", &D0KchCand::gamLat1, _tuple);
  ntupleVar(cands, "GamLat2", &D0KchCand::gamLat2, _tuple);
  ntupleVar(cands, "pi0Hel",&D0KchCand::pi0Hel,_tuple);
  ntupleVar(cands, "pi0HelDfm",&D0KchCand::pi0HelDfm,_tuple);
  ntupleVar(cands, "pi0CMmom", &D0KchCand::pi0CMmom,_tuple);
  ntupleVar(cands, "pi0Pcm", &D0KchCand::pi0CMmom0,_tuple);
  ntupleVar(cands, "pi0EAsy", &D0KchCand::pi0EAsy,_tuple);
  // one the photon from pi0 combined with other photon from the event
  ntupleVar(cands, "Best2GmaMass1", &D0KchCand::pi0BestMass1, _tuple);
  ntupleVar(cands, "Best2GmaMass2", &D0KchCand::pi0BestMass2, _tuple);
  ntupleVar(cands, "BestCos2GmaHel1", &D0KchCand::Bestcospi0Hel1, _tuple);
  ntupleVar(cands, "BestCos2GmaHel2", &D0KchCand::Bestcospi0Hel2, _tuple);
  ntupleVar(cands, "Best2gPcm1", &D0KchCand::Best2gPcm1,_tuple);
  ntupleVar(cands, "Best2gPcm2", &D0KchCand::Best2gPcm2,_tuple);	
  //track kaon pid info
  ntupleVarInt(cands, "KpPidBit", &D0KchCand::KpPidBit,_tuple);
  ntupleVarInt(cands, "KmPidBit", &D0KchCand::KmPidBit,_tuple);
  
  //track pion pid info
  ntupleVarInt(cands, "PipPidBit", &D0KchCand::PipPidBit,_tuple);
  ntupleVarInt(cands, "PimPidBit", &D0KchCand::PimPidBit,_tuple);

  //trk info
  ntupleVar(cands, "D0DauTrkP1", &D0KchCand::D0DauTrkP1,_tuple);
  ntupleVar(cands, "D0DauTrkP2", &D0KchCand::D0DauTrkP2,_tuple);

  ntupleVar(cands, "D0DauTrkTheta1", &D0KchCand::D0DauTrkTheta1,_tuple);
  ntupleVar(cands, "D0DauTrkTheta2", &D0KchCand::D0DauTrkTheta2,_tuple);

  ntupleVar(cands, "D0DauTrkPhi1", &D0KchCand::D0DauTrkPhi1,_tuple);
  ntupleVar(cands, "D0DauTrkPhi2", &D0KchCand::D0DauTrkPhi2,_tuple);

  ntupleVar(cands, "pipiDecLen", &D0KchCand::pipiDecLen, _tuple);
  ntupleVar(cands, "pipiMass", &D0KchCand::pipiMass, _tuple);
  ntupleVarInt(cands, "pipiKsTrue", &D0KchCand::pipiKsTrue, _tuple);

  if( isMC() ) {
     ntupleVarInt(cands, "piPlusTru", &D0KchCand::pipTru, _tuple); 
     ntupleVarInt(cands, "piMinusTru", &D0KchCand::pimTru, _tuple); 
  }

  //Dalitz varaibles from truth
  if( isMC() ) {
    ntupleVar(cands, "d0kkMCMass", &D0KchCand::d0kkMCMass, _tuple);
    ntupleVar(cands, "d0kkMCHelic", &D0KchCand::d0kkMCHelic, _tuple);
    ntupleVar(cands, "d0piKpMCMass", &D0KchCand::d0pKpMCMass, _tuple);
    ntupleVar(cands, "d0piKpMCHelic", &D0KchCand::d0pKpMCHelic, _tuple);
    ntupleVar(cands, "d0piKmMCMass", &D0KchCand::d0pKmMCMass, _tuple);
    ntupleVar(cands, "d0piKmMCHelic", &D0KchCand::d0pKmMCHelic, _tuple);
    ntupleVar(cands, "d0ppMCMass", &D0KchCand::d0ppMCMass, _tuple);
    ntupleVar(cands, "d0ppMCHelic", &D0KchCand::d0ppMCHelic, _tuple);
    ntupleVar(cands, "d0pPpMCMass", &D0KchCand::d0pPpMCMass, _tuple);
    ntupleVar(cands, "d0pPpMCHelic", &D0KchCand::d0pPpMCHelic, _tuple);
    ntupleVar(cands, "d0pPmMCMass", &D0KchCand::d0pPmMCMass, _tuple);
    ntupleVar(cands, "d0pPmMCHelic", &D0KchCand::d0pPmMCHelic, _tuple); 
    ntupleVarInt(cands, "d0trueKsPP", &D0KchCand::d0trueKsPP, _tuple); 
  }

  // tow body compound mass for D0 daughters
  // 1. K+K-, K K_S0 compound mass and helicity
  ntupleVar(cands, "d0kkMass", &D0KchCand::d0kkMass, _tuple);
  ntupleVar(cands, "d0kkHelic", &D0KchCand::d0kkHelic, _tuple);
  
  // 2. K+pi0, K+pi-, pi+Ks compound mass and helicity
  ntupleVar(cands, "d0piKpMass", &D0KchCand::d0pKpMass, _tuple);
  ntupleVar(cands, "d0piKpHelic", &D0KchCand::d0pKpHelic, _tuple);

  //3. K-pi0, K-pi+, pi-Ks compound mass and helicity
  ntupleVar(cands, "d0piKmMass", &D0KchCand::d0pKmMass, _tuple);
  ntupleVar(cands, "d0piKmHelic", &D0KchCand::d0pKmHelic, _tuple);

  //4. pi+pi- compound mass and helicity
  ntupleVar(cands, "d0ppMass", &D0KchCand::d0ppMass, _tuple);
  ntupleVar(cands, "d0ppHelic", &D0KchCand::d0ppHelic, _tuple);
  
  //5. pi+pi0 compound mass and helicity
  ntupleVar(cands, "d0pPpMass", &D0KchCand::d0pPpMass, _tuple);
  ntupleVar(cands, "d0pPpHelic", &D0KchCand::d0pPpHelic, _tuple);

  //6. pi-pi0 compound mass and helicity
  ntupleVar(cands, "d0pPmMass", &D0KchCand::d0pPmMass, _tuple);
  ntupleVar(cands, "d0pPmHelic", &D0KchCand::d0pPmHelic, _tuple); 


  //more D0 mass constraint
  ntupleVar(cands, "d0kkUPMass", &D0KchCand::d0kkUPMass, _tuple);
  ntupleVar(cands, "d0kkUPHelic", &D0KchCand::d0kkUPHelic, _tuple);
  ntupleVar(cands, "d0piKpUPMass", &D0KchCand::d0pKpUPMass, _tuple);
  ntupleVar(cands, "d0piKpUPHelic", &D0KchCand::d0pKpUPHelic, _tuple);
  ntupleVar(cands, "d0piKmUPMass", &D0KchCand::d0pKmUPMass, _tuple);
  ntupleVar(cands, "d0piKmUPHelic", &D0KchCand::d0pKmUPHelic, _tuple);
  ntupleVar(cands, "d0ppUPMass", &D0KchCand::d0ppUPMass, _tuple);
  ntupleVar(cands, "d0ppUPHelic", &D0KchCand::d0ppUPHelic, _tuple);
  ntupleVar(cands, "d0pPpUPMass", &D0KchCand::d0pPpUPMass, _tuple);
  ntupleVar(cands, "d0pPpUPHelic", &D0KchCand::d0pPpUPHelic, _tuple);
  ntupleVar(cands, "d0pPmUPMass", &D0KchCand::d0pPmUPMass, _tuple);
  ntupleVar(cands, "d0pPmUPHelic", &D0KchCand::d0pPmUPHelic, _tuple); 


  //7. HardTrack combined with neutral or opposite track of D0 daughters
  ntupleVar(cands, "mixNeuMass", &D0KchCand::MixNeuMass, _tuple);
  ntupleVar(cands, "mixChgMass", &D0KchCand::MixChgMass, _tuple);
  ntupleVar(cands, "D0WrongMass", &D0KchCand::D0WrongMass, _tuple);
  //D0 parameters
  //D0 vertexing fit chisq:
  ntupleVar(cands, "d0VtxFitChisq", &D0KchCand::d0VtxFitChisq, _tuple);
  ntupleVar(cands, "d0VtxDY",  &D0KchCand::d0VtxDY,  _tuple);
  ntupleVar(cands, "d0VtxDYChisq",  &D0KchCand::d0VtxDYChisq,  _tuple);
  //ntupleVar(cands, "cosPbD0Thr", &D0KchCand::cosBMomD0Thr, _tuple);

  //bachelor Kaon to D0 vtx doca and chisquare
  ntupleVar(cands, "KbD0doca", &D0KchCand::KbD0Doca, _tuple);	 
//  ntupleVar(cands, "ChisqKbD0doca", &D0KchCand::KbD0DocaChisq, _tuple);  

  //D0 decay mod, invarant mass, helicity and momentum in CM:
  ntupleVar(cands, "d0Mass", &D0KchCand::d0Mass, _tuple);
  ntupleVar(cands, "D0FkPromptMass",&D0KchCand::D0FkPromptMass, _tuple);
  ntupleVar(cands, "d0Pcm", &D0KchCand::d0CMMomtum, _tuple);
  ntupleVar(cands, "cosBDThXDFm",&D0KchCand::cosBDThXDFm, _tuple);
  ntupleVar(cands, "cosBmomDThrDFm",&D0KchCand::cosBmomDThrDFm, _tuple);
  ntupleVarInt(cands, "d0DecMode", &D0KchCand::D0DecMode, _tuple);

  if( isMC() ) {
    // to decide the good D0 candidate by mc 
    ntupleVarInt(cands, "TrueD0flg", &D0KchCand::OrigD0, _tuple); 
    // reco structed d0 dec truth with mc 
    ntupleVarInt(cands, "recD0CaDec", &D0KchCand::recD0CaDec, _tuple);
    ntupleVarInt(cands, "d0RecDec", &D0KchCand::d0RecDec, _tuple);
  }
	
  //fake D0 parameters
  //D0 -> K+K-pi0, K+ replaced by pi+
  //D0 -> pi+pi-pi0, pi+ replaced by K+
  //D0 -> K+ pi- K_S0, K+ replaced by pi+
  ntupleVar(cands, "d0FakeMasspKp", &D0KchCand::d0FkMpKp, _tuple);
 
  //D0 -> K+K-pi0, K- replaced by pi-
  //D0 -> pi+pi-pi0, pi- replaced by K-
  //D0 ->K+pi-K_S0,  pi- replaced by K-
  ntupleVar(cands, "d0FakeMasspKm", &D0KchCand::d0FkMpKm, _tuple);

  //wrong D0 mass, two tracks are swaped
  ntupleVar(cands, "d0SwpMass", &D0KchCand::d0SwpMass, _tuple);

  // B+ parameters
  //B candidate theta angle in CM frame
  ntupleVar(cands, "BCMcosTheta", &D0KchCand::BCMcosTheta, _tuple);
  
  if( isMC() ) { 
     //B decay truth
     ntupleVarInt(cands, "exclTruth", &D0KchCand::exclTruth, _tuple);
     // is a true B
     ntupleVarInt(cands, "trueBflg", &D0KchCand::truBflg, _tuple); 		
     // B decay into a Charm stat???
     ntupleVarInt(cands, "decInChrm", &D0KchCand::decInChrm, _tuple);  
  }	

  //B vertexing chisq:
  ntupleVar(cands,"bVtxFitChisq", &D0KchCand::bVtxFitChisq, _tuple);

  // the cosine of the angle between B, D vtx vector and D momentum
  ntupleVar(cands,"cosBDVtxDmom", &D0KchCand::cosBDVtxDmom, _tuple);

  if(_debug.value() && isMC() ){
     ntupleVar(cands,"TruCosBDVtxDmom", &D0KchCand::TruCosBDVtxDmom, _tuple);
  }
	  
  // Thanking:
  ntupleVar(cands,"probB0E",      &D0KchCand::probB0E,      _tuple);
  ntupleVar(cands,"probB0Mu",     &D0KchCand::probB0Mu,     _tuple);
  
  //B, D0 vertex distance chisq:
  ntupleVar(cands,"bdVtxVctChisq", &D0KchCand::bdVtxVctChisq, _tuple);
  ntupleVar(cands,"d0FlightDist", &D0KchCand::d0FlightDist, _tuple);

  //z direction of distance between B vtx and vtx of rest of Event && error
  ntupleVar(cands,"BandBRoeY", &D0KchCand::BandBRoeY, _tuple); 	
  ntupleVar(cands,"BandBRoeZ", &D0KchCand::BandBRoeZ, _tuple); 	
  ntupleVar(cands,"BandBRoeYerr", &D0KchCand::BandBRoeYerr, _tuple); 	
  ntupleVar(cands,"BandBRoeZerr", &D0KchCand::BandBRoeZerr, _tuple); 	
  if( isMC() )   ntupleVar(cands,"dzdiff", &D0KchCand::dZdf, _tuple); 	
  ntupleVar(cands,"BandBRoEnoKz", &D0KchCand::BandBRoEnoKz, _tuple);
  ntupleVar(cands,"BandBRoEnoKzErr", &D0KchCand::BandBRoEnoKzErr, _tuple);

  // energy substitute mass
  ntupleVar(cands, "mES", &D0KchCand::Mes, _tuple);
 
  // energy difference in CM
  ntupleVar(cands, "DeltaE", &D0KchCand::exclDE, _tuple);

  // pi0 variables required by Abi
  ntupleVarInt(cands, "good", &D0KchCand::good, _tuple);
  ntupleVar(cands, "mass", &D0KchCand::mass, _tuple);
  ntupleVar(cands, "helic", &D0KchCand::helic, _tuple);
  ntupleVar(cands, "mom", &D0KchCand::mom, _tuple);
  ntupleVar(cands, "massH1", &D0KchCand::massH1, _tuple);
  ntupleVar(cands, "massH2", &D0KchCand::massH2, _tuple);
  ntupleVar(cands, "massS1", &D0KchCand::massS1, _tuple);
  ntupleVar(cands, "massS2", &D0KchCand::massS2, _tuple);
  ntupleVar(cands, "helicH1", &D0KchCand::helicH1, _tuple);
  ntupleVar(cands, "helicH2", &D0KchCand::helicH2, _tuple);
  ntupleVar(cands, "helicS1", &D0KchCand::helicS1, _tuple);
  ntupleVar(cands, "helicS2", &D0KchCand::helicS2, _tuple);
  ntupleVar(cands, "asymH1", &D0KchCand::asymH1, _tuple);
  ntupleVar(cands, "asymH2", &D0KchCand::asymH2, _tuple);
  ntupleVar(cands, "asymS1", &D0KchCand::asymS1, _tuple);
  ntupleVar(cands, "asymS2", &D0KchCand::asymS2, _tuple);
  ntupleVar(cands, "momH1", &D0KchCand::momH1, _tuple);
  ntupleVar(cands, "momH2", &D0KchCand::momH2, _tuple);
  ntupleVar(cands, "momS1", &D0KchCand::momS1, _tuple);
  ntupleVar(cands, "momS2", &D0KchCand::momS2, _tuple);
  ntupleVarInt(cands, "goodH1", &D0KchCand::goodH1, _tuple);
  ntupleVarInt(cands, "goodH2", &D0KchCand::goodH2, _tuple);
  ntupleVarInt(cands, "goodS1", &D0KchCand::goodS1, _tuple);
  ntupleVarInt(cands, "goodS2", &D0KchCand::goodS2, _tuple);
}
