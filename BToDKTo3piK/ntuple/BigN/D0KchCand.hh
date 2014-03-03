//==========================================================================
// File and Version Information:
// 	$Id: D0KchCand.hh,v 1.2 2006/01/17 17:25:13 zhangjl Exp $
//
//--------------------------------------------------------------------------
// Description:
//	Class D0KchCand holds variables of a D* rho candidate,
// allowing us to sort candidates, etc., before ntupling them.
//
//--------------------------------------------------------------------------
// Collaborating classes:
//
//--------------------------------------------------------------------------
// Sample User Code:
//
//--------------------------------------------------------------------------
// Compiler Default Member Functions:
//
//--------------------------------------------------------------------------
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
//--------------------------------------------------------------------------
// Author List:
//	Qinglin Zeng              (Original author)
//
//--------------------------------------------------------------------------
// Copyright Information:
//	Copyright (C) 2000	Colorado State University
//
//==========================================================================

#ifndef D0KchCand_HH
#define D0KchCand_HH

// A macro to define a member variable and a static accessor that can be 
// passed to a function:
//#define D0KchVarDef(A, B) D0KchVar<A,0> _B; A B() const {return _B.var();}

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

//-----------------
// BaBar Headers --
//-----------------
#include "BaBar/BaBar.hh"

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//--------------------------------------------
// Collaborating Class Forward Declarations --
//--------------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

// All members of this class are initialized to the value I:
template <class T, int I>
class D0KchVar {
public:
  D0KchVar() {_var = I;}
  virtual ~D0KchVar() {}
  D0KchVar & operator=(const T t) {_var = t; return *this;}
  T var() const {return _var;}

private:
  T _var;
};


class D0KchCand {

public:
  // Constructors
  D0KchCand();

  // Destructor
  virtual ~D0KchCand();

  // operators:
  bool operator==(const D0KchCand &) const {return false;}

  // Data and accessors:

  //Fisher Discriminant
  D0KchVar<float, -99>  _BFisher;
  float BFisher() const { return _BFisher.var(); }

  //CLEO fisher
  D0KchVar<float, -99>  _BFCLeo;
  float BFCLeo() const { return _BFCLeo.var(); }

  D0KchVar<float, -99>  _BFCLeo2;
  float BFCLeo2() const { return _BFCLeo2.var(); }

  D0KchVar<float, -99>  _BFLegendre;
  float BFLegendre() const { return _BFLegendre.var(); }
 
  //B candidate thrust, only a shape parameter
  D0KchVar<float, -99> _BCanThrust;
  float BCanThrust() const { return _BCanThrust.var(); }

  //cosine of the angle between the thrust axis of the B
  //candidate and the thrust axis of the rest of the event
  D0KchVar<float, -99> _BCanCosThR;
  float BCanCosThR() const { return _BCanCosThR.var(); }

  //cosine of the angle between the B direction and the beam axis
  D0KchVar<float, -99> _cosBmomAx;
  float cosBmomAx() const { return _cosBmomAx.var(); }

  //cosine of the angle between the thrust axis of the B candidate
  //and the beam axis
  D0KchVar<float, -99> _cosBThRAx;
  float cosBThRAx() const { return _cosBThRAx.var(); }

  //the angle between B mom and D0 thrust axis in CM
  D0KchVar<float, -99> _cosBMomD0Thr;
  float cosBMomD0Thr() const { return _cosBMomD0Thr.var(); }

  //the z distance between B vtx and rest of Event Vtx
  D0KchVar<float, -1000> _BandBRoeX;
  float BandBRoeX() const { return _BandBRoeX.var(); }
 
  D0KchVar<float, -1000> _BandBRoeXerr;
  float BandBRoeXerr() const { return _BandBRoeXerr.var(); }

  D0KchVar<float, -1000> _BandBRoeY;
  float BandBRoeY() const { return _BandBRoeY.var(); }

  D0KchVar<float, -1000> _BandBRoeYerr;
  float BandBRoeYerr() const { return _BandBRoeYerr.var(); }

  D0KchVar<float, -1000> _BandBRoeZ;
  float BandBRoeZ() const { return _BandBRoeZ.var(); }

  D0KchVar<float, -1000> _BandBRoeZerr;
  float BandBRoeZerr() const { return _BandBRoeZerr.var(); }

  D0KchVar<float, -1000> _dZdf;
  float dZdf() const { return _dZdf.var(); }
  
  //the z dir distance between beam spot and vtx of K bach and non k trk
  D0KchVar<float, -1000> _BandBRoEnoKzErr;
  float BandBRoEnoKzErr() const { return _BandBRoEnoKzErr.var(); }

  D0KchVar<float, -1000> _BandBRoEnoKz;
  float BandBRoEnoKz() const { return _BandBRoEnoKz.var(); }
 
  //Hemisphere charge diff
  D0KchVar<float, -100> _QDiff;
  float QDiff() const { return _QDiff.var(); }	

  // HardTrack quantities:
  D0KchVar<int,0> _HdTrkChg; // the charge of the hard track
  int HdTrkChg() const { return _HdTrkChg.var(); }

  D0KchVar<float,-999> _HdTrkTheta; // theta angle of the track
  float HdTrkTheta() const { return _HdTrkTheta.var(); }	

  D0KchVar<float,-999> _HdTrkCosTheta; // cosine of theta angle
  float HdTrkCosTheta() const {return _HdTrkCosTheta.var();}

  D0KchVar<float,-999> _HdTrkPhi;  // phi angle
  float HdTrkPhi() const {return _HdTrkPhi.var();}

  D0KchVar<float,-999> _HdTrkMom;  // lab momentum  
  float HdTrkMom() const {return _HdTrkMom.var();}
  
  D0KchVar<float,-999> _HdTrkThetaC;// cherenkov angle
  float HdTrkThetaC() const {return _HdTrkThetaC.var();}

  D0KchVar<float,-999> _HdTrkThcErr; // cherenkov angle err 
  float HdTrkThcErr() const { return _HdTrkThcErr.var();}

  D0KchVar<float,-99> _HdTrkSigOffProt; //kaon hypo off prot hypo;
  float HdTrkSigOffProt() const { return _HdTrkSigOffProt.var(); }

  D0KchVar<float, 0> _HdTrkNPhot; //dirc phot num for kaon;
  float HdTrkNPhot() const { return _HdTrkNPhot.var(); }

  D0KchVar<float, 0> _HdTrkNBkPhot; //dirc phot num for kaon;
  float HdTrkNBkPhot() const { return _HdTrkNBkPhot.var(); }; 
  
  D0KchVar<float, 0> _HdTrkKaonNN; //dirc phot num for kaon;
  float HdTrkKaonNN() const { return _HdTrkKaonNN.var(); };

  D0KchVar<int, 0> _HdTrkPid;
  int HdTrkPid() const { return _HdTrkPid.var(); }

  D0KchVar<int, 0> _HdTrkElecPid; //HdTrk is a  elec, default false
  int HdTrkElecPid() const { return _HdTrkElecPid.var(); }

  D0KchVar<int, 0> _HdTrkMuonPid; //HdTrk is a  muon, default is false
  int HdTrkMuonPid() const { return _HdTrkMuonPid.var(); }

  D0KchVar<float,-999> _HdTrkCMTheta; // cosine of theta in CM
  float HdTrkCMTheta() const { return _HdTrkCMTheta.var(); }

  D0KchVar<float,-999> _HdTrkCMCosTheta; // cosine of theta in CM
  float HdTrkCMCosTheta() const { return _HdTrkCMCosTheta.var(); }

  D0KchVar<float,-999> _HdTrkCMPhi;
  float HdTrkCMPhi() const { return _HdTrkCMPhi.var(); }
 
  D0KchVar<float,-999> _HdTrkCMMomtum; // momentum in CM
  float HdTrkCMMomtum() const {return _HdTrkCMMomtum.var();}

  D0KchVar<int,-1> _HdTrkTruth; // truth of Hard Track
  int HdTrkTruth() const { return _HdTrkTruth.var(); }

  D0KchVar<int, -1> _TotBDau;
  int TotBDau() const { return _TotBDau.var() ; }
	
  D0KchVar<int, -1> _TotChrmDau[2];
  int TotChrmDau1() const { return _TotChrmDau[0].var(); }
  int TotChrmDau2() const { return _TotChrmDau[1].var(); }

  D0KchVar<int, -1> _ChrmDau[2][4];
  int ChrmDau11() const { return _ChrmDau[0][0].var(); }
  int ChrmDau12() const { return _ChrmDau[0][1].var(); }
  int ChrmDau13() const { return _ChrmDau[0][2].var(); }
  int ChrmDau14() const { return _ChrmDau[0][3].var(); }
  int ChrmDau21() const { return _ChrmDau[1][0].var(); }
  int ChrmDau22() const { return _ChrmDau[1][1].var(); }
  int ChrmDau23() const { return _ChrmDau[1][2].var(); }
  int ChrmDau24() const { return _ChrmDau[1][3].var(); } 	 

  D0KchVar<int, -1> _BDau[6];
  int BDau1() const { return _BDau[0].var(); }
  int BDau2() const { return _BDau[1].var(); }
  int BDau3() const { return _BDau[2].var(); }
  int BDau4() const { return _BDau[3].var(); }
  int BDau5() const { return _BDau[4].var(); }
  int BDau6() const { return _BDau[5].var(); }

  // bachelor kaon with another kaon from D0 and its trk vtx	
  D0KchVar<float, -999> _KbKUpMass;
  float KbKUpMass() const { return _KbKUpMass.var(); }
 
  D0KchVar<float, -999> _KbKLowMass;
  float KbKLowMass() const { return _KbKLowMass.var(); }

  D0KchVar<int, -10> _nRoEKp;
  int nRoEKp() const { return _nRoEKp.var(); }

  D0KchVar<int, -10> _nRoEKm;
  int nRoEKm() const { return _nRoEKm.var(); }

  D0KchVar<float, -999> _KbKsUpMass;
  float KbKsUpMass() const { return _KbKsUpMass.var(); }

  D0KchVar<float, -999> _KbKsLowMass;
  float KbKsLowMass() const { return _KbKLowMass.var(); }

  D0KchVar<int, -10> _nRoEKs;
  int nRoEKs() const { return _nRoEKs.var(); }

  D0KchVar<float, -999> _ChisqKbKUpVtx;
  float ChisqKbKUpVtx() const { return _ChisqKbKUpVtx.var(); }

  D0KchVar<float, -999> _ChisqKbKLowVtx;
  float ChisqKbKLowVtx() const { return _ChisqKbKLowVtx.var(); }

  D0KchVar<float, -999> _ChisqKbLepVtx;
  float ChisqKbLepVtx() const { return _ChisqKbLepVtx.var(); }

  //vtx of bachelor Kaon with a trk from ROE and its distance to beam spot in y
  D0KchVar<float, -999> _KbRoETrkChisq;
  float KbRoETrkChisq() const { return _KbRoETrkChisq.var(); }

  D0KchVar<float, -999> _VtxKbRoETrkDistY;
  float VtxKbRoETrkDistY() const { return _VtxKbRoETrkDistY.var(); }

  D0KchVar<float, -999> _VtxKbRoETrkDistChisqY;
  float VtxKbRoETrkDistChisqY() const { return _VtxKbRoETrkDistChisqY.var(); }		      	

  //doca of a Kaon from ROE to B(signal) vtx
  D0KchVar<float, -999> _KroeBvtxDoca;
  float KroeBvtxDoca() const { return _KroeBvtxDoca.var(); }

  // dist between vtx of Kaon of ROE and beamspot and B(signal) vtx
  D0KchVar<float, -999> _DistKroeBsptBvtx;
  float DistKroeBsptBvtx() const { return _DistKroeBsptBvtx.var(); }
 	
  //doca of a trk(no ks) from ROE to B(signal) vtx
  D0KchVar<float, -999> _TrkRoeBvtxDoca;
  float TrkRoeBvtxDoca() const { return _TrkRoeBvtxDoca.var(); }

  // dist between vtx of a trk(no ks) of ROE and beamspot and B(signal) vtx
  D0KchVar<float, -999> _DistRoeTrknoKsBsptBvtx;
  float DistRoeTrknoKsBsptBvtx() const { return _DistRoeTrknoKsBsptBvtx.var(); }

  // bachelor kaon with a lepton
  D0KchVar<float, -999> _KbLepMass;
  float KbLepMass() const { return _KbLepMass.var(); }

  //doca of bachelor kaon to D0 vtx
  D0KchVar<float, -1000> _KbD0Doca;
  float KbD0Doca() const { return _KbD0Doca.var(); }

  D0KchVar<float, -1000> _KbD0DocaChisq;
  float KbD0DocaChisq() const { return _KbD0DocaChisq.var(); }
	
  // K_S0 mass 
  D0KchVar<float,0> _KsMass;
  float KsMass() const { return _KsMass.var();}
 
  D0KchVar<float,-100> _KsHel;
  float KsHel() const { return _KsHel.var();}
 
  // K_S0 fit vtx chisquare:
  D0KchVar<float,-999> _KsVtxFitChisq;
  float KsVtxFitChisq() const { return _KsVtxFitChisq.var();}
 
 
  // origin of Ks 1. K+ pi- Ks, 2. K- pi+ Ks
  D0KchVar<int, -1> _OrigKs;
  int OrigKs() const { return _OrigKs.var(); }

  // K_S0 decay distance
  D0KchVar<float,-1> _KsDecLen;
  float KsDecLen() const { return _KsDecLen.var();}

  // K_S0 angle
  D0KchVar<float,-100> _KsTheta;
  float KsTheta() const { return _KsTheta.var(); }
 
  // Ks momentum
  D0KchVar<float, -100> _pKs;
  float pKs() const { return _pKs.var(); }
 
  // K_S0 decay length in chisq expression:
  D0KchVar<float, -1> _dksVtxVctChisq;
  float dksVtxVctChisq() const { return _dksVtxVctChisq.var(); }
 
  // cosine of the angle between the vector link decay point 
  // of ks and D0 and the momentum of Ks.
  D0KchVar<float,-99> _cosDKsVtx;
  float cosDKsVtx() const { return _cosDKsVtx.var(); }

  // pi0 Mass
  D0KchVar<float,0> _pi0Mass;
  float pi0Mass() const { return _pi0Mass.var();}

  //pi0 origin 1. K+ K- pi0, 2. pi+ pi- pi0 
  D0KchVar<int, -1> _OrigPi0;
  int OrigPi0() const { return _OrigPi0.var(); }  

  //pi+ pi- track
  D0KchVar<int, -10> _pipTru; 
  int pipTru() const { return _pipTru.var(); }

  D0KchVar<int, -10> _pimTru;
  int pimTru() const { return _pimTru.var(); }
  
  // gamma energy 
  D0KchVar<float,0> _gamEng1;
  float gamEng1() const { return _gamEng1.var();}

  D0KchVar<float,0> _gamEng2;
  float gamEng2() const { return _gamEng2.var();}

  //gamma energy from MC
  D0KchVar<float, 0> _truPi0E1;
  float truPi0E1() const { return _truPi0E1.var(); }

  D0KchVar<float, 0> _truPi0E2; 
  float truPi0E2() const { return _truPi0E2.var(); }	

  //gamma Latency variables
  D0KchVar<float,-1> _gamLat1;
  float gamLat1() const { return _gamLat1.var(); }

  D0KchVar<float, -1> _gamLat2;
  float gamLat2() const { return _gamLat2.var(); }

  //calerimeter bum raw energy for two gamma shower;
  D0KchVar<float, -1> _gamBumRawE1;
  float gamBumRawE1() const { return _gamBumRawE1.var(); }

  D0KchVar<float, -1> _gamBumRawE2;
  float gamBumRawE2() const { return _gamBumRawE2.var(); }

  D0KchVar<float, -999> _pi0Hel;
  float pi0Hel() const { return _pi0Hel.var(); }

  D0KchVar<float, -999> _pi0HelDfm;
  float pi0HelDfm() const { return _pi0HelDfm.var(); }

  D0KchVar<float, -10> _pi0CMmom0;
  float pi0CMmom0() const { return _pi0CMmom0.var(); }

  D0KchVar<float, -10> _pi0CMmom;
  float pi0CMmom() const { return _pi0CMmom.var(); }

  D0KchVar<float, -10> _pi0EAsy;
  float pi0EAsy() const { return _pi0EAsy.var(); }

  //one of photon from pi0 combines with any other photon from the event	
  D0KchVar<float, -100> _pi0BestMass1;
  float pi0BestMass1() const { return _pi0BestMass1.var(); }
  D0KchVar<float, -100> _pi0BestMass2;
  float pi0BestMass2() const { return _pi0BestMass2.var(); }
  
  D0KchVar<float, -100> _Bestcospi0Hel1;
  float Bestcospi0Hel1() const { return _Bestcospi0Hel1.var(); }

  D0KchVar<float, -100> _Bestcospi0Hel2;
  float Bestcospi0Hel2() const { return _Bestcospi0Hel2.var(); }

  D0KchVar<float, -100> _Best2gPcm1;
  float Best2gPcm1() const { return _Best2gPcm1.var(); }
  
  D0KchVar<float, -100> _Best2gPcm2;
  float Best2gPcm2() const { return _Best2gPcm2.var(); }
  //trk kaon pid info
  D0KchVar<int,0> _KpPidBit; //default is false
  int KpPidBit() const { return _KpPidBit.var(); }

  D0KchVar<int,0> _KmPidBit; //default is false
  int KmPidBit() const { return _KmPidBit.var(); }

  //trk pion pid info
  D0KchVar<int,0> _PipPidBit; //default is false
  int PipPidBit() const { return _PipPidBit.var(); }

  D0KchVar<int,0> _PimPidBit; //default is false
  int PimPidBit() const { return _PimPidBit.var(); }

  //trk geometry info
  D0KchVar<float, -999> _D0DauTrkP1;
  float D0DauTrkP1() const { return _D0DauTrkP1.var(); }

  D0KchVar<float, -999> _D0DauTrkP2;
  float D0DauTrkP2() const { return _D0DauTrkP2.var(); }

  D0KchVar<float, -999> _D0DauTrkTheta1;
  float D0DauTrkTheta1() const { return _D0DauTrkTheta1.var(); }
  
  D0KchVar<float, -999> _D0DauTrkTheta2;
  float D0DauTrkTheta2() const { return _D0DauTrkTheta2.var(); }

  D0KchVar<float, -999> _D0DauTrkPhi1;
  float D0DauTrkPhi1() const { return _D0DauTrkPhi1.var(); }

  D0KchVar<float, -999> _D0DauTrkPhi2;
  float D0DauTrkPhi2() const { return _D0DauTrkPhi2.var(); }
  // D0 decay mode:
  D0KchVar<int, 0> _D0DecMode; 
  // D0 decay mode 
  //1. D0 -> K+ K-  pi0;  2. D0 -> K+  pi-  K_S0 
  //3. D0 -> K- pi+ K_S0 4. D0 -> pi+ pi- pi0
  int D0DecMode() const { return _D0DecMode.var(); }

  // D0 Mass
  D0KchVar<float,0> _d0Mass;
  float d0Mass() const { return _d0Mass.var();}

  // D0 truth flag
  D0KchVar<int, -1> _OrigD0; //1. good, 0 bad
  int OrigD0() const { return _OrigD0.var(); }

  //all rec D0 dec mode
  D0KchVar<int,0> _d0RecDec;
  int d0RecDec() const { return _d0RecDec.var(); }
  //rec D0 dec Cabibbo mode
  D0KchVar<int, 0> _recD0CaDec;
  int recD0CaDec() const { return _recD0CaDec.var(); }

  //d0 vertexing fit chisqare
  D0KchVar<float, -1> _d0VtxFitChisq;
  float d0VtxFitChisq() const { return _d0VtxFitChisq.var(); }

  //d0 vertexing fit chisq with mass constraint --JZ
  D0KchVar<float, -1> _d0VtxFitMSChisq;
  float d0VtxFitMSChisq() const { return _d0VtxFitMSChisq.var(); }

  //d0 flight Distance(from B)
  D0KchVar<float, -1> _d0FlightDist;
  float d0FlightDist() const { return _d0FlightDist.var(); }

  //d0 vtx to beam spot distance and chisquare in y direction
  D0KchVar<float, -100> _d0VtxDY;
  float d0VtxDY() const { return _d0VtxDY.var(); }

  D0KchVar<float, -100> _d0VtxDYChisq;
  float d0VtxDYChisq() const { return _d0VtxDYChisq.var(); }
 
  //fake D0 Mass, D0 daughter pi is replaced by the prompt track K
  D0KchVar<float,0> _D0FkPromptMass;
  float D0FkPromptMass() const { return _D0FkPromptMass.var(); }

  // fake D0 Mass, D0 daughter K/pi is replaced by a pi/k 
  D0KchVar<float,0> _d0FkMpKm;
  //D0->K+K-K_S0,  K- replaced by pi-
  //D0->pi+pi-pi0, pi- replaced by K- 
  //D0->K+pi-pi0,  pi- replaced by K-
  float d0FkMpKm() const { return _d0FkMpKm.var();}

  D0KchVar<float,0> _d0FkMpKp;
  //D0->K+ K- pi0,   K+ repalced by pi+
  //D0->K+ pi- K_S0, K+ replaced by pi+
  //D0->pi+ pi- pi0 , pi+ replaced by K+
  float d0FkMpKp() const { return _d0FkMpKp.var();}

  //fake D0 mass, two charge trks are swaped
  D0KchVar<float,0> _d0SwpMass;
  float d0SwpMass() const { return _d0SwpMass.var(); }

  //compound mass:
  //two body compound mass D0 -> K+ K- pi0
  // two body compound mass D0 -> K pi Ks
  //truth  
  D0KchVar<float,0> _d0kkMCMass; //K+K-, KKs
  float d0kkMCMass() const { return _d0kkMCMass.var();}

  D0KchVar<float,0> _d0pKpMCMass; //K+pi0,K+pi-,pi+Ks
  float d0pKpMCMass() const { return _d0pKpMCMass.var();}

  D0KchVar<float,0> _d0pKmMCMass; //K-pi0,K-pi+,pi-Ks
  float d0pKmMCMass() const { return _d0pKmMCMass.var();}

  D0KchVar<float,0> _d0kkMass; //K+K-, KKs
  float d0kkMass() const { return _d0kkMass.var();}

  D0KchVar<float,0> _d0pKpMass; //K+pi0,K+pi-,pi+Ks
  float d0pKpMass() const { return _d0pKpMass.var();}

  D0KchVar<float,0> _d0pKmMass; //K-pi0,K-pi+,pi-Ks
  float d0pKmMass() const { return _d0pKmMass.var();}

  //compond mass with D0 mass constraint
  D0KchVar<float,0> _d0kkUPMass; //K+K-, KKs
  float d0kkUPMass() const { return _d0kkUPMass.var();}

  D0KchVar<float,0> _d0pKpUPMass; //K+pi0,K+pi-,pi+Ks
  float d0pKpUPMass() const { return _d0pKpUPMass.var();}

  D0KchVar<float,0> _d0pKmUPMass; //K-pi0,K-pi+,pi-Ks
  float d0pKmUPMass() const { return _d0pKmUPMass.var();}

  D0KchVar<float, -100> _pipiDecLen;
  float pipiDecLen() const {return _pipiDecLen.var(); }

  D0KchVar<float, -100> _pipiMass;
  float pipiMass() const {return _pipiMass.var(); }

  D0KchVar<int, 0> _pipiKsTrue;
  int pipiKsTrue() const {return _pipiKsTrue.var(); }

  //hdtrk combined with a neutral or opposite charge trk of D0 daughters
  D0KchVar<float,0> _MixNeuMass;
  float MixNeuMass() const { return _MixNeuMass.var(); }

  D0KchVar<float,0> _MixChgMass;
  float MixChgMass() const { return _MixChgMass.var(); }

  D0KchVar<float,0> _D0WrongMass;
  float D0WrongMass() const { return _D0WrongMass.var(); }

  // tow body compound helicity D0 -> K+ K- pi0
  // two body compound helicity D0 -> K pi Ks
  //truth
  D0KchVar<float,-999> _d0kkMCHelic; //K+K-, KKs
  float d0kkMCHelic() const { return _d0kkMCHelic.var();}

  D0KchVar<float,-999> _d0pKpMCHelic; //K+pi0,K+pi-,pi+Ks
  float d0pKpMCHelic() const { return _d0pKpMCHelic.var();}

  D0KchVar<float,-999> _d0pKmMCHelic; //K-pi0,K-pi+,pi-Ks
  float d0pKmMCHelic() const { return _d0pKmMCHelic.var();}

  D0KchVar<float,-999> _d0kkHelic; //K+K-, KKs
  float d0kkHelic() const { return _d0kkHelic.var();}

  D0KchVar<float,-999> _d0pKpHelic; //K+pi0,K+pi-,pi+Ks
  float d0pKpHelic() const { return _d0pKpHelic.var();}

  D0KchVar<float,-999> _d0pKmHelic; //K-pi0,K-pi+,pi-Ks
  float d0pKmHelic() const { return _d0pKmHelic.var();}

  // D0 mass constraint
  D0KchVar<float,-999> _d0kkUPHelic; //K+K-, KKs
  float d0kkUPHelic() const { return _d0kkUPHelic.var();}

  D0KchVar<float,-999> _d0pKpUPHelic; //K+pi0,K+pi-,pi+Ks
  float d0pKpUPHelic() const { return _d0pKpUPHelic.var();}

  D0KchVar<float,-999> _d0pKmUPHelic; //K-pi0,K-pi+,pi-Ks
  float d0pKmUPHelic() const { return _d0pKmUPHelic.var();}

  // two body compound mass D0 -> pi+ pi- pi0
  //KS or not 
  D0KchVar<int, -1> _d0trueKsPP;
  int d0trueKsPP() const { return _d0trueKsPP.var();}

  //truth
  D0KchVar<float,0> _d0ppMCMass;
  float d0ppMCMass() const { return _d0ppMCMass.var();}

  D0KchVar<float,0> _d0pPmMCMass;
  float d0pPmMCMass() const { return _d0pPmMCMass.var();}

  D0KchVar<float,0> _d0pPpMCMass;
  float d0pPpMCMass() const { return _d0pPpMCMass.var();}

  D0KchVar<float,0> _d0ppMass;
  float d0ppMass() const { return _d0ppMass.var();}

  D0KchVar<float,0> _d0pPmMass;
  float d0pPmMass() const { return _d0pPmMass.var();}

  D0KchVar<float,0> _d0pPpMass;
  float d0pPpMass() const { return _d0pPpMass.var();}

  //D0 mass constraint
  D0KchVar<float,0> _d0ppUPMass;
  float d0ppUPMass() const { return _d0ppUPMass.var();}

  D0KchVar<float,0> _d0pPmUPMass;
  float d0pPmUPMass() const { return _d0pPmUPMass.var();}

  D0KchVar<float,0> _d0pPpUPMass;
  float d0pPpUPMass() const { return _d0pPpUPMass.var();}


  // two body compound Helic D0 -> pi+ pi- pi0
  //truth
  D0KchVar<float,-999> _d0ppMCHelic;
  float d0ppMCHelic() const { return _d0ppMCHelic.var();}

  D0KchVar<float,-999> _d0pPmMCHelic;
  float d0pPmMCHelic() const { return _d0pPmMCHelic.var();}

  D0KchVar<float,-999> _d0pPpMCHelic;
  float d0pPpMCHelic() const { return _d0pPpMCHelic.var();}

  D0KchVar<float,-999> _d0ppHelic;
  float d0ppHelic() const { return _d0ppHelic.var();}

  D0KchVar<float,-999> _d0pPmHelic;
  float d0pPmHelic() const { return _d0pPmHelic.var();}

  D0KchVar<float,-999> _d0pPpHelic;
  float d0pPpHelic() const { return _d0pPpHelic.var();}


  //D0 mass constraint
  D0KchVar<float,-999> _d0ppUPHelic;
  float d0ppUPHelic() const { return _d0ppUPHelic.var();}

  D0KchVar<float,-999> _d0pPmUPHelic;
  float d0pPmUPHelic() const { return _d0pPmUPHelic.var();}

  D0KchVar<float,-999> _d0pPpUPHelic;
  float d0pPpUPHelic() const { return _d0pPpUPHelic.var();}


  //D0 momentum
  D0KchVar<float,-999> _d0CMMomtum;
  float d0CMMomtum() const { return _d0CMMomtum.var();}

  //D0 helic
  D0KchVar<float,-999> _D0Helic; //D0 and Kch helicity
  float D0Helic() const { return _D0Helic.var(); }

  //cosine of the angle between the thrust axis/Mom of B and the thrust of 
  //D0 daughters in D0 CMS
  D0KchVar<float, -999> _cosBDThXDFm;
  float cosBDThXDFm() const { return _cosBDThXDFm.var(); }
  
  D0KchVar<float, -999> _cosBmomDThrDFm;
  float cosBmomDThrDFm() const { return _cosBmomDThrDFm.var(); }

  //wrong D0 helicity B->D0 Kch, Kch replaced by pi or rho
  D0KchVar<float,-999> _D0fakePiHelic;  
  float D0fakePiHelic() const { return _D0fakePiHelic.var();}
     
  //B, pep parameters, Mes & DeltaE
  D0KchVar<float,-999> _Mes;
  float Mes() const { return _Mes.var();}

  D0KchVar<float,-999> _exclDE;
  float exclDE() const { return _exclDE.var();}

  D0KchVar<float,-999> _hMes;
  float hMes() const { return _hMes.var();}

  D0KchVar<float,-999> _exclhDE;
  float exclhDE() const { return _exclhDE.var();}
  
  D0KchVar<float, -999> _mBhat;
  float mBhat() const { return _mBhat.var(); }
  
  D0KchVar<float,-999> _mBhatPull;
  float mBhatPull() const { return _mBhatPull.var(); }

  D0KchVar<float, -999> _deltaEPull;
  float deltaEPull() const { return _deltaEPull.var(); }

  D0KchVar<float, -999> _BCMcosTheta;
  float BCMcosTheta() const { return _BCMcosTheta.var(); }

  D0KchVar<int, -1> _exclTruth;
  int exclTruth() const { return _exclTruth.var(); }		  

  // is a true B
  D0KchVar<int, -100> _truBflg;
  int truBflg() const { return _truBflg.var(); }

  //B decay into a D charm state?
  D0KchVar<int, -10> _decInChrm;
  int decInChrm() const { return _decInChrm.var(); }
 
  //B, D0 vertex point distance chisquare:
  D0KchVar<float, -1> _bdVtxVctChisq;
  float bdVtxVctChisq() const { return _bdVtxVctChisq.var(); }
 
  // D B vtx and d mom angle cosine value
  D0KchVar<float, -100> _cosBDVtxDmom;
  float cosBDVtxDmom() const { return _cosBDVtxDmom.var(); }

  D0KchVar<float, -100> _TruCosBDVtxDmom;
  float TruCosBDVtxDmom() const { return _TruCosBDVtxDmom.var(); }

  //B vertexing fit chisquare
  D0KchVar<float, -1> _bVtxFitChisq;
  float bVtxFitChisq() const { return _bVtxFitChisq.var(); }

  //slow pion tag
  D0KchVar<float, -100> _probB0SlowPi;
  float probB0SlowPi() const { return _probB0SlowPi.var(); }

  //lepton tag
  D0KchVar<float, -100> _probB0E;
  float probB0E() const { return _probB0E.var(); }

  D0KchVar<float, -100> _probB0Mu;
  float probB0Mu() const { return _probB0Mu.var(); }
 
  D0KchVar<float, -100> _probB0KinLep;
  float probB0KinLep() const { return _probB0KinLep.var(); }

  //pi0 variables required by Abi 
    // truth:
  D0KchVar<int, -99>  _good;
  int good() const { return _good.var(); }

  // pi0 mass:
  D0KchVar<float, -99>  _mass;
  float mass() const { return _mass.var(); }

  // helicity:
  D0KchVar<float, -99>  _helic;
  float helic() const { return _helic.var(); }

  // Energy asymmetry:
  D0KchVar<float, -99>  _asym;
  float asym() const { return _asym.var(); }

  // CM momentum:
  D0KchVar<float, -99>  _mom;
  float mom() const { return _mom.var(); }

  // Best mass:
  D0KchVar<float, -99>  _massH[2];
  float massH1() const { return _massH[0].var(); }
  float massH2() const { return _massH[1].var(); }

  D0KchVar<float, -99>  _massS[2];
  float massS1() const { return _massS[0].var(); }
  float massS2() const { return _massS[1].var(); }

  // Best mass helic:
  D0KchVar<float, -99>  _helicH[2];
  float helicH1() const { return _helicH[0].var(); }
  float helicH2() const { return _helicH[1].var(); }

  D0KchVar<float, -99>  _helicS[2];
  float helicS1() const { return _helicS[0].var(); }
  float helicS2() const { return _helicS[1].var(); }

  // Best mass asym:
  D0KchVar<float, -99>  _asymH[2];
  float asymH1() const { return _asymH[0].var(); }
  float asymH2() const { return _asymH[1].var(); }

  D0KchVar<float, -99>  _asymS[2];
  float asymS1() const { return _asymS[0].var(); }
  float asymS2() const { return _asymS[1].var(); }

  // Best mass pair momentum:
  D0KchVar<float, -99>  _momH[2];
  float momH1() const { return _momH[0].var(); }
  float momH2() const { return _momH[1].var(); }

  D0KchVar<float, -99>  _momS[2];
  float momS1() const { return _momS[0].var(); }
  float momS2() const { return _momS[1].var(); }

  // Best mass truth:
  D0KchVar<int, -99>  _goodH[2];
  int goodH1() const { return _goodH[0].var(); }
  int goodH2() const { return _goodH[1].var(); }

  D0KchVar<int, -99>  _goodS[2];
  int goodS1() const { return _goodS[0].var(); }
  int goodS2() const { return _goodS[1].var(); }

 };

#endif








