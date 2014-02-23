//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: GqaMCAnalysis.cc,v 1.2 2006/05/03 20:52:19 abi Exp $
//
// Description:
//      Module to work in the BaBar software framework for
//      analysing MC truth information to study the performance 
//      of the generators.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Anders Ryd                    Original Author
//
// Copyright Information:
//	Copyright (C)
//
//------------------------------------------------------------------------
//--------------
// BaBar Header:
//--------------
#include "BaBar/BaBar.hh"

//-------------
// C Headers --
//-------------
#include <assert.h>

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <iomanip>
#include <strstream>
#include <string>
#include <vector>
#include <math.h>

//-----------------------
// This Class's Header --
//-----------------------

#include "GeneratorsQA/GqaMCAnalysis.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Alist/ConstAList.h"
#include "CLHEP/Alist/ConstAIterator.h"
#include "HepTuple/TupleManager.h"
#include "HepTuple/Tuple.h"
#include "HepTuple/Histogram.h"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/getTmpAList.hh"
#include "CLHEP/Alist/AList.h"
#include "ProxyDict/IfdStrKey.hh"
#include "ProxyDict/IfdKey.hh"
#include "ProxyDict/Ifd.hh"
#include "AbsEnv/AbsEnv.hh"
#include "GenEnv/GenEnv.hh"
#include "PDT/Pdt.hh"
#include "PDT/PdtEntry.hh"

#include "PepEnv/PepEnv.hh"
#include "PepEnvData/PepBeams.hh"
#include "PepData/PepCollision.hh"

#include "Beta/EventInfo.hh"
#include "Beta/BtaCandidate.hh"
#include "Beta/BtaAbsVertex.hh"

#include "ErrLogger/ErrLog.hh"
#include "BaBar/Constants.hh"
using std::cout;
using std::endl;
using std::istrstream;
using std::ostream;
using std::setprecision;
using std::setw;
using std::strstream;


//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//----------------
// Constructors --
//----------------

// in general, a module constructor should not do much.  The beginJob or
// beginRun members are better places to put initialization
GqaMCAnalysis::GqaMCAnalysis(const char* const theName, 
			     const char* const theDescription,
			     PdtLund::LundType DFlavor) : 
  AppModule( theName, theDescription ),
  _eventInfoList("eventInfoList",this,"Default"),
  _btaTruthList("truthCandidates",this,"MCTruth"),
  _applyBoost("applyBoost", this, true),
  _pepCollisionKey("pepCollisionKey", this, "Default"),
  _makeNtuple("makeNtuple", this, false),
  _checkMomenta("checkMomenta", this, false),
  _massTolerance("massTolerance", this, 1.e-06),
  _printCands("printCands", this, ""),
  _printTerminals("printTerminals", this, "pi0 eta"),
  _numPrintEvents("numPrintEvents", this, 10),
  _evtNum(0),
  _printTree(BtaPrintTree::printName, BtaPrintTree::printP4),
  _DFlavor(DFlavor)
{
  commands()->append(&_eventInfoList);
  commands()->append(&_btaTruthList);
  commands()->append(&_applyBoost);
  commands()->append(&_pepCollisionKey);
  commands()->append(&_makeNtuple);
  commands()->append(&_checkMomenta);
  commands()->append(&_massTolerance);
  commands()->append(&_printCands);
  commands()->append(&_printTerminals);
  commands()->append(&_numPrintEvents);
}

//--------------
// Destructor --
//--------------

// The destructor should be limited to undoing the work of the constructor
GqaMCAnalysis::~GqaMCAnalysis( )
{
}

//--------------
// Operations --
//--------------

// The beginJob member function is run before any events are
// processed.  In this analysis, it opens the output histogram file
// and then books a number of histograms and a ntuple.

AppResult
GqaMCAnalysis::beginJob( AbsEvent* anEvent )
{
  ErrMsg(routine) << name()<< " beginJob "<< endmsg; 
  
  // Add terminal particles to the BtaPrintTree:
  istrstream terminalsStream(_printTerminals.value().c_str());
  std::string terminal;
  while (terminalsStream >> terminal){
    const PdtEntry * entry = Pdt::lookup(terminal);
    if (entry != 0){
      _printTree.addTerminal(entry->lundId());
    }
  }
  
  HepTupleManager* manager = gblEnv->getGen()->ntupleManager();
  assert(manager != 0);
  
  // book histograms  (arguments are nBins, low edge, high edge,
  //                  the PAW histogram number)
  float lower = 0.0;
  float upper = 3.0;
  int channels = 60;
  
  float upperd = 5.0;
  int channelsd = 100; 
  
  _tuple = manager->ntuple(name());
  
  // Get the average beam vertex to offset the vertex plots:
  if (true == _verbose.value()){
    ErrMsg(routine) << name() << " Beams = " << endl 
		    << *(gblEnv->getPep()->pepBeams()) << endl;
  }  
  
  const HepLorentzVector avgVtx = 
    gblEnv->getPep()->pepBeams()->interactionPoint();
  
  _collisionX = manager->histogram("collision xVtx", channels, 
				   avgVtx.x() - 0.15, avgVtx.x() + 0.15, 981); 
  _collisionY = manager->histogram("collision yVtx", channels, 
				   avgVtx.y() -0.003, avgVtx.y() + 0.003, 982);
  _collisionZ = manager->histogram("collision zVtx", channels, 
				   avgVtx.z() - 4.0, avgVtx.z() + 4.0, 983 ); 
  
  _collisionPX  = manager->histogram("collision Px", 60, -0.15, 0.05, 985); 
  _collisionPY  = manager->histogram("collision Py", 60, -0.05, 0.05, 986); 
  _collisionPZ  = manager->histogram("collision Pz", 60,  5.8, 6.0, 987); 
  _collisionE   = manager->histogram("collision E", 60, 12.05, 12.15, 988); 
  _collisionECM = manager->histogram("collision ECM", 60, 10.55, 10.61, 989); 
  
  _totalPX  =manager->histogram("total X momentum", 60, -0.15, 0.05, 995); 
  _totalPY  =manager->histogram("total Y momentum", 60, -0.05, 0.05, 996); 
  _totalPZ  =manager->histogram("total Z momentum", 60,  5.8, 6.0, 997); 
  _totalE   =manager->histogram("total energy", 60, 12.05, 12.15, 998); 
  _totalECM =manager->histogram("total energy CM", 60, 10.55, 10.61, 999); 
  
  _BP  =manager->histogram("B+ Momentum", channels, lower,  0.6 , 1000 ); 
  _BM  =manager->histogram("B- Momentum", channels, lower,  0.6 , 1001 ); 
  _B0  =manager->histogram("B0 Momentum",  channels, lower,  0.6 , 1002 ); 
  _B0B  =manager->histogram("B0B Momentum", channels, lower,  0.6 , 1003 ); 
  
  _Btheta =manager->histogram("B0(B) cos(theta)", channels, 0.95, 1.0, 1004);
  _BthetaCM =manager->histogram("B0(B) cos(theta) CM", channels, -1.0, 1.0, 1005);
  _Bphi =manager->histogram("B0(B) phi/PI", channels, -1.0, 1.0, 1006);
  _BphiCM =manager->histogram("B0(B) phi/PI CM", channels, -1.0, 1.0, 1007);
  
  _PIP  =manager->histogram("pi+ not Ks/lambda Momentum", channels, lower,  
			    upper,1010  ); 
  _PIM  =manager->histogram("pi- not Ks/lambda Momentum", channels, lower,  
			    upper,1011  ); 
  _PI0  =manager->histogram("pi0 not Ks/lambda Momentum", channels, lower,  
			    upper,1012  ); 
  
  _PIPDST  =manager->histogram("pi+ from D*", channels, lower,  
			       upper,2013  ); 
  _PIMDST  =manager->histogram("pi- from D*", channels, lower,  
			       upper,2014  ); 
  _PI0DST  =manager->histogram("pi0 from D*", channels, lower,  
			       upper,2015  ); 
  
  _BxVtx =manager->histogram("B0(B) xVtx", channels, 
			     avgVtx.x() - 0.15, avgVtx.x() + 0.15, 1016 ); 
  _ByVtx =manager->histogram("B0(B) yVtx", channels, 
			     avgVtx.y() - 0.01, avgVtx.y() + 0.01, 1017 ); 
  _BzVtx =manager->histogram("B0(B) zVtx", channels, 
			     avgVtx.z() -4.0, avgVtx.z() + 4.0, 1018 ); 
  
  
  _KP  =manager->histogram("K+ Momentum", channels, lower,  upper, 1020 ); 
  _KM  =manager->histogram("K- Momentum", channels, lower,  upper, 1021 ); 
  _KS  =manager->histogram("KS Momentum", channels, lower,  upper, 1022  ); 
  _KL  =manager->histogram("KL Momentum", channels, lower,  upper, 1023  ); 
  
  _KPB0  =manager->histogram("K+ Momentum in B0 decay", channels, lower,  upper,1100  ); 
  _KMB0  =manager->histogram("K- Momentum in B0 decay", channels, lower,  upper,1101  ); 
  
  _KPB0B  =manager->histogram("K+ Momentum in B0B decay",channels, lower,  upper,1200  ); 
  _KMB0B  =manager->histogram("K- Momentum in B0B decay",channels, lower,  upper,1201  ); 
  
  _KPDIR  =manager->histogram("Prompt K+ momentum from B", channels, lower,  upper , 1030 ); 
  _KMDIR  =manager->histogram("Prompt K- momentum from B", channels, lower,  upper , 1031 ); 
  
  _KPDIRB0  =manager->histogram("Prompt K+ momentum in B0 decay", channels, lower,  upper , 1110 ); 
  _KMDIRB0  =manager->histogram("Prompt K- momentum in B0 decay",channels, lower,  upper , 1111 ); 
  
  _KPDIRB0B  =manager->histogram("Prompt K+ momentum in B0B decay",channels, lower,  upper , 1210 ); 
  _KMDIRB0B  =manager->histogram("Prompt K- momentum in B0B decay",channels, lower,  upper , 1211 ); 
  
  
  _DPB0  =manager->histogram("D+ Momentum in B0 decays",channels, lower,  upper,1340 ); 
  _DMB0  =manager->histogram("D- Momentum in B0 decays",channels, lower,  upper,1341 ); 
  _D0B0  =manager->histogram("D0 Momentum in B0 decays",channels, lower,  upper,1342 ); 
  _DBB0  =manager->histogram("DB Momentum in B0 decays",channels, lower,  upper,1343  ); 
  
  _DPB0B  =manager->histogram("D+ Momentum in B0B decays",channels, lower,  upper,1350); 
  _DMB0B  =manager->histogram("D- Momentum in B0B decays",channels, lower,  upper,1351); 
  _D0B0B  =manager->histogram("D0 Momentum in B0B decays",channels, lower,  upper,1352); 
  _DBB0B  =manager->histogram("DB Momentum in B0B decays",channels, lower,  upper,1353); 
  
  _DSTPB0  =manager->histogram("D*+ Momentum in B0 decays",channels, lower,  upper,1360); 
  _DSTMB0  =manager->histogram("D*- Momentum in B0 decays",channels, lower,  upper,1361); 
  _DST0B0  =manager->histogram("D*0 Momentum in B0 decays",channels, lower,  upper,1362); 
  _DSTBB0  =manager->histogram("D*B Momentum in B0 decays",channels, lower,  upper,1363); 
  
  _DSTPB0B  =manager->histogram("D*+ Momentum in B0B decays",channels, lower,  upper,1370); 
  _DSTMB0B  =manager->histogram("D*- Momentum in B0B decays",channels, lower,  upper,1371); 
  _DST0B0B  =manager->histogram("D*0 Momentum in B0B decays",channels, lower,  upper,1372); 
  _DSTBB0B  =manager->histogram("D*B Momentum in B0B decays",channels, lower,  upper,1373); 
  
  
  _DSP  =manager->histogram("DS+ Momentum",channelsd, lower,  upperd,1048); 
  _DSM  =manager->histogram("DS- Momentum",channelsd, lower,  upperd,1049); 
  
  _DSPB0  =manager->histogram("DS+ Momentum in B0 decays",channels, lower,  upper,1384); 
  _DSMB0  =manager->histogram("DS- Momentum in B0 decays",channels, lower,  upper,1349); 
  
  _DSPB0B  =manager->histogram("DS+ Momentum in B0B decays",channels, lower,  upper,1358); 
  _DSMB0B  =manager->histogram("DS- Momentum in B0B decays",channels, lower,  upper,1359); 
  
  _PSI  =manager->histogram("PSI Momentum",channels, lower,  upper,1050); 
  _PSIP  =manager->histogram("PSIP Momentum",channels, lower,  upper,1051); 
  
  
  _CHIC1  =manager->histogram("CHIC1 Momentum",channels, lower,  upper,1052); 
  _CHIC2  =manager->histogram("CHIC2 Momentum",channels, lower,  upper,1053); 
  _CHIC0  =manager->histogram("CHIC0 Momentum",channels, lower,  upper,1054); 
  _ETAC   =manager->histogram("ETAC Momentum",channels, lower,  upper,1058); 
  
  _ETA  =manager->histogram("ETA Momentum",channels, lower,  upper,1055); 
  
  
  _ETAP  =manager->histogram("ETAP Momentum",channels, lower,  upper,1056); 
  
  _PHI  =manager->histogram("PHI Momentum",channels, lower,  upper,1057); 
  
  
  _PP  =manager->histogram("P+ Momentum",channels, lower,  upper,1060); 
  _PM  =manager->histogram("P- Momentum",channels, lower,  upper,1061); 
  
  
  
  _LAM  =manager->histogram("LAM Momentum",channelsd, lower,  upperd,1062); 
  _ALAM  =manager->histogram("ALAM Momentum",channelsd, lower,  upperd,1063); 
  
  
  
  _LAMC  =manager->histogram("LAMC Momentum",channelsd, lower,  upperd ,1064); 
  _ALAMC  =manager->histogram("ALAMC Momentum",channelsd, lower,  upperd ,1065 ); 
  
  
  
  _EM  =manager->histogram("e- Momentum",channels, lower,  6. , 1070 ); 
  _EP  =manager->histogram("e+ Momentum",channels, lower,  6. , 1071 ); 
  
  _MUM  =manager->histogram("mu- Momentum",channels, lower,  upper , 1072 ); 
  _MUP  =manager->histogram("mu+ Momentum",channels, lower,  upper , 1073 ); 
  
  _MUX =manager->histogram("muon cos(x)",channels, -1., 1., 3730); 
  _MUY =manager->histogram("muon cos(y)",channels, -1., 1., 3731); 
  _MUZ =manager->histogram("muon cos(z)",channels, -1., 1., 3732); 
  
  _TAUM  =manager->histogram("tau- Momentum CM",channels, lower,  5.5 , 1074); 
  _TAUP  =manager->histogram("tau+ Momentum CM",channels, lower,  5.5 , 1075); 
  
  _TAUMLAB  =manager->histogram("tau- Momentum lab",channels, 0., 10. , 3740);
  _TAUPLAB  =manager->histogram("tau+ Momentum lab",channels, 0., 10. , 3750);
  
  _TAUMTHETA  =manager->histogram("tau- cos(theta)",channels, -1., 1., 1076); 
  _TAUPTHETA  =manager->histogram("tau+ cos(theta)",channels, -1., 1., 1077); 
  
  _TAUMTHETACM  =manager->histogram("tau- cos(theta) CM",channels, -1., 1., 3076); 
  _TAUPTHETACM  =manager->histogram("tau+ cos(theta) CM",channels, -1., 1., 3077); 
  
  _TAUMPHI  =manager->histogram("tau- phi/PI",channels, -1., 1., 3078); 
  _TAUPPHI  =manager->histogram("tau+ phi/PI",channels, -1., 1., 3079); 
  
  _TAULHELIC  =manager->histogram("tau->l cos helicity",channels, -1., 1., 3781); 
  _TAUPIHELIC  =manager->histogram("tau->pi cos helicity",channels, -1., 1., 3782); 
  _TAUKHELIC  =manager->histogram("tau->K cos helicity",channels, -1., 1., 3783); 
  
  _EMDIR  =manager->histogram("Prompt e- momentum in B decays",channels, lower,  upper , 1080 ); 
  _EPDIR  =manager->histogram("Prompt e+ momentum in B decays",channels, lower,  upper , 1081 ); 
  
  
  _MUMDIR  =manager->histogram("Prompt mu- momentum in B decays",channels, lower,  upper , 1082 ); 
  _MUPDIR  =manager->histogram("Prompt mu+ momentum in B decays",channels, lower,  upper , 1083 ); 
  
  _TAUMDIR  =manager->histogram("Prompt tau- momentum in B decays",channels, lower,  upper , 1084 ); 
  _TAUPDIR  =manager->histogram("Prompt tau+ momentum in B decays",channels, lower,  upper , 1085 ); 
  
  _EMDDIR  =manager->histogram("Prompt e- momentum in D decays",channels, lower,  upper , 1086 ); 
  _EPDDIR  =manager->histogram("Prompt e+ momentum in D decays",channels, lower,  upper , 1087 ); 
  
  _MUMDDIR  =manager->histogram("Prompt mu- momentum in D decays",channels, lower,  upper , 1088 ); 
  _MUPDDIR  =manager->histogram("Prompt mu+ momentum in D decays",channels, lower,  upper , 1089 ); 
  
  _TAUMDDIR  =manager->histogram("Prompt tau- momentum in D decays",channels, lower,  upper , 1090 ); 
  _TAUPDDIR  =manager->histogram("Prompt tau+ momentum in D decays",channels, lower,  upper , 1091 ); 
  
  _EMDSDIR  =manager->histogram("Prompt e- momentum in Ds decays",channels, lower,  upper , 1092 ); 
  _EPDSDIR  =manager->histogram("Prompt e+ momentum in Ds decays",channels, lower,  upper , 1093 ); 
  
  _MUMDSDIR  =manager->histogram("Prompt mu- momentum in Ds decays",channels, lower,  upper , 1094 ); 
  _MUPDSDIR  =manager->histogram("Prompt mu+ momentum in Ds decays",channels, lower,  upper , 1095 ); 
  
  _TAUMDSDIR  =manager->histogram("Prompt tau- momentum in Ds decays",channels, lower,  upper , 1096 ); 
  _TAUPDSDIR  =manager->histogram("Prompt tau+ momentum in Ds decays",channels, lower,  upper , 1097 ); 
  
  _EMTDIR  =manager->histogram("Prompt e- momentum in tau decays",channels, 0., 6. , 1066 ); 
  _EPTDIR  =manager->histogram("Prompt e+ momentum in tau decays",channels, 0., 6. , 1067 ); 
  
  _MUMTDIR  =manager->histogram("Prompt mu- momentum in tau decays",channels, 0., 6. , 1068 ); 
  _MUPTDIR  =manager->histogram("Prompt mu+ momentum in tau decays",channels, 0., 6. , 1069 ); 
  
  _EMB0  =manager->histogram("e- Momentum in B0 decays",channels, lower,  upper , 1120 ); 
  _EPB0  =manager->histogram("e+ Momentum in B0 decays",channels, lower,  upper , 1121 ); 
  
  _MUMB0  =manager->histogram("mu- Momentum in B0 decays",channels, lower,  upper , 1122 ); 
  _MUPB0  =manager->histogram("mu+ Momentum in B0 decays",channels, lower,  upper , 1123 ); 
  
  
  _TAUMB0  =manager->histogram("tau- Momentum in B0 decays",channels, lower,  upper , 1124 ); 
  _TAUPB0  =manager->histogram("tau+ Momentum in B0 decays",channels, lower,  upper , 1125 ); 
  
  
  _EMDIRB0  =manager->histogram("Prompt e- Momentum in B0 decays",channels, lower,  upper , 1130 ); 
  _EPDIRB0  =manager->histogram("Prompt e+ Momentum in B0 decays",channels, lower,  upper , 1131 ); 
  
  _MUMDIRB0  =manager->histogram("Prompt mu- Momentum in B0 decays",channels, lower,  upper , 1132 ); 
  _MUPDIRB0  =manager->histogram("Prompt mu+ Momentum in B0 decays",channels, lower,  upper , 1133 ); 
  
  _TAUMDIRB0  =manager->histogram("Prompt tau- Momentum in B0 decays",channels, lower,  upper , 1134 ); 
  _TAUPDIRB0  =manager->histogram("Prompt tau+ Momentum in B0 decays",channels, lower,  upper , 1135 ); 
  
  
  
  _EMDDIRB0  =manager->histogram("Prompt e- momentum from D in B0 decays",channels, lower,  upper , 1140 ); 
  _EPDDIRB0  =manager->histogram("Prompt e+ Momentum from D in B0 decays",channels, lower,  upper , 1141 ); 
  
  _MUMDDIRB0  =manager->histogram("Prompt mu- Momentum from D in B0 decays",channels, lower,  upper , 1142 ); 
  _MUPDDIRB0  =manager->histogram("Prompt mu+ Momentum from D in B0 decays",channels, lower,  upper , 1143 ); 
  
  _TAUMDDIRB0  =manager->histogram("Prompt tau- Momentum from D in B0 decays",channels, lower,  upper , 1144 ); 
  _TAUPDDIRB0  =manager->histogram("Prompt tau+ Momentum from D in B0 decays",channels, lower,  upper , 1145 ); 
  
  _EMDSDIRB0  =manager->histogram("Prompt e- Momentum from Ds in B0 decays",channels, lower,  upper , 1150 ); 
  _EPDSDIRB0  =manager->histogram("Prompt e+ Momentum from Ds in B0 decays",channels, lower,  upper , 1151 ); 
  
  _MUMDSDIRB0  =manager->histogram("Prompt mu- Momentum from Ds in B0 decays",channels, lower,  upper , 1152 ); 
  _MUPDSDIRB0  =manager->histogram("Prompt mu+ Momentum from Ds in B0 decays",channels, lower,  upper , 1153 ); 
  
  _TAUMDSDIRB0  =manager->histogram("Prompt tau- Momentum from Ds in B0 decays",channels, lower,  upper , 1154 ); 
  _TAUPDSDIRB0  =manager->histogram("Prompt tau+ Momentum frm Ds in B0 decays",channels, lower,  upper , 1155 ); 
  
  
  _EMTDIRB0  =manager->histogram("Prompt e- Momentum from tau in B0 decays",channels, lower,  upper , 1160 ); 
  _EPTDIRB0  =manager->histogram("Prompt e+ Momentum from tau in B0 decays",channels, lower,  upper , 1161 ); 
  
  
  _MUMTDIRB0  =manager->histogram("Prompt mu- Momentum from tau in B0 decays",channels, lower,  upper , 1162 ); 
  _MUPTDIRB0  =manager->histogram("Prompt mu+ Momentum from tau in B0 decays",channels, lower,  upper , 1163 ); 
  
  
  _EMB0B  =manager->histogram("e- Momentum in B0B decays",channels, lower,  upper , 1220 ); 
  _EPB0B  =manager->histogram("e+ Momentum in B0B decays",channels, lower,  upper , 1221 ); 
  
  
  _MUMB0B  =manager->histogram("mu- Momentum in B0B decays",channels, lower,  upper , 1222 ); 
  _MUPB0B  =manager->histogram("mu+ Momentum in B0B decays",channels, lower,  upper , 1223 ); 
  
  
  _TAUMB0B  =manager->histogram("tau- Momentum in B0B decays",channels, lower,  upper , 1224 ); 
  _TAUPB0B  =manager->histogram("tau+ Momentum in B0B decays",channels, lower,  upper , 1225 ); 
  
  
  _EMDIRB0B  =manager->histogram("Prompt e- Momentum in B0B decays",channels, lower,  upper , 1230 ); 
  _EPDIRB0B  =manager->histogram("Prompt e+ Momentum in B0B decays",channels, lower,  upper , 1231 ); 
  
  
  _MUMDIRB0B  =manager->histogram("Prompt mu- Momentum in B0B decays",channels, lower,  upper , 1232 ); 
  _MUPDIRB0B  =manager->histogram("Prompt mu+ Momentum in B0B decays",channels, lower,  upper , 1233 ); 
  
  
  _TAUMDIRB0B  =manager->histogram("Prompt tau- Momentum in B0B decays",channels, lower,  upper , 1234 ); 
  _TAUPDIRB0B  =manager->histogram("Prompt tau+ Momentum in B0B decays",channels, lower,  upper , 1235 ); 
  
  
  _EMDDIRB0B  =manager->histogram("Prompt e- Momentum from D in B0B decays",channels, lower,  upper , 1240 ); 
  _EPDDIRB0B  =manager->histogram("Prompt e+ Momentum from D in B0B decays",channels, lower,  upper , 1241 ); 
  
  
  _MUMDDIRB0B  =manager->histogram("Prompt mu- Momentum from D in B0B decays",channels, lower,  upper , 1242 ); 
  _MUPDDIRB0B  =manager->histogram("Prompt mu+ Momentum from D in B0B decays",channels, lower,  upper , 1243 ); 
  
  
  _TAUMDDIRB0B  =manager->histogram("Prompt tau- Momentum from D in B0B decays",channels, lower,  upper , 1244 ); 
  _TAUPDDIRB0B  =manager->histogram("Prompt tau+ Momentum from D in B0B decays",channels, lower,  upper , 1245 ); 
  
  
  _EMDSDIRB0B  =manager->histogram("Prompt e- Momentum from Ds in B0B decays",channels, lower,  upper , 1250 ); 
  _EPDSDIRB0B  =manager->histogram("Prompt e+ Momentum from Ds in B0B decays",channels, lower,  upper , 1251 ); 
  
  
  _MUMDSDIRB0B  =manager->histogram("Prompt mu- Momentum from Ds in B0B decays",channels, lower,  upper , 1252 ); 
  _MUPDSDIRB0B  =manager->histogram("Prompt mu+ Momentum from Ds in B0B decays",channels, lower,  upper , 1253 ); 
  
  _TAUMDSDIRB0B  =manager->histogram("Prompt tau- Momentum from Ds in B0B decays",channels, lower,  upper , 1254 ); 
  _TAUPDSDIRB0B  =manager->histogram("Prompt tau+ Momentum from Ds in B0B decays",channels, lower,  upper , 1255 ); 
  
  
  _EMTDIRB0B  =manager->histogram("Prompt e- Momentum from tau in B0B decays",channels, lower,  upper , 1260 ); 
  _EPTDIRB0B  =manager->histogram("Prompt e+ Momentum from tau in B0B decays",channels, lower,  upper , 1261 ); 
  
  _MUMTDIRB0B  =manager->histogram("Prompt mu- Momentum from tau in B0B decays",channels, lower,  upper , 1262 ); 
  _MUPTDIRB0B  =manager->histogram("Prompt mu+ Momentum from tau in B0B decays",channels, lower,  upper , 1263 ); 
  
  _CISRgamma  =manager->histogram("ISR gamma Momentum",channelsd, lower,  upperd,2050  ); 
  
  
  // those are for charm production
  
  
  
  _CDCOS  =manager->histogram("D cos",channelsd, -1.0,  1.0,2039  ); 
  
  _CDP  =manager->histogram("D+ Momentum",channelsd, lower,  upperd,2040  ); 
  _CDM  =manager->histogram("D- Momentum",channelsd, lower,  upperd,2041  ); 
  _CD0  =manager->histogram("D0 Momentum",channelsd, lower,  upperd,2042  ); 
  _CDB  =manager->histogram("DB Momentum",channelsd, lower,  upperd,2043  ); 
  
  _CDSTP  =manager->histogram("D*+ Momentum",channelsd, lower,  upperd,2044); 
  _CDSTM  =manager->histogram("D*- Momentum",channelsd, lower,  upperd,2045); 
  _CDST0  =manager->histogram("D*0 Momentum",channelsd, lower,  upperd,2046); 
  _CDSTB  =manager->histogram("D*B Momentum",channelsd, lower,  upperd,2047); 
  
  _CDD1p  =manager->histogram("DDstar+ 2420 1+",channelsd, lower,  upperd,2100  ); 
  _CDD1m  =manager->histogram("DDstar- 2420 1+",channelsd, lower,  upperd,2101  ); 
  _CDD2p  =manager->histogram("DDstar+ 2460 2+",channelsd, lower,  upperd,2102  ); 
  _CDD2m  =manager->histogram("DDstar- 2460 2+",channelsd, lower,  upperd,2103  ); 
  _CDD3p  =manager->histogram("DDstar+ 2420 1-",channelsd, lower,  upperd,2104  ); 
  _CDD3m  =manager->histogram("DDstar- 2420 1-",channelsd, lower,  upperd,2105  ); 
  _CDD4p  =manager->histogram("DDstar+ 2400 0-",channelsd, lower,  upperd,2106  ); 
  _CDD4m  =manager->histogram("DDstar- 2400 0-",channelsd, lower,  upperd,2107  ); 
  
  
  _CDD1  =manager->histogram("DDstar0 2420 1+",channelsd, lower,  upperd,2108 ); 
  _CDD1b  =manager->histogram("DDstar0b 2420 1+",channelsd, lower,  upperd,2109  ); 
  _CDD2  =manager->histogram("DDstar0 2460 2+",channelsd, lower,  upperd,2110  ); 
  _CDD2b  =manager->histogram("DDstar0b 2460 2+",channelsd, lower,  upperd,2111  ); 
  _CDD3  =manager->histogram("DDstar0 2420 1-",channelsd, lower,  upperd,2112  ); 
  _CDD3b  =manager->histogram("DDstar0b 2420 1-",channelsd, lower,  upperd,2113  ); 
  _CDD4  =manager->histogram("DDstar0 2400 0-",channelsd, lower,  upperd,2114  ); 
  _CDD4b  =manager->histogram("DDstar0b 2400 0-",channelsd, lower,  upperd,2115  ); 
  
  _CDPDstar1 =manager->histogram(" D+ from D*+ ",channelsd, lower,  upperd,2120  );
  _CDPDstar2 =manager->histogram(" D+ from D2*+ ",channelsd, lower,  upperd,2121  );
  _CDPDstar3 =manager->histogram(" D+ from D1+ ",channelsd, lower,  upperd,2122 );
  _CDPDstar4 =manager->histogram(" D+ from D'1+ ",channelsd, lower,  upperd,2123  );
  _CDPDstar5 =manager->histogram(" D+ from D0*+ ",channelsd, lower,  upperd,2124  );
  _CDPDstar6 =manager->histogram(" D+ from D0*0 ",channelsd, lower,  upperd,2125  );
  _CDPDstar7 =manager->histogram(" D+ from D2*0 ",channelsd, lower,  upperd,2126  );
  _CDPDstar8 =manager->histogram(" D+ from D10 ",channelsd, lower,  upperd,2127  );
  _CDPDstar9 =manager->histogram(" D+ from D'10 ",channelsd, lower,  upperd,2128  );
  
  _CDMDstar1 =manager->histogram(" D- from D*- ",channelsd, lower,  upperd,2130  );
  _CDMDstar2 =manager->histogram(" D- from D2*+ ",channelsd, lower,  upperd,2131  );
  _CDMDstar3 =manager->histogram(" D- from D1+ ",channelsd, lower,  upperd,2132 );
  _CDMDstar4 =manager->histogram(" D- from D'1+ ",channelsd, lower,  upperd,2133  );
  _CDMDstar5 =manager->histogram(" D- from D0*+ ",channelsd, lower,  upperd,2134  );
  _CDMDstar6 =manager->histogram(" D- from D0*0 ",channelsd, lower,  upperd,2135  );
  _CDMDstar7 =manager->histogram(" D- from D2*0 ",channelsd, lower,  upperd,2136  );
  _CDMDstar8 =manager->histogram(" D- from D10 ",channelsd, lower,  upperd,2137  );
  _CDMDstar9 =manager->histogram(" D- from D'10 ",channelsd, lower,  upperd,2138  );
  
  _CD0Dstar0 =manager->histogram(" D0 from D*+ ",channelsd, lower,  upperd,2140  );
  _CD0Dstar1 =manager->histogram(" D0 from D*0 ",channelsd, lower,  upperd,2141  );
  _CD0Dstar2 =manager->histogram(" D0 from D2*+ ",channelsd, lower,  upperd,2142  );
  _CD0Dstar3 =manager->histogram(" D0 from D1+ ",channelsd, lower,  upperd,2143  );
  _CD0Dstar4 =manager->histogram(" D0 from D'1+ ",channelsd, lower,  upperd,2144  );
  _CD0Dstar5 =manager->histogram(" D0 from D0*+ ",channelsd, lower,  upperd,2145  );
  _CD0Dstar6 =manager->histogram(" D0 from D0*0 ",channelsd, lower,  upperd,2146  );
  _CD0Dstar7 =manager->histogram(" D0 from D2*0 ",channelsd, lower,  upperd,2147  );
  _CD0Dstar8 =manager->histogram(" D0 from D10 ",channelsd, lower,  upperd,2148  );
  _CD0Dstar9 =manager->histogram(" D0 from D'10 ",channelsd, lower,  upperd,2149  );
  
  _CDBDstar0 =manager->histogram(" D0b from D*- ",channelsd, lower,  upperd,2180  );
  _CDBDstar1 =manager->histogram(" D0b from D*0b ",channelsd, lower,  upperd,2181  );
  _CDBDstar2 =manager->histogram(" D0b from D2*+ ",channelsd, lower,  upperd,2182  );
  _CDBDstar3 =manager->histogram(" D0b from D1+ ",channelsd, lower,  upperd,2183  );
  _CDBDstar4 =manager->histogram(" D0b from D'1+ ",channelsd, lower,  upperd,2184  );
  _CDBDstar5 =manager->histogram(" D0b from D0*+ ",channelsd, lower,  upperd,2185  );
  _CDBDstar6 =manager->histogram(" D0b from D0*0 ",channelsd, lower,  upperd,2186  );
  _CDBDstar7 =manager->histogram(" D0b from D2*0 ",channelsd, lower,  upperd,2187  );
  _CDBDstar8 =manager->histogram(" D0b from D10 ",channelsd, lower,  upperd,2188 );
  _CDBDstar9 =manager->histogram(" D0b from D'10 ",channelsd, lower,  upperd,2189  );
  
  _CDSDstarP2 =manager->histogram(" D*+ from D2*+ ",channelsd, lower,  upperd,2190  );
  _CDSDstarP3 =manager->histogram(" D*+ from D1+ ",channelsd, lower,  upperd,2191  );
  _CDSDstarP4 =manager->histogram(" D*+ from D'1+ ",channelsd, lower,  upperd,2192  );
  _CDSDstarP5 =manager->histogram(" D*+ from D0*+ ",channelsd, lower,  upperd,2193  );
  _CDSDstarP6 =manager->histogram(" D*+ from D0*0 ",channelsd, lower,  upperd,2194  );
  _CDSDstarP7 =manager->histogram(" D*+ from D2*0 ",channelsd, lower,  upperd,2195  );
  _CDSDstarP8 =manager->histogram(" D*+ from D10 ",channelsd, lower,  upperd,2196  );
  _CDSDstarP9 =manager->histogram(" D*+ from D'10 ",channelsd, lower,  upperd,2197  );
  
  _CDSDstarM2 =manager->histogram(" D*- from D2*+ ",channelsd, lower,  upperd,2150  );
  _CDSDstarM3 =manager->histogram(" D*- from D1+ ",channelsd, lower,  upperd,2151  );
  _CDSDstarM4 =manager->histogram(" D*- from D'1+ ",channelsd, lower,  upperd,2152  );
  _CDSDstarM5 =manager->histogram(" D*- from D0*+ ",channelsd, lower,  upperd,2153  );
  _CDSDstarM6 =manager->histogram(" D*- from D0*0 ",channelsd, lower,  upperd,2154  );
  _CDSDstarM7 =manager->histogram(" D*- from D2*0 ",channelsd, lower,  upperd,2155  );
  _CDSDstarM8 =manager->histogram(" D*- from D10 ",channelsd, lower,  upperd,2156  );
  _CDSDstarM9 =manager->histogram(" D*- from D'10 ",channelsd, lower,  upperd,2157  );
  
  _CDSDstarO2 =manager->histogram(" D*0 from D2*+ ",channelsd, lower,  upperd,2160  );
  _CDSDstarO3 =manager->histogram(" D*0 from D1+ ",channelsd, lower,  upperd,2161  );
  _CDSDstarO4 =manager->histogram(" D*0 from D'1+ ",channelsd, lower,  upperd,2162  );
  _CDSDstarO5 =manager->histogram(" D*0 from D0*+ ",channelsd, lower,  upperd,2163  );
  _CDSDstarO6 =manager->histogram(" D*0 from D0*0 ",channelsd, lower,  upperd,2164  );
  _CDSDstarO7 =manager->histogram(" D*0 from D2*0 ",channelsd, lower,  upperd,2165  );
  _CDSDstarO8 =manager->histogram(" D*0 from D10 ",channelsd, lower,  upperd,2166  );
  _CDSDstarO9 =manager->histogram(" D*0 from D'10 ",channelsd, lower,  upperd,2167  );
  
  _CDSDstarB2 =manager->histogram(" D*0b from D2*+ ",channelsd, lower,  upperd,2170  );
  _CDSDstarB3 =manager->histogram(" D*0b from D1+ ",channelsd, lower,  upperd,2171  );
  _CDSDstarB4 =manager->histogram(" D*0b from D'1+ ",channelsd, lower,  upperd,2172  );
  _CDSDstarB5 =manager->histogram(" D*0b from D0*+ ",channelsd, lower,  upperd,2173  );
  _CDSDstarB6 =manager->histogram(" D*0b from D0*0 ",channelsd, lower,  upperd,2174  );
  _CDSDstarB7 =manager->histogram(" D*0b from D2*0 ",channelsd, lower,  upperd,2175  );
  _CDSDstarB8 =manager->histogram(" D*0b from D10 ",channelsd, lower,  upperd,2176  );
  _CDSDstarB9 =manager->histogram(" D*0b from D'10 ",channelsd, lower,  upperd,2177  );
  
  _collisionTuple = manager->ntuple("collision Tuple");
  
  return AppResult::OK;
}

AppResult
GqaMCAnalysis::event(AbsEvent* anEvent){
  ++_evtNum;   // count events
  
  if ((_evtNum%1000)==0) ErrMsg(routine) << "Processing event:"
					 << _evtNum << endmsg;
  
  
  static const double PI = Constants::pi;
  
  HepAList< EventInfo >* infoList=NULL;
  getTmpAList(anEvent, infoList, _eventInfoList.value());
  EventInfo* eventInfo = infoList->first();
  
  // Plot the parameters of the collision:
  const PepCollision * collision = 
    Ifd<PepCollision>::get(anEvent, _pepCollisionKey.value());
  
  HepLorentzVector collision4Mom;
  
  if (0 != collision){
    if (_verbose.value()){
      ErrMsg(routine) << *collision << endmsg;
    }
    
    HepPoint collisionVertex = collision->vertex();
    _collisionX->accumulate(collisionVertex.x());
    _collisionY->accumulate(collisionVertex.y());
    _collisionZ->accumulate(collisionVertex.z());
    
    collision4Mom = collision->total4Momentum();
    _collisionPX->accumulate(collision4Mom.x());
    _collisionPY->accumulate(collision4Mom.y());
    _collisionPZ->accumulate(collision4Mom.z());
    _collisionE->accumulate(collision4Mom.t());
    
    _collisionECM->accumulate(collision->energyCM());
    
    if (_makeNtuple.value()){
      _collisionTuple->column("x", collisionVertex.x());
      _collisionTuple->column("y", collisionVertex.y());
      _collisionTuple->column("z", collisionVertex.z());
      _collisionTuple->column("px", collision4Mom.x());
      _collisionTuple->column("py", collision4Mom.y());
      _collisionTuple->column("pz", collision4Mom.z());
      _collisionTuple->column("e", collision4Mom.t());
      _collisionTuple->column("ecm", collision->energyCM());
      _collisionTuple->dumpData();
    }
  }
  
  HepAList<BtaCandidate>* mcList;
  getTmpAList (anEvent, mcList, _btaTruthList.value());
  
  HepAListIterator<BtaCandidate> iter(*(mcList));
  
  if (true == _verbose.value()){
    ErrMsg(routine)
      << name() << " # of tracks read = " << mcList->length() << endmsg;
  }
  
  // Plot the total 4-momentum of observed partcles:
  BtaCandidate* MCptr;  

  HepLorentzVector p40, p4p, p4m;

  bool pipF = false;
  bool pimF = false;
  bool pi0F = false;
  while (MCptr = iter()) {
    BtaCandidate * motherPart = MCptr->theMother();
    if (0 != motherPart && motherPart->nDaughters() == 3) {
      PdtLund::LundType mother = motherPart->pdtEntry()->lundId();
      if (_DFlavor == mother) {
	if (PdtLund::pi_plus == MCptr->pdtEntry()->lundId()) {
	  p4p = MCptr->p4();
	  if (false == pipF) {
	    pipF = true;
	  }
	}
	else if (PdtLund::pi_minus == MCptr->pdtEntry()->lundId()) {
	  p4m = MCptr->p4();
	  if (false == pimF) {
	    pimF = true;
	  }
	}
	else if (PdtLund::pi0 == MCptr->pdtEntry()->lundId()) {
	  p40 = MCptr->p4();
	  if (false == pi0F) {
	    pi0F = true;
	  }
	}
      }
    }
  }

  double M2p0 = (p4p + p40).mag();
  double M2m0 = (p4m + p40).mag();
  double M2pm = (p4p + p4m).mag();

  if (pipF && pimF && pi0F) {
    _tuple->column("d0pppupmass", M2p0);
    _tuple->column("d0ppmupmass", M2m0);
    _tuple->column("d0ppupmass" , M2pm);
    _tuple->column("M12", M2p0*M2p0);
    _tuple->column("M13", M2m0*M2m0);
    _tuple->column("M23" , M2pm*M2pm);
    _tuple->dumpData();
  }


  return AppResult::OK;
}

//--------------------------------------------------------------------
AppResult
GqaMCAnalysis::endJob( AbsEvent* anEvent )
{
  ErrMsg(routine) << name( ) << " endJob" << endmsg;
  return AppResult::OK;
}
 
//--------------------------------------------------------------------
bool 
GqaMCAnalysis::printCand(const char * particle) const {
  if (_printCands.value().length() > 0){
    istrstream candsStream(_printCands.value().c_str());
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


