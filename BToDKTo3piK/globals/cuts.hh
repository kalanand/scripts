// $Id: cuts.hh,v 1.19 2006/06/15 18:21:45 abi Exp $
// Define cuts
// Cuts should be defined using the TCut constructor given a name AND a title
// The name can be used in scripts to label graphs 
// Note: If you do cut3 = cut2 && cut1, cut3 will by default inherit the name from cut2

#ifndef CUTS_HH
#define CUTS_HH

#include "TCut.h"

TCut cutPiVT("very tight pi PID", "Pippidbit==15&&Pimpidbit==15");
TCut cutPiT("tight pi PID", "Pippidbit>6&&Pimpidbit>6");
TCut cutNNq("cutNNq","nnout>0.1");
TCut cutNNd("cutNNd","bknnout>0.25");
TCut cutDtoKpi("cutDtoKpi", "(mixneumass<1.84||mixneumass>1.89)");
TCut cutNoKs("noKs", "d0truekspp<=0");
TCut cutKsVeto("KsVeto", "((0.489>d0ppupmass)||(d0ppupmass>0.508))&&(d0flightdist<1.5)");

TCut cutGoodD("GoodD","(Trued0flg>0&&d0recdec==3)");
TCut cutBadD("BadD","(Trued0flg<=0||d0recdec!=3)");

TCut cutDKBadD("DKBadD","(B1decmode==16 || B2decmode==16) && excltruth!=16");
TCut cutDKGoodD("DKGoodD","excltruth==16");

TCut cutDPiBadD("DPiBadD",TCut("(B1decmode==156 || B2decmode==156)") && cutBadD);
TCut cutDPiGoodD("DPiGoodD",TCut("(B1decmode==156 || B2decmode==156)") && cutGoodD);

TCut cutDPiX("DPiX","B1decmode==50 || (B1decmode>79 && B1decmode<109) || (B1decmode>149 && B1decmode<156) || (B1decmode>156 && B1decmode<199) || B2decmode==50 || (B2decmode>79 && B2decmode<109) || (B2decmode>149 && B2decmode<156) || (B2decmode>156 && B2decmode<199)");

TCut cutDKX("DKX","(B1decmode>9 && B1decmode<16) || (B1decmode>16 && B1decmode<19) || (B1decmode>59 && B1decmode<79) || (B1decmode>109 && B1decmode<149) || (B2decmode>9 && B2decmode<16) || (B2decmode>16 && B2decmode<19) || (B2decmode>59 && B2decmode<79) || (B2decmode>109 && B2decmode<149)");

TCut cutBBBadD("BBBadD",TCut("(B1decmode>199 || B2decmode>199) || (B1decmode<=0 && B2decmode<=0)")&&cutBadD);

TCut cutBBGoodD("BBGoodD",TCut("B1decmode<=0 && B2decmode<=0")&&cutGoodD);

TCut cutqqBadD("qqBadD",TCut("B1decmode<=0 && B2decmode<=0") && cutBadD);
TCut cutqqGoodD("qqGoodD",TCut("B1decmode<=0 && B2decmode<=0") && cutGoodD);

// cut on runnumber
TCut cutRun14("run14","runnumber<52000");
TCut cutRun5("run5","runnumber>=52000");


// signal region cuts
TCut cutDeltaE("#DeltaE signal cut","-0.07<Deltae&&Deltae<0.06");
TCut cutmES("m_{ES} signal cut","mes>5.272");
TCut cutMD("m_{D} signal cut","1.83<d0mass&&d0mass<1.895");

// Cut and weight for PID efficiencies
// Use this with the "*" operator, i.e. TCut cut = pidEff*cutSigReg;
TCut pidEff("pid eff","(Pippeff>=0 && Pimpeff>=0 && Hdtrkpeff>=0)*((1-Pippeff)*(1-Pimpeff)*Hdtrkpeff)");

// Charge cuts:
TCut cutPlus("B+", "Hdtrkchge>0");
TCut cutMinus("B-", "Hdtrkchge<0");


// Composite cuts. Initialized in setupCuts()
TCut cutNN;
TCut cutBasic;
TCut cutSigReg;
TCut readCut;

// sideband cuts
TCut cutSBUpperDE;
TCut cutSBLowerDE;
TCut cutSBmES;
TCut cutSBUpperMD;
TCut cutSBLowerMD;
TCut cutSBMD;


#endif
