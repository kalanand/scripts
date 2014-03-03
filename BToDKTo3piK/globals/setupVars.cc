// $Id: setupVars.cc,v 1.35 2006/08/04 22:27:24 fwinkl Exp $
// Sets up the global variables

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooCategory.hh"
#include "RooFitCore/RooFormulaVar.hh"
#include "RooFitCore/RooArgSet.hh"

#include "BToDKTo3piK/BdkPdfAbsBase.hh"

#include "../BToDKTo3piK/globals/vars.hh"

// forward
void setBinning(RooRealVar *r, Double_t gev);

void setupVars() {
  cout << "--- START setup() ---" << endl;

  // initialize the variables:
  randAdd = new RooRealVar("randAdd", "randAdd", 0, 1); 
  mes = new RooRealVar("mes","m_{ES}", 5.272, 5.30,"GeV/c^{2}");
  Deltae = new RooRealVar("Deltae", "#Delta E", -0.07, 0.06, "GeV");
  d0mass = new RooRealVar("d0mass", "m_{D}",1.83,1.895, "GeV/c^{2}");
  nnout = new RooRealVar("nnout", "q", 0.1, 1.0);
  bknnout = new RooRealVar("bknnout", "d", 0.25, 1.0);
  Hdtrkkaonnn = new RooRealVar("Hdtrkkaonnn", "trknn", 0.0, 1.0); 
  excltruth = new RooRealVar("excltruth", "excltruth", -100, 200);
  Dalitz1 = new RooRealVar("Dalitz1", "Dalitz variable 1", 0.0, 1.0); 
  Dalitz2 = new RooRealVar("Dalitz2", "Dalitz variable 2", 0.0, 1.0); 
  Pippidbit = new RooRealVar("Pippidbit", "pi+ PID bit", 0 , 16);
  Pimpidbit = new RooRealVar("Pimpidbit", "p- PID bit", 0 , 16);
  d0ppupmass = new RooRealVar("d0ppupmass", "two pion trk mass", 0, 2.9);
  Trued0flg = new RooRealVar("Trued0flg", "MC truth for D0", -200, 100);
  d0recdec = new RooRealVar("d0recdec", "d0 decay bit", -2, 10);
  d0flightdist = new RooRealVar("d0flightdist", "d0 flight distance", 0, 1000);

  // PID data efficencies
  Pippeff = new RooRealVar("Pippeff","pi+ PID eff",0,1);
  Pimpeff = new RooRealVar("Pimpeff","pi+ PID eff",0,1);
  Hdtrkpeff = new RooRealVar("Hdtrkpeff","pi+ PID eff",0,1);
  eventPidEff = new RooFormulaVar("eventPidEff","event PID eff",
                                  "(1-@0)*(1-@1)*@2",
                                  RooArgList(*Pippeff,*Pimpeff,*Hdtrkpeff));
  
  // Unsquared Dalitz masses
  mass12 = new RooRealVar("d0pppupmass", "M(#pi^{+}#pi^{0})", 0, 2, "GeV/c^{2}");
  mass13 = new RooRealVar("d0ppmupmass", "M(#pi^{-}#pi^{0})", 0, 2, "GeV/c^{2}");
  // RooFormularVar's to calculate squared Dalitz masses
  s12 = new RooFormulaVar("m12","M(#pi^{+}#pi^{0})^{2}","@0*@0",*mass12);
  s13 = new RooFormulaVar("m13","M(#pi^{-}#pi^{0})^{2}","@0*@0",*mass13);
  s12->setUnit("GeV^{2}/c^{4}");
  s13->setUnit(s12->getUnit());
  // Can be used to access squared Dalitz masses in a RooDataset
  m12 = new RooRealVar("m12","M(#pi^{+}#pi^{0})^{2}", 0, 3, "GeV^{2}/c^{4}");
  m13 = new RooRealVar("m13","M(#pi^{-}#pi^{0})^{2}", 0, 3, "GeV^{2}/c^{4}");
  m23 = new RooRealVar("m23","M(#pi^{+}#pi^{-})^{2}", 0, 3, "GeV^{2}/c^{4}");

  // Same for true Dalitz variables:
  mass12mc = new RooRealVar("d0pppmcmass", "M(#pi^{+}#pi^{0})", 0, 2, "GeV/c^{2}");
  mass13mc = new RooRealVar("d0ppmmcmass", "M(#pi^{-}#pi^{0})", 0, 2, "GeV/c^{2}");

  s12mc = new RooFormulaVar("m12mc","M(#pi^{+}#pi^{0})^{2}","@0*@0",*mass12);
  s13mc = new RooFormulaVar("m13mc","M(#pi^{-}#pi^{0})^{2}","@0*@0",*mass13);
  s12mc->setUnit("GeV^{2}/c^{4}");
  s13mc->setUnit(s12->getUnit());

  m12mc = new RooRealVar("m12mc","M(#pi^{+}#pi^{0})^{2}", 0, 3, "GeV^{2}/c^{4}");
  m13mc = new RooRealVar("m13mc","M(#pi^{-}#pi^{0})^{2}", 0, 3, "GeV^{2}/c^{4}");

  // nnout transformation
  qprimeF = new RooFormulaVar("qprime","","atanh((@0-.55)/0.45)",RooArgList(*nnout));
  qprime = new RooRealVar("qprime","q'",-20,20);
  qprime->setRange("plot",-5,5);    // range for making plots

  dprimeF = new RooFormulaVar("dprime","","atanh((@0-.625)/0.375)",RooArgList(*bknnout));
  dprime = new RooRealVar("dprime","d'",-20,20);
  dprime->setRange("plot",-5,5);    // range for making plots

  B1decmode = new RooRealVar("B1decmode","1st B decmode", -200, 400); 
  B2decmode = new RooRealVar("B2decmode","2nd B decmode", -200, 400); 

  R2 = new RooRealVar("R2","R2", 0.0, 0.5); 

  mixneumass = new RooRealVar("mixneumass", "mass of Kpi pair", 0, 5);
 
  Hdtrkchge = new RooCategory("Hdtrkchge", "K charge"); 
  Hdtrkchge->defineType("+",  1);
  Hdtrkchge->defineType("-", -1);

  blindMode = new RooCategory("blindMode", "Blind Status"); 
  blindMode->defineType("unblinded", 0);   
  blindMode->defineType("blinded",   1);
  blindMode->setIndex(0);
  blindMode->setConstant();

  noBlinding = new RooCategory("noBlinding", "No blinding"); 
  noBlinding->defineType("unblinded", 0);
  noBlinding->setIndex(0);
  noBlinding->setConstant();


  // set some binnings
  setBinning(Deltae, 0.005);
  setBinning(mes, 0.001);
  setBinning(bknnout, 0.05);
  setBinning(nnout, 0.05);
  setBinning(d0mass, 0.0025);

  usedRepCont = kFALSE;
  noNNPdf = kFALSE;
  noN2Pdf = kFALSE;
  noN1Pdf = kFALSE;
  useBBoAll = kFALSE;
  TightMesDe = kFALSE;
 
  allVars = new RooArgSet(*mes, *Deltae, *d0mass,
			  *nnout, *bknnout, *Hdtrkkaonnn,
			  *d0ppupmass,*B1decmode, 
		          *B2decmode, "allVars");

  allVars->add(*R2);			
  allVars->add(*Trued0flg);
  allVars->add(*d0recdec); 
  allVars->add(*Hdtrkchge);
  allVars->add(*mass12);
  allVars->add(*mass13);
  allVars->add(*mass12mc);
  allVars->add(*mass13mc);
  allVars->add(*mixneumass);
  allVars->add(*Pippidbit);
  allVars->add(*Pimpidbit);
  allVars->add(*d0ppupmass);
  allVars->add(*d0flightdist);
  allVars->add(*excltruth);

  allVars->Print("V");

  pidEffVars = new RooArgSet(*Pippeff,*Pimpeff,*Hdtrkpeff,"pidEffVars");

  // make sure nothing in allVars gets returned by parameters():
  BdkPdfAbsBase::setAllAnalysisVars(allVars);

  // Configure integrator for 1D projections 
  plot1dIntCfg.setEpsAbs(1E-5);
  plot1dIntCfg.setEpsRel(1E-5);
  plot1dIntCfg.method1D().setLabel("RooSegmentedIntegrator1D");
  
  // Read Dalitz plot configuration
  dalitzCfg = new BdkDalitzCfg("dalitzCfg","dalitzCfg");
  const char* cfgFile = "../BToDKTo3piK/params/dalitzCfg.par";
  cout << "Reading Daliz configuration from "<<cfgFile<<endl;
  dalitzCfg->getParameters(RooArgSet())->readFromFile(cfgFile);
  dalitzCfg->getParameters(RooArgSet())->Print("v");

  // setup m23 variable
  mtotal = new RooRealVar("mtotal","",dalitzCfg->M()*dalitzCfg->M() + 
                                      dalitzCfg->m1()*dalitzCfg->m1() +
                                      dalitzCfg->m2()*dalitzCfg->m2() + 
                                      dalitzCfg->m3()*dalitzCfg->m3());

  s23 = new RooFormulaVar("m23","@0-@1-@2",RooArgList(*mtotal,*m12,*m13));
  s23->setUnit("GeV^{2}/c^{4}");

  // m13(m13,m23)
  s13_23 = new RooFormulaVar("m13","@0-@1-@2",RooArgList(*mtotal,*m12,*m23));
  s13_23->setUnit("GeV^{2}/c^{4}");

  cout << "--- END setup() ---" << endl;
}


void setBinning(RooRealVar *r, Double_t gev)
{
  if (!r) return;
  // This avoids rounding problems
  r->setBins((Int_t)((10000*r->getMax()-10000*r->getMin())/(10000*gev)));
}
