// declares global variables


#ifndef VARS_HH
#define VARS_HH

#include "TString.h"
#include "TFile.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooCategory.hh"
#include "RooFitCore/RooFormulaVar.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooFitResult.hh"
#include "RooFitCore/RooNumIntConfig.hh"
#include "RooFitCore/RooNLLVar.hh"
#include "RooFitCore/RooMinuit.hh"

// The variables:
RooRealVar * randAdd = 0;
RooRealVar * mes = 0;
RooRealVar * Deltae = 0;
RooRealVar * d0mass = 0;
RooRealVar * nnout = 0;
RooRealVar * bknnout = 0;
RooRealVar * Hdtrkkaonnn = 0;
RooRealVar * excltruth = 0;
RooRealVar * Pippidbit = 0;
RooRealVar * Pimpidbit = 0;
RooRealVar * Dalitz1 = 0;
RooRealVar * Dalitz2 = 0;
RooRealVar * d0ppmass = 0; 
RooRealVar * d0ppupmass = 0;
RooRealVar * d0recdec = 0;
RooRealVar * d0flightdist = 0;
RooRealVar * Trued0flg = 0;
RooRealVar * B1decmode = 0;
RooRealVar * B2decmode = 0;
RooCategory * Hdtrkchge = 0;
RooRealVar * Pippeff = 0;
RooRealVar * Pimpeff = 0;
RooRealVar * Hdtrkpeff = 0;
RooRealVar * mass12 = 0;
RooRealVar * mass13 = 0;
RooFormulaVar *s12 = 0;
RooFormulaVar *s13 = 0;
RooFormulaVar *s23 = 0;
RooFormulaVar *s13_23 = 0;
RooRealVar *m12 = 0;
RooRealVar *m13 = 0;
RooRealVar *m23 = 0;
RooRealVar * mtotal = 0;
RooRealVar * mass12mc = 0;
RooRealVar * mass13mc = 0;
RooFormulaVar *s12mc = 0;
RooFormulaVar *s13mc = 0;
RooRealVar *m12mc = 0;
RooRealVar *m13mc = 0;
RooRealVar * d0pPpMass = 0;
RooRealVar * d0pPmMass = 0;
RooRealVar * R2 = 0;
RooRealVar * mixneumass = 0;
RooFormulaVar *qprimeF = 0;
RooFormulaVar *dprimeF = 0;
RooRealVar *qprime = 0;
RooRealVar *dprime = 0;
RooFormulaVar* eventPidEff = 0;

Bool_t usedRepCont;
Bool_t noNNPdf;
Bool_t noN2Pdf;
Bool_t noN1Pdf;
Bool_t useBBoAll;
Bool_t TightMesDe;

RooCategory * blindMode = 0;
RooCategory * noBlinding = 0;

// some arg sets to hold them:
RooArgSet * allVars = 0;
RooArgSet * pidEffVars = 0;

// the data set for fitting:
RooDataSet * data = 0;

// the verbosity of the pdf for debugging:
TString verbosePdf = 0; // "vcidg+" for full verbosity

// A TFile to write stuff into:
TFile * tFile = 0;

// A global fit result:
RooFitResult * fitResult = 0;

// # of events to generate for toy MC. If 0 then figure out the # to generate
// from the sum of the numEvts:
int numEvtsToGenerate = 0;

// Reduced accuracy for 1D projections
RooNumIntConfig plot1dIntCfg(*RooAbsReal::defaultIntegratorConfig());

// These are used by utils/fit.cc
RooMinuit* minuit = 0;
RooNLLVar* nllVar = 0;

#endif
