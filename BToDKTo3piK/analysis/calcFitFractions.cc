// $Id: calcFitFractions.cc,v 1.5 2006/06/03 02:33:12 fwinkl Exp $
// Script to calculate the fixed fit parameters from MC

#include "TTree.h"
#include "TCut.h"

#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooRealVar.hh"

#include "../BToDKTo3piK/globals/globals.hh"


Double_t sqr(Double_t x)
{
  return x*x;
}

// Count events in tree (Poisson errors)
void count(RooRealVar& result, TCut cut, TTree* tree)
{
  TTree *t = tree->CopyTree(cut);
  result.setVal(t->GetEntries());
  result.setError(sqrt(result.getVal()));
  delete t;
}

// Calculate expected number of events
void expected(RooRealVar& result, TCut cut, TTree* tree)
{
  count(result, cut, tree);
  result.setVal(result.getVal()*getTreeWeight(tree));
  result.setError(result.getError()*getTreeWeight(tree));
}

// Calculate expected number of events from two trees
void expected(RooRealVar& result, TCut cut, TTree* tree1, TTree* tree2)
{
  RooRealVar r1("r1","",0);
  RooRealVar r2("r2","",0);
  expected(r1,cut,tree1);
  expected(r2,cut,tree2);
  add(result,r1,r2);
}


// Default: Both R16 and R18
// For R16 do: calcFitFractions(sigTree16,b0Tree16,bpTree16,udsTree16,ccTree16,dpiTree16)
// For R18 do: calcFitFractions(sigTree18,b0Tree18,bpTree18,udsTree18,ccTree18,dpiTree18)
void calcFitFractions(TCut cut = cutSigReg, TTree* sig = sigTree,
		      TTree* b0 = b0Tree, TTree* bp = bpTree,
		      TTree* uds = udsTree, TTree* cc = ccTree,
		      TTree* dpi = dpiTree)
{
  gROOT->cd();

  RooArgSet set;
  RooRealVar sigGoodD0("sigGoodD0","",0); set.add(sigGoodD0);
  RooRealVar sigBadD0("sigBadD0","",0); set.add(sigBadD0);
  RooRealVar sigBadD0Frac("sigBadD0Frac","",0); set.add(sigBadD0Frac);
  expected(sigGoodD0, cut+cutDKGoodD, sig);
  expected(sigBadD0, cut+cutDKBadD, sig);
  divide(sigBadD0Frac, sigBadD0, sigGoodD0);

  RooRealVar DpiGoodD0("DpiGoodD0","",0); set.add(DpiGoodD0);
  RooRealVar DpiBadD0("DpiBadD0","",0); set.add(DpiBadD0);
  RooRealVar DpiBadD0Frac("DpiBadD0Frac","",0); set.add(DpiBadD0Frac);
  expected(DpiGoodD0, cut+cutDPiGoodD, dpi);
  expected(DpiBadD0, cut+cutDPiBadD, dpi);
  divide(DpiBadD0Frac, DpiBadD0, DpiGoodD0);

  RooRealVar DKX("DKX","",0); set.add(DKX);
  RooRealVar DPiX("DPiX","",0); set.add(DPiX);
  expected(DKX, cut+cutDKX, b0, bp);
  expected(DPiX, cut+cutDPiX, b0, bp);

  RooRealVar BBGoodD0("BBGoodD0","",0); set.add(BBGoodD0);
  RooRealVar BBBadD0("BBBadD0","",0); set.add(BBBadD0);
  expected(BBGoodD0, cut+cutBBGoodD, b0, bp);
  expected(BBBadD0, cut+cutBBBadD, b0, bp);

  RooRealVar totBBNumEvts("totBBNumEvts","",0); set.add(totBBNumEvts);
  add(totBBNumEvts, DKX, DPiX, BBBadD0);

  RooRealVar DKXFrac("DKXFrac","",0); set.add(DKXFrac);
  RooRealVar DPiXFrac("DPiXFrac","",0); set.add(DPiXFrac);
  divide(DKXFrac, DKX, DPiX);
  divide(DPiXFrac, DPiX, totBBNumEvts);

  RooRealVar BBGoodD0Frac("BBGoodD0Frac","",0); set.add(BBGoodD0Frac);
  divide(BBGoodD0Frac, BBGoodD0, totBBNumEvts);

  RooRealVar qqGoodD0("qqGoodD0","",0); set.add(qqGoodD0);
  RooRealVar qqBadD0("qqBadD0","",0); set.add(qqBadD0);
  RooRealVar qqGoodD0Frac("qqGoodD0Frac","",0); set.add(qqGoodD0Frac);
  expected(qqGoodD0, cut+cutqqGoodD, uds, cc);
  expected(qqBadD0, cut+cutqqBadD, uds, cc);
  divide(qqGoodD0Frac, qqGoodD0, qqBadD0);

  set.Print("v");
  set.writeToFile("calcFitFractions.par");
  set.printLatex(Format("e", AutoPrecision(1)));
}
