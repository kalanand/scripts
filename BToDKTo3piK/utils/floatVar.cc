// Floats a specific variable by name. If floatIt = kFALSE then fixes it:

#ifndef FLOATVAR_CC
#define FLOATVAR_CC

#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooRealVar.hh"

void floatVar(RooArgSet & argSet, const char * varName, Bool_t floatIt = kTRUE)
{
  RooRealVar * var = (RooRealVar *)argSet.find(varName);
  if (0 == var) {
    cout << "floatvar(): RooArgSet \"" << argSet.GetName() 
	 << "\" has no variable \"" << varName << "\". Doing nothing." << endl;
  }
  else {
    var->setConstant(!floatIt);
  }
}

void floatVar(RooArgSet * argSet, const char * varName, Bool_t floatIt = kTRUE)
{
  if (0 != argSet) {
    floatVar(*argSet, varName, floatIt);
  }
}

void fixVar(RooArgSet & argSet, const char * varName, Bool_t fixIt = kTRUE) {
  floatVar(argSet, varName, !fixIt);
}

void fixVar(RooArgSet * argSet, const char * varName, Bool_t fixIt = kTRUE) {
  floatVar(argSet, varName, !fixIt);
}

#endif
