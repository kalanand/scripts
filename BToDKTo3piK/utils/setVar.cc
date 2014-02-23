// Set the value of a specific variable by name:

#ifndef SETVAR_CC
#define SETVAR_CC

#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooRealVar.hh"

void setVar(RooArgSet & argSet, const char * varName, double value)
{
  RooRealVar * var = (RooRealVar *)argSet.find(varName);
  if (0 == var) {
    cout << "setVar(): RooArgSet \"" << argSet.GetName() 
	 << "\" has no variable \"" << varName << "\". Doing nothing." << endl;
  }
  else {
    var->setVal(value);
  }
}

void setVar(RooArgSet * argSet, const char * varName, double value) 
{
  if (0 != argSet) {
    setVar(*argSet, varName, value);
  }

}

void incVar(RooArgSet & argSet, const char * varName, double value)
{
  RooRealVar * var = (RooRealVar *)argSet.find(varName);
  if (0 == var) {
    cout << "incVar(): RooArgSet \"" << argSet.GetName() 
	 << "\" has no variable \"" << varName << "\". Doing nothing." << endl;
  }
  else {
    var->setVal(var->getVal()+value);
  }
}

void incVar(RooArgSet * argSet, const char * varName, double value)
{
  if (0 != argSet) {
    incVar(*argSet, varName, value);
  }
}

// set the value of the var and fix it if it's 0:
void setVarFix0(RooArgSet & argSet, const char * varName, double value)
{
  RooRealVar * var = (RooRealVar *)argSet.find(varName);
  if (0 == var) {
    cout << "floatvar(): RooArgSet \"" << argSet.GetName() 
	 << "\" has no variable \"" << varName << "\". Doing nothing." << endl;
  }
  else {
    var->setVal(value);
    if (0 == value) {
      var->setConstant();
    }
  }
}

void setVarFix0(RooArgSet * argSet, const char * varName, double value)
{
  if (0 != argSet) {
    setVarFix0(*argSet, varName, value);
  }
}

// set the value of the var and fix it if it's 0 or 1:
void setVarFix0or1(RooArgSet & argSet, const char * varName, double value)
{
  RooRealVar * var = (RooRealVar *)argSet.find(varName);
  if (0 == var) {
    cout << "floatvar(): RooArgSet \"" << argSet.GetName() 
	 << "\" has no variable \"" << varName << "\". Doing nothing." << endl;
  }
  else {
    var->setVal(value);
    if (0 == value || 1 == value) {
      var->setConstant();
    }
  }
}

void setVarFix0or1(RooArgSet * argSet, const char * varName, double value)
{
  if (0 != argSet) {
    setVarFix0or1(*argSet, varName, value);
  }
}

#endif
