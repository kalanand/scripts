// A class for deciding if to replace variables in MC with
// toy-generated variables for a specific set of event types given by
// typeBit.
//

#ifndef REPLACER_HH
#define REPLACER_HH

#include "RooFitCore/RooRealVar.hh"

struct Replacer {
  Replacer() :
    _var(0),
    _typeBit(0)
  {}

    
  Replacer(RooRealVar * var, int typeBit) :    
    _var(var),
    _typeBit(typeBit)
  {}


  // Data:
  
  RooRealVar * _var;
  int _typeBit;
};

#endif
