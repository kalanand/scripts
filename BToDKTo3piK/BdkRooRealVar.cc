/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package:
 *    File: $Id: BdkRooRealVar.cc,v 1.1 2006/06/12 22:42:02 fwinkl Exp $
 * Authors:
 *   Frank Winklmeier, Colorado State University, fwinkl@slac.stanford.edu
 * Description:
 *   Sortable RooRealVar
 *
 * History:
 *   12-Jun-2006 fwinkl  Created initial version
 *
 * Copyright (C) 2006 Colorado State University and SLAC
 *****************************************************************************/

#include "BToDKTo3piK/BdkRooRealVar.hh"

#include "RooFitCore/RooAbsArg.hh"

ClassImp(BdkRooRealVar)

BdkRooRealVar::BdkRooRealVar(const char *name, const char *title,
                             Double_t value, const char *unit) :
  RooRealVar(name, title, value, unit)
{
}

BdkRooRealVar::BdkRooRealVar(const char *name, const char *title, Double_t minValue, 
                             Double_t maxValue, const char *unit) :
  RooRealVar(name, title, minValue, maxValue, unit)
{
}

BdkRooRealVar::BdkRooRealVar(const char *name, const char *title, Double_t value, 
                             Double_t minValue, Double_t maxValue, const char *unit) :
  RooRealVar(name, title, value, minValue, maxValue, unit)
{
}

BdkRooRealVar::BdkRooRealVar(const RooRealVar& other, const char* name) :
  RooRealVar(other, name)
{
}

// Compare ourself to obj
// If obj is not a RooRealVar call parent Compare()
Int_t BdkRooRealVar::Compare(const TObject* obj) const
{
  if (const RooRealVar* r = dynamic_cast<const RooRealVar*>(obj)) {
    if (getVal() > r->getVal()) return 1;
    else if (getVal() < r->getVal()) return -1;
    else return 0;
  }
  else return RooAbsArg::Compare(obj);
}
