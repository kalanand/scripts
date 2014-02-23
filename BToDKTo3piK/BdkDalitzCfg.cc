/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkDalitzCfg.cc,v 1.1 2006/02/20 03:44:08 fwinkl Exp $
 * Description:
 *    Singleton class to hold Dalitz plot configuration
 * History:
 *   Feb 19 2006, created, Frank Winklmeier
 *
 * Copyright (C) 2006 Colorado State University and SLAC
 *****************************************************************************/

#include "BToDKTo3piK/BdkDalitzCfg.hh"
#include "RooFitCore/RooListProxy.hh"
#include "RooFitCore/RooRealVar.hh"

/// constructor
BdkDalitzCfg::BdkDalitzCfg(const char *name, const char *title) :
  BdkDalitzBase(name,title,BdkDalitzBase::D0,BdkDalitzBase::PPP0)
{

  _m23VetoMass = new RooRealVar(TString(GetName())+".m23VetoMass",
				TString(GetName())+" m23 mass veto window",0);

  // link our parameters to the pdf:
  RooListProxy * proxyList =  new RooListProxy(TString(GetName()) + ".proxyList",
                                               TString(GetTitle()) + " proxyList",
                                               this);

  proxyList->add(*_m23VetoMass);
}

/// copy constructor
BdkDalitzCfg::BdkDalitzCfg(const BdkDalitzCfg& other, const char* name) :
  BdkDalitzBase(other, name)
{
}

/// destructor
BdkDalitzCfg::~BdkDalitzCfg()
{
}

/// set m23 mass veto window
/// min and max are the unsquared m23 limits
void BdkDalitzCfg::setM23VetoMass(Double_t min, Double_t max)
{
  _m23VetoMass->setRange(min,max);
}
