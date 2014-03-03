/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkDalitzEff.cc,v 1.3 2006/02/20 03:44:08 fwinkl Exp $
 * Description:
 *   Base class for efficiency functions
 * History:
 *   Dec 1 2005, created, Abi Soffer
 *
 * Copyright (C) 2005 Colorado State University and SLAC
 *****************************************************************************/

#include "BToDKTo3piK/BdkDalitzEff.hh"
#include "RooFitCore/RooListProxy.hh"

/// constructor
BdkDalitzEff::BdkDalitzEff(const char *name, const char *title, 
                           BdkDalitzBase::Flavor flavor, 
                           BdkDalitzBase::Mode DdecMode) :
  BdkDalitzBase(name,title,flavor,DdecMode)
{
}

/// copy constructor
BdkDalitzEff::BdkDalitzEff(const BdkDalitzEff& other, const char* name) :
  BdkDalitzBase(other, name)
{
}

/// destructor
BdkDalitzEff::~BdkDalitzEff()
{
}
