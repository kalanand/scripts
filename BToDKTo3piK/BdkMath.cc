/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkMath.cc,v 1.1 2007/05/11 12:31:07 fwinkl Exp $
 * Description:
 *   Some math functions and macros
 * History:
 *   04 Apr 2007, created, Frank Winklmeier, 
 *
 * Copyright (C) 2007 Colorado State University and SLAC
 *****************************************************************************/

#include "BToDKTo3piK/BdkMath.hh"

#include "TMatrix.h"
#include "RooFitCore/RooFitResult.hh"
#include "RooFitCore/RooRealVar.hh"

//--------------------------------------------------------------------
// Fill the correlation and covariance matrix from the fit result
void BdkMath::getCovCorMatrix(RooFitResult* res, TMatrix& cov, TMatrix& cor)
{
  const int size = res->floatParsFinal().getSize();
  cov.ResizeTo(size, size);
  cor.ResizeTo(size, size);
  for (int i = 0; i < size; ++i) {
    RooRealVar * vari = (RooRealVar *)(res->floatParsFinal().at(i));
    for (int j = 0; j < size; ++j) {
      RooRealVar * varj = (RooRealVar *)(res->floatParsFinal().at(j));
      cor(i,j) = res->correlation(*vari, *varj);
      cov(i,j) = cor(i,j) * vari->getError() * varj->getError();
    }
  }
}
