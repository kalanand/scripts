/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfSmear.cc,v 1.2 2006/03/17 01:13:57 fwinkl Exp $
 * Description:
 *   BdkSmear PDF wrapper
 * History:
 *   15 Mar 2006, created, Frank Winklmeier
 *
 * Copyright (C) 2006 Colorado State University
 *****************************************************************************/

#include "TString.h"

#include "BToDKTo3piK/BdkPdfSmear.hh"

ClassImp(BdkPdfSmear)
  
/// Constructors:
BdkPdfSmear::BdkPdfSmear() {
  setIsValid(kFALSE);
}


BdkPdfSmear::BdkPdfSmear(const char * theName, const char * theDesc, BdkPdfAbsBase& pdf)
{
  baseInit(theName, theDesc);

  _pdf = &pdf;
  _pdfVars.add(pdf.dependents());
  
  setIsValid(kFALSE);
}

/// destructor:
BdkPdfSmear::~BdkPdfSmear()
{
}

/// Add a resolution model
Bool_t BdkPdfSmear::addResolutionModel(BdkPdfAbsBase& model)
{
  if (model.dependents().getSize()!=1) {
    cout << "BdkPdfSmear::addResolutionModel(" << GetName()
	 << ") pdf " << model.GetName()
	 << " has more or less than one dependent:"
	 << endl;
    model.dependents().Print("v");
    return kFALSE;
  }
  _resModels.add(*model.getPdf());
  _resVars.add(model.dependents());

  setIsValid(kFALSE);
  
  return kTRUE;
}

/// Remove all resolution models
void BdkPdfSmear::removeResolutionModels()
{
  _resModels.removeAll();
  _resVars.removeAll();
  
  setIsValid(kFALSE);
}

/// Build the PDF using the symmetric/asymmetric parametrization.
void BdkPdfSmear::createPdf() {

  RooArgList models;
  RooArgList modelVars;
  
  BdkSmear * thePdf = new BdkSmear(TString(GetName())+".pdf",
				   TString(GetTitle())+" Pdf",
				   *_pdf->getPdf(), _pdfVars,
				   _resModels, _resVars); 
	  
  setPdf(*thePdf);
  setIsValid(kTRUE);
}
