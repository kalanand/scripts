/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfDalitzSmear.cc,v 1.1 2006/03/20 19:04:31 fwinkl Exp $
 * Description:
 *   BdkDalitzSmear PDF wrapper
 * History:
 *   17 Mar 2006, created, Frank Winklmeier
 *
 * Copyright (C) 2006 Colorado State University
 *****************************************************************************/

#include "TString.h"

#include "BToDKTo3piK/BdkPdfDalitzSmear.hh"

ClassImp(BdkPdfDalitzSmear)
  
/// Constructors:
BdkPdfDalitzSmear::BdkPdfDalitzSmear() {
  setIsValid(kFALSE);
}


BdkPdfDalitzSmear::BdkPdfDalitzSmear(const char * theName, const char * theDesc,
                                     BdkPdfDalitzBase& pdf,
                                     BdkPdfAbsBase& res12, BdkPdfAbsBase& res13)
{
    
  _pdf = &pdf;
  _res12 = &res12;
  _res13 = &res13;
  
  if (!res12.dependents().equals(RooArgSet(*pdf.m12()))) {
    cout << "BdkPdfDalitzSmear::BdkPdfDalitzSmear(" << GetName()
         << "m12 resolution pdf "<<res12.GetName()<<" does not depend on "
         << pdf.m12()->GetName()
         << endl;
  }
  if (!res13.dependents().equals(RooArgSet(*pdf.m13()))) {
    cout << "BdkPdfDalitzSmear::BdkPdfDalitzSmear(" << GetName()
         << "m13 resolution pdf "<<res13.GetName()<<" does not depend on "
         << pdf.m13()->GetName()
         << endl;
  }
  
  // base class initialization:
  BdkPdfDalitzBase::init(theName, theDesc, *pdf.m12(), *pdf.m13(), pdf.flavor());
  
  setIsValid(kFALSE);
}

/// destructor:
BdkPdfDalitzSmear::~BdkPdfDalitzSmear()
{
}


/// Build the PDF using the symmetric/asymmetric parametrization.
void BdkPdfDalitzSmear::createPdf() {

  BdkDalitzSmear* thePdf = new BdkDalitzSmear(TString(GetName())+".pdf",
                                              TString(GetTitle())+" Pdf",
                                              *(BdkDalitzBase*)_pdf->getPdf(),
                                              *_res12->getPdf(),
                                              *_res13->getPdf(),
                                              *_pdf->m12(), *_pdf->m13());
  setPdf(*thePdf);
  setIsValid(kTRUE);
}
