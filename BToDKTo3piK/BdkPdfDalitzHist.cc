/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfDalitzHist.cc,v 1.4 2006/05/22 21:53:06 fwinkl Exp $
 * Description:
 *   Wrapper class for BdkDalitzHist
 * History:
 *   Mar 22, 2006     Frank Winklmeier, created
 *
 * Copyright (C) 2006 Colorado State University and SLAC
 *****************************************************************************/

#include "RooFitCore/RooDataHist.hh"
#include "RooFitCore/RooArgSet.hh"

#include "BToDKTo3piK/BdkPdfDalitzHist.hh"
#include "BToDKTo3piK/BdkDalitzBase.hh"
#include "BToDKTo3piK/BdkDalitzHist.hh"

ClassImp(BdkPdfDalitzHist);

BdkPdfDalitzHist::BdkPdfDalitzHist() {
}

BdkPdfDalitzHist::BdkPdfDalitzHist(const char * theName,
                                   const char * theDesc,
                                   RooAbsReal &m12, RooAbsReal& m13,
                                   BdkDalitzBase::Flavor flavor,
				   const TH2& hist)
{
  init(theName, theDesc, m12, m13, flavor, hist);
}


void BdkPdfDalitzHist::init(const char * theName,
                            const char * theDesc,
                            RooAbsReal &m12, RooAbsReal& m13,
                            BdkDalitzBase::Flavor flavor,
			    const TH2& hist)
{
  BdkPdfDalitzBase::init(theName, theDesc, m12, m13, flavor);
  setDataHist(hist);
}


BdkPdfDalitzHist::~BdkPdfDalitzHist() {} //Destructor

void BdkPdfDalitzHist::setDataHist(const TH2& hist) {
  _hist =(TH2*) &hist;
  setIsValid(kFALSE);
} 

void BdkPdfDalitzHist::createPdf() {
//
// Build the HistPdf
//
  _thePdf = new BdkDalitzHist(TString(GetName())+".pdf",
			      TString(GetTitle())+" Pdf",
			      flavor(), BdkDalitzBase::PPP0,
			      *_m12, *_m13, *_hist);
                             
   setIsValid(kTRUE);
}
