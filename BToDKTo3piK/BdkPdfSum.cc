/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfSum.cc,v 1.1 2005/10/09 22:50:37 abi Exp $
 * Authors:
 *   Abi Soffer, Colorado State University, abi@slac.stanford.edu
 * Description:
 *   See .rdl file.
 * History:
 *   17-Apr-2002 abi Created initial version
 *   18-Sep-2002 WW Adjusted syntax to work with SunOS 5.8 compiler.
 *
 * Copyright (C) 2002 Colorado State University and SLAC
 *****************************************************************************/

// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
//
// This class provides a wrapper for the sum of pdfs.
//

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooRealVar.hh"
#include "BToDKTo3piK/BdkPdfSum.hh"

ClassImp(BdkPdfSum)

//-------------------------------------------------------------------
BdkPdfSum::BdkPdfSum() : _linkFracs(0) {}

//-------------------------------------------------------------------
BdkPdfSum::BdkPdfSum(const char * theName,
		     const char * theDesc,
		     BdkPdfAbsBase & pdf1,
		     BdkPdfAbsBase & pdf2,
		     NameScheme scheme,
		     BdkPdfSum * linkFracs) {
  // Constructor with 2 pdfs
  init(theName, theDesc, pdf1, pdf2, scheme, linkFracs);
}

//-------------------------------------------------------------------
BdkPdfSum::BdkPdfSum(const char * theName,
		     const char * theDesc,
		     TList & pdfList,
		     NameScheme scheme,
		     BdkPdfSum * linkFracs) {
  // Constructor with a list of pdfs
  init(theName, theDesc, pdfList, scheme, linkFracs);
}
  
//-------------------------------------------------------------------
BdkPdfSum::~BdkPdfSum() {
  _fractions.removeAll();
}
  
//-------------------------------------------------------------------
void BdkPdfSum::init(const char * theName,
		     const char * theDesc,
		     BdkPdfAbsBase & pdf1,
		     BdkPdfAbsBase & pdf2,
		     NameScheme scheme,
		     BdkPdfSum * linkFracs) {
  // Initializer with 2 PDFs:
  TList pdfList;
  pdfList.Add(&pdf1);
  pdfList.Add(&pdf2);
  init(theName, theDesc, pdfList, scheme, linkFracs);
}
  
//-------------------------------------------------------------------
void BdkPdfSum::init(const char * theName,
		     const char * theDesc,
		     TList & pdfList,
		     NameScheme scheme,
		     BdkPdfSum * linkFracs) {
  // Initializer with a list of PDFs:
  baseInit(theName, theDesc);
  setComponents(pdfList);
  _nameScheme = scheme;
  setLinkFracs(linkFracs);
}

//-------------------------------------------------------------------
void BdkPdfSum::initParameters() {
  // Not used, replaced essentially by makeFracms(NameScheme).
}

//-------------------------------------------------------------------
void BdkPdfSum::setLinkFracs(BdkPdfSum * link) {
  // Use the fractions of link. If link==0, makes its own fractions:
  _linkFracs = link;
  makeFracs(_linkFracs);
}

//-------------------------------------------------------------------
void BdkPdfSum::makeFracs(BdkPdfSum * theLinkFracs) {
  // Make the _fractions list.
  // First, empty the _fractions list:
  _fractions.removeAll();

  const int nFracs = numFractions();  // # of fractions needed

  // Check that _linkFracs has enough fractions:
  if (0 != _linkFracs) {
    if (_linkFracs->numComponents() < numComponents()) {
      cerr << GetName() 
	   << "::makeFracs(): Cannot use fractions of lincFracs \""
	   << _linkFracs->GetName() << "\"" << endl
	   << "  which has only "
	   << _linkFracs->numComponents() << " components while "
	   << numComponents() << " are needed." << endl;

      assert(0);
    }

    if (_linkFracs->numComponents() > numComponents()) {
      cerr << GetName() 
	   << "::makeFracs(): Warning: lincFracs \""
	   << _linkFracs->GetName() << "\"" << endl
	   << "  which has "
	   << _linkFracs->numComponents() << " components while only "
	   << numComponents() << " are needed. Will proceed." << endl;
    }
  } // end of checking _linkFracs->numComponents()

  // Loop over fractions needed:
  for (int c = 0; c < nFracs; ++c){
    if (0 != _linkFracs) {
      // Take the fractions from _linkFracs, and add it to the list as
      // non-owned:
      _fractions.add(*(_linkFracs->fraction(c)));
    }
    else {
      // Make and maintain our own fractions. First, initialize the
      // fraction's name:
      TString fracName = TString(GetName()) + ".";
      
      // The scheme determines the rest of the fracName:
      TString suffix;
      TString pdfName(component(c)->GetName());
      
      switch(_nameScheme) {
      case LAST:
	suffix = pdfName(pdfName.Last('.') + 1,  pdfName.Length() - 1);
	break;
      case FULL:
	suffix = pdfName;
	break;
      case NUMBER:
	suffix += c;
	break;
      }
      
      // Append the suffix to the name and create the title:
      fracName += suffix + "Frac";
      
      TString fracTitle = 
	TString(GetTitle()) + TString(" frac of ") + component(c)->GetTitle();
      
      // Create the fraction. The default value gives all components 
      // equal fractions:
      RooRealVar * frac = 
	new RooRealVar(fracName, fracTitle, 1.0/numComponents(), 0, 1);

      frac->setError(0.01);

      // Add the fraction to the list as owned:
      _fractions.addOwned(*frac);
    } // end of if no _linkFractions
  } // end loop on fractions needed
}

//-------------------------------------------------------------------
void BdkPdfSum::createPdf() {
  // Create _thePdf:
  _thePdf = new RooAddPdf(TString(GetName()) + ".pdf",
			  TString(GetTitle()) + " pdf",
			  componentPdfs(),
			  _fractions);

  setIsValid(kTRUE);
}
			  
//-------------------------------------------------------------------
BdkPdfSum::BdkPdfSum(const BdkPdfSum &) {
  // Forbidden copy constructor:
  assert(0);
}

//-------------------------------------------------------------------
BdkPdfSum & BdkPdfSum::operator=(const BdkPdfSum &) {
  // Forbidden assignment operator:
  assert(0);
  BdkPdfSum *pdf;
  return *pdf;
}

//-------------------------------------------------------------------
void BdkPdfSum::fixFractions(Bool_t fix) {
  for (int f = 0; f < numFractions(); ++f){
    fraction(f)->setConstant(fix);
  }
}
  












