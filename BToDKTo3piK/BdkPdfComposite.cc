/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfComposite.cc,v 1.1 2005/10/09 22:50:36 abi Exp $
 * Authors:
 *   Abi Soffer, Colorado State University, abi@slac.stanford.edu
 * Description:
 *   Base class for composite Pdf wrappers. Subclasses need to use the 
 *   addComponent function to add their components. They also need to 
 *   implement BdkPdfAbsBase functionalities.
 * History:
 *   21-Mar-2002 abi Created initial version
 *   28-Mar-2002 ww  add dependents()
 *
 * Copyright (C) 2002 Colorado State University and SLAC
 *****************************************************************************/

// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
// 
// Wrapper for a composite pdf
// 

#include <iostream>
using std::cout;
using std::endl;

#include "BToDKTo3piK/BdkPdfComposite.hh"

#include "RooFitCore/RooAbsPdf.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooArgList.hh"
#include "RooFitCore/RooRealVar.hh"

ClassImp(BdkPdfComposite)


//--------------------------------------------------------------
BdkPdfComposite::BdkPdfComposite() {
  // Default constructor
}

//--------------------------------------------------------------
BdkPdfComposite::BdkPdfComposite(BdkPdfComposite & other) {
  // Copy constructor. Argument must be non-const in order for 
  // the components of its PDF to be used as non-const.  

  *this = other;
}

//--------------------------------------------------------------
BdkPdfComposite::BdkPdfComposite(const char * name, const char * desc) {
  // Constructor with strings only:
  baseInit(name, desc);
}

//--------------------------------------------------------------
BdkPdfComposite::BdkPdfComposite(const char * name, const char * desc, 
				 BdkPdfAbsBase & component) {
  // Constructor with 1 component:
  baseInit(name, desc);
  addComponent(component);
}

//--------------------------------------------------------------
BdkPdfComposite::BdkPdfComposite(const char * name, const char * desc, 
				 BdkPdfAbsBase & component1, 
				 BdkPdfAbsBase & component2) {
  // Constructor with 2 components:
  baseInit(name, desc);
  addComponent(component1);
  addComponent(component2);
}

//--------------------------------------------------------------
BdkPdfComposite::BdkPdfComposite(const char * name, const char * desc, 
				 TList & theComponents) {
  // Constructor with a list of components:
  baseInit(name, desc);
  setComponents(theComponents);
}

//--------------------------------------------------------------
BdkPdfComposite::~BdkPdfComposite() {
  // Destructor
}

//--------------------------------------------------------------
BdkPdfComposite & BdkPdfComposite::operator=(BdkPdfComposite & other) {
  // Assignment operator. Argument must be non-const in order for 
  // the components of its PDF to be used as non-const.
  setComponents(other.components());
  return *this;
}

//--------------------------------------------------------------
void BdkPdfComposite::setComponents(TList & theComponents) {
  // First remove all components:
  removeComponents();

  // Add all the components:
  TIterator * iter = theComponents.MakeIterator();
  while (BdkPdfAbsBase * comp = (BdkPdfAbsBase *)iter->Next()) {
    addComponent(*comp);
  }
}

//--------------------------------------------------------------
void BdkPdfComposite::addComponent(BdkPdfAbsBase & component) {
  // Add a component PDF and set not valid:
  _components.Add(&component);
  setIsValid(kFALSE);

  // Set verbosity of this component if recursive verbosity requested:
  if (verbose().Contains("+")) {
    component.setVerbose(verbose());
  }
}

//--------------------------------------------------------------
void BdkPdfComposite::addComponent(BdkPdfAbsBase * component) {
  if (0 == component) {
    cout << "BdkPdfComposite::addComponent(BdkPdfAbsBase * component):"
	 << endl << "  ignoring null component pointer." << endl;
  }
  else {
    addComponent(*component);
  }
}

//--------------------------------------------------------------
Bool_t BdkPdfComposite::getIsValid() const {
  // Overrides base class function. Checks its own _isValid, as well
  // as that of all components, returns true only if all are true:
  if (kFALSE == BdkPdfAbsBase::getIsValid()) {
    return kFALSE;
  }

  TIterator * iter = _components.MakeIterator();
  BdkPdfAbsBase * component;

  while(component = (BdkPdfAbsBase*)iter->Next()) {
    if (kFALSE == component->getIsValid()) {
      return kFALSE;
    }
  }

  return kTRUE;
}

//--------------------------------------------------------------
RooArgSet BdkPdfComposite::subDependents() {
  // Does nothing. May be overridden by subclasses:
  return RooArgSet();
}

//--------------------------------------------------------------
RooArgSet BdkPdfComposite::dependents() {
  // Returns all the dependents of this object:
  // First take the subclass's dependents:
  RooArgSet result = subDependents();

  // Then loop over the components:
  TIterator * iter = _components.MakeIterator();
  BdkPdfAbsBase * component;

  while(0 != (component = (BdkPdfAbsBase*)iter->Next())) {
    // Get the dependents of the component:
    RooArgSet deps = component->dependents();

    // And add them one by one to the result (repetitions are
    // automatically removed):
    TIterator * depsIter = deps.createIterator();
    depsIter->Reset();
    RooRealVar * var;
    while(0 != (var= (RooRealVar*)depsIter->Next())) {
      result.add(*var);
    }
    delete depsIter;
  }
      
  return result;
}

//--------------------------------------------------------------
RooArgList BdkPdfComposite::componentPdfs() {
  // Loop over the components and make a list from their RooAbsPdfs:
  RooArgList result;

  TIterator * iter = _components.MakeIterator();
  BdkPdfAbsBase * component;
  
  while(0 != (component = (BdkPdfAbsBase*)iter->Next())) {
    RooAbsPdf * pdf = component->getPdf();
    if (0 != pdf) {
      result.add(*pdf);
    }
  }

  return result;
}

//-------------------------------------------------------------------
int BdkPdfComposite::numComponents() const {
  TIterator * iter = ((BdkPdfComposite&)*this).compIter();
  int counter = 0;
  while (iter->Next()){
    ++counter;
  }
  return counter;
}

//-------------------------------------------------------------------
void BdkPdfComposite::setVerbose(const char * val) {
  // Set verbose for this PDF:
  BdkPdfAbsBase::setVerbose(val);
  
  // Determine if to set for all components too:
  if (verbose().Contains("+")) {
    TIterator * iter = _components.MakeIterator();
    while (BdkPdfAbsBase * comp = (BdkPdfAbsBase *)iter->Next()) {
      comp->setVerbose(val);
      
      if (verbose().Contains("i")) {
	cout << "setting verbose of " << comp->GetName() 
	     << " to " << val << endl;
      }
    }
  }
}
      
//-------------------------------------------------------------------
void BdkPdfComposite::removeComponents() {
  // Remove all components:
  if (verbose().Contains("d")) {
    cout << GetName() << "::removeComponents() called." << endl;
  }

  while (1) {
    TIterator * myIter = _components.MakeIterator();  
    if (BdkPdfAbsBase * comp = (BdkPdfAbsBase *)myIter->Next()) {
      _components.Remove(comp);
    }
    else {
      break;
    }
  }
}

