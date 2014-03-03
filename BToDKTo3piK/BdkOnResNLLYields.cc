/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkOnResNLLYields.cc,v 1.7 2006/06/23 20:44:52 fwinkl Exp $
 * Authors:
 *   Abi Soffer, Colorado State University, abi@slac.stanford.edu
 * Description:
 *   Applies the yields penalty on a BdkPdfOnRes
 *
 * History:
 *   23-may-2006 abi  Created initial version
 *
 * Copyright (C) 2006 Colorado State University and SLAC
 *****************************************************************************/

#include <iostream>
using namespace std;
#include <math.h>
#include "TMatrixD.h"
#include "RooFitCore/RooFitResult.hh"
#include "BToDKTo3piK/BdkOnResNLLYields.hh"
#include "BToDKTo3piK/BdkDKDalitz.hh"
#include "BToDKTo3piK/BdkEvtTypes.hh"

ClassImp(BdkOnResNLLYields)

//------------------------------------------
BdkOnResNLLYields::BdkOnResNLLYields(BdkPdfOnRes & pdf, 
				     double corrNsigAsym,
				     int systBit) :
  _pdf(&pdf),
  _cpProxies(GetName()+TString(".cpProxies"),
             GetTitle()+TString(" cpProxies"), this)
{
  SetName (TString(_pdf->GetName()) + ".NLLYields");
  SetTitle(TString(_pdf->GetTitle()) + ".NLLYields");

  // add all the dependents to the proxy list
  _cpProxies.add(_pdf->cpParams());
  _cpProxies.add(*_pdf->numEvt(BdkEvtTypes::SIG_GOOD_D));
  _cpProxies.add(*_pdf->typeAsym(BdkEvtTypes::SIG_GOOD_D));
                 

  // the expected yield PER CHARGE for no CP:
  _noCPNsig = 0.5 * pdf.nBB()->getVal() * BRBtoDK() * BRDto3piOverK2pi() * BRDtoK2pi() * efficiency();

  // the relative error on this expected yield is an incoherent sum
  // of the relative errors of the factors:
  _noCPNsigRelErr = 0;

  if (NBB & systBit) {
    _noCPNsigRelErr += nBBRelErr() * nBBRelErr();
  }

  if (BR3PIoverK2PI & systBit) {
    _noCPNsigRelErr += BRDto3piOverK2piRelErr() * BRDto3piOverK2piRelErr();
  }

  if (BRK2PI & systBit) {
    _noCPNsigRelErr += BRDtoK2piRelErr() * BRDtoK2piRelErr();
  }

  if (BRDK & systBit) {
    _noCPNsigRelErr += BRBtoDKRelErr() * BRBtoDKRelErr();
  }

  if (EFF & systBit) {
    _noCPNsigRelErr += efficiencyRelErr() * efficiencyRelErr();
  }

  _noCPNsigRelErr = sqrt(_noCPNsigRelErr);

  // Get parameters of the yields fit:
  RooRealVar* nsigVar= (RooRealVar *)(_pdf->numEvt(BdkEvtTypes::SIG_GOOD_D));
  RooRealVar* asymVar= (RooRealVar *)(_pdf->typeAsym(BdkEvtTypes::SIG_GOOD_D));

  _nsig = nsigVar->getVal();
  _asym = asymVar->getVal();

  // Use parabolic or average asymmetric errors for nSig and asym
  if (nsigVar->hasError()) _nsigErr = nsigVar->getError();
  else {
    cout << GetName() << ": Using average of asymmetric errors for "
         <<nsigVar->GetName() << endl;
    _nsigErr = (fabs(nsigVar->getAsymErrorLo())+fabs(nsigVar->getAsymErrorHi()))/2;
  }

  if (asymVar->hasError()) _asymErr = asymVar->getError();
  else {
    cout << GetName() << ": Using average of asymmetric errors for "
         <<asymVar->GetName() << endl;
    _asymErr = (fabs(asymVar->getAsymErrorLo())+fabs(asymVar->getAsymErrorHi()))/2;
  }

  // From the correlation matrix, we reconstruct the error matrix element:
  _errNsigAsym = corrNsigAsym * _asymErr * _nsigErr;

  // Make sure our value isn't cached but recalculated each time (it's cheap):
  setOperMode(RooAbsReal::ADirty);
}

//------------------------------------------
BdkOnResNLLYields::~BdkOnResNLLYields() {
}

//------------------------------------------
Double_t BdkOnResNLLYields::evaluate() const {
  // The error in the expected total yield comes from the relative
  // error on the nBB, etc. But these errors cancel in the expected asym:
  const double expectedNsigErr = expectedNsig() * _noCPNsigRelErr;
  const double expectedAsymErr = 0;

  // The observed yield, asymmetry, and their errors from the fit:
  TMatrixDSym errMatrix(2);
  errMatrix(0,0) = _nsigErr * _nsigErr + expectedNsigErr * expectedNsigErr;
  errMatrix(1,1) = _asymErr * _asymErr + expectedAsymErr * expectedAsymErr;
  errMatrix(1,0) = _errNsigAsym;
  errMatrix(0,1) = _errNsigAsym;

  // Invert the matrix:
  double det;
  TMatrixDSym invErrMatrix = errMatrix;
  invErrMatrix.InvertFast(&det);

  // Put the differences between the expected and measured nsig and
  // asymmetry into a vector:
  TMatrixD nsigAsymVector(2,1);
  nsigAsymVector(0,0) = _nsig - expectedNsig();
  nsigAsymVector(1,0) = _asym - expectedAsym();

  // Calculate the chi^2. We copy the invErrMatrix so as not to change it:
  TMatrixDSym chi2Matrix = 
    TMatrixDSym(invErrMatrix).SimilarityT(nsigAsymVector);

  assert(chi2Matrix.GetNoElements()==1);
  double result = chi2Matrix(0,0) / 2;   // divide by 2 for NLL

  if (verbose().Contains("v")) {
    cout << GetName() << ": errMatrix = ";
    errMatrix.Print();
    cout << "  invErrMatrix = ";
    invErrMatrix.Print();
    cout << " diffVector = ";
    nsigAsymVector.Print();
    cout << " value=" << result 
	 << " for ";
    printParams(cout);
  }
  return result;
}

//------------------------------------------
double BdkOnResNLLYields::expectedNsigP() const {
  BdkDKDalitz * pdfP = 
    (BdkDKDalitz *)(_pdf->sigGoodD0P().getDalitzPdf()->getPdf());

  return _noCPNsig * pdfP->normOverNoCP();
}

//------------------------------------------
double BdkOnResNLLYields::expectedNsigN() const {
  BdkDKDalitz * pdfN = 
    (BdkDKDalitz *)(_pdf->sigGoodD0N().getDalitzPdf()->getPdf());

  return _noCPNsig * pdfN->normOverNoCP();
}

//------------------------------------------
TObject* BdkOnResNLLYields::clone(const char* newname) const {
  return new BdkOnResNLLYields(*this);
}


//------------------------------------------
// The expected total yield:
double BdkOnResNLLYields::expectedNsig() const {
  return expectedNsigN() + expectedNsigP();
}

//------------------------------------------
// The expected total yield:
double BdkOnResNLLYields::expectedAsym() const {
  return (expectedNsigN() - expectedNsigP()) / expectedNsig();
}

//------------------------------------------
void BdkOnResNLLYields::printToStream(ostream& os, 
				      PrintOption opt,
				      TString indent) const {
  
  os << "BdkOnResNLLYields " << GetName() << ":" << endl
     << "  errNsigAsym=\t" << errNsigAsym() << endl
     << "  noCPNsig=\t" << noCPNsig() << endl
     << "  noCPNsigRelErr=\t" << noCPNsigRelErr() << endl
     << "  nsig=\t" << nsig() << endl
     << "  asym=\t" << asym() << endl
     << "  nsigErr=\t" << nsigErr() << endl
     << "  asymErr=\t" << asymErr() << endl
     << "  expectedNsigP=\t" << expectedNsigP() << endl
     << "  expectedNsigN=\t" << expectedNsigN() << endl  
     << "  expectedNsig=\t" << expectedNsig() << endl
     << "  expectedAsym=\t" << expectedAsym() << endl
     << "  getVal=\t" << getVal() << endl;
  printParams(os);
}

//------------------------------------------
void BdkOnResNLLYields::printParams(ostream& os) const {

  TIterator* iter = _cpProxies.createIterator();
  while (RooAbsReal* r = (RooAbsReal*)iter->Next()) {
    os << r->GetName() << "=" << r->getVal() << endl;
  }
  delete iter;
}

