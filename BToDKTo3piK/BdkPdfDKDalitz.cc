/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkPdfDKDalitz.cc,v 1.27 2007/08/22 07:59:05 fwinkl Exp $
 * Description:
 *   Signal Dalitz PDF wrapper
 * History:
 *   18 Oct 2005, created, Abi soffer
 *
 * Copyright (C) 2005 Colorado State University and SLAC
 *****************************************************************************/
// -- CLASS DESCRIPTION [BDKPDFWRAPPER] --
// 
// Wrapper for B->DK Dalitz plot PDF
// 

#include "TString.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooFormulaVar.hh"
#include "RooFitCore/RooArgList.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitModels/RooUnblindPrecision.hh"

#include "BToDKTo3piK/BdkMath.hh"
#include "BToDKTo3piK/BdkPdfAbsBase.hh"
#include "BToDKTo3piK/BdkDKDalitz.hh"
#include "BToDKTo3piK/BdkDKNonCDalitz.hh"
#include "BToDKTo3piK/BdkPdfDKDalitz.hh"
#include "BToDKTo3piK/BdkAbsDDalitzAmp.hh"


ClassImp(BdkPdfDKDalitz)


// Constructors:
BdkPdfDKDalitz::BdkPdfDKDalitz() {
  setIsValid(kFALSE);
}


BdkPdfDKDalitz::BdkPdfDKDalitz(const char * theName, const char * theDesc,
                               RooAbsReal & m12, RooAbsReal & m13,
                               BdkDalitzBase::Flavor flavor,
                               BdkPdfDKDalitz::COORD coord,
                               RooCategory* blindMode,
                               BdkAbsDDalitzAmp * externalAmp,
                               BdkDalitz::Mode DdecMode,
                               BdkAbsDDalitzAmp * externalAmpBar) :
  _x(0), _y(0),
  _rho(0), _theta(0),
  _x0(0), _y0(0),
  _deltaD(0)
{
  init(theName, theDesc, m12, m13, flavor, coord,
       blindMode, externalAmp, DdecMode, externalAmpBar);
}
  

// destructor:
BdkPdfDKDalitz::~BdkPdfDKDalitz() 
{
}

// initializer:  
void BdkPdfDKDalitz::init(const char * theName, const char * theDesc,
                          RooAbsReal & m12, RooAbsReal & m13,
                          BdkDalitzBase::Flavor flavor,
                          BdkPdfDKDalitz::COORD coord,
                          RooCategory* blindMode,
                          BdkAbsDDalitzAmp * externalAmp,
                          BdkDalitz::Mode DdecMode,
                          BdkAbsDDalitzAmp * externalAmpBar) 
{
  // base class initialization:
  BdkPdfDalitzBase::init(theName, theDesc, m12, m13, flavor, DdecMode);

  _amp = externalAmp;
  _ampBar = externalAmpBar;
  _blindMode = blindMode;

  if (_ampBar) _deltaD = new RooRealVar(TString(GetName()) + ".deltaD",
                                        TString(GetTitle()) + ".deltaD",
                                        0);
  
  // Need to initialize x0/y0 after init() manually
  // because dalitzAmp is not available at this point yet
  if (coord==BdkPdfDKDalitz::POLAR) usePolar(0,0);
  else useCartesian();
}


// Use cartesian coordiantes
void BdkPdfDKDalitz::useCartesian()
{
  _coord = CART;

  Double_t x = 0.1;
  Double_t y = 0.1;

  // Save a pointer to dalitzAmp
  // This is needed in createPdf to not duplicate BdkAbsDDalitzAmp
  if (_thePdf) _amp = ((BdkDKDalitz*)_thePdf)->dalitzAmp();

  // If we are coming from polar coordinates we can use 
  // the old values as initial values.
  if (_x!=0) x = _x->getVal();
  if (_y!=0) y = _y->getVal();
    
  delete _x;
  delete _y;
  // Create the RRV's:
  _x = new RooRealVar(TString(GetName()) + ".x",
                      TString(GetTitle()) + " x",
                      x);
  ((RooRealVar*)_x)->setConstant(kFALSE);
  ((RooRealVar*)_x)->setError(0.005);
  
  _y = new RooRealVar(TString(GetName()) + ".y",
                      TString(GetTitle()) + " y",
                      y);
  ((RooRealVar*)_y)->setConstant(kFALSE);
  ((RooRealVar*)_y)->setError(0.005);   

  // Delete polar coordinates
  delete _rho; _rho = 0;
  delete _theta; _theta = 0;
  delete _x0; _x0 = 0;
  delete _y0; _y0 = 0;

  setupBlinding();
  setIsValid(kFALSE);
}


// Use polar coordinates with origin at minimum of yield
void BdkPdfDKDalitz::usePolar()
{
  Double_t x0 = 0;
  Double_t y0 = 0;
  
  if (_thePdf) {
    x0 = pdfType()->dalitzAmp()->x0();    
    y0 = pdfType()->dalitzAmp()->y0();

    cout << GetName() << ": Using polar coordinates with origin at x0 = "<<x0
       << ", y0 = "<<y0 << endl;
  }

  usePolar(x0, y0);
}

// Use polar coordinates with origin at (x0,y0)
void BdkPdfDKDalitz::usePolar(Double_t x0, Double_t y0)
{
  _coord = POLAR;

  // Save a pointer to dalitzAmp
  // This is needed in createPdf to not duplicate BdkAbsDDalitzAmp
  if (_thePdf) _amp = ((BdkDKDalitz*)_thePdf)->dalitzAmp();

  // Center of the polar coordinates
  _x0 = new RooRealVar(TString(GetName()) + ".x0",
                       TString(GetTitle()) + " x0",
                       x0);

  _y0 = new RooRealVar(TString(GetName()) + ".y0",
                       TString(GetTitle()) + " y0",
                       y0);

  // Polar coordinates
  _rho = new RooRealVar(TString(GetName()) + ".rho",
                        TString(GetTitle()) + " rho",
                        0.1, 0, RooNumber::infinity);

  _theta = new RooRealVar(TString(GetName()) + ".theta",
                          TString(GetTitle()) + " theta",
                          0.1, 0, 360);

  // Initialize polar coordinates from cartesian coordinates
  if (_x!=0 && _y!=0) {
    Double_t dx = _x->getVal() - x0;
    Double_t dy = _y->getVal() - y0;
    ((RooRealVar*)_rho)->setVal(sqrt(dx*dx + dy*dy));
    Double_t theta = atan2(dy,dx);
    if (theta<0) theta += 2*TMath::Pi();
    ((RooRealVar*)_theta)->setVal(theta*TMath::RadToDeg());
  }

  delete _x;
  delete _y;
  _x = new RooFormulaVar(TString(GetName()) + ".x",
                         TString(GetTitle()) + " x",
                         "@0+@1*cos(@2*pi/180.0)", 
                         RooArgList(*_x0,*_rho,*_theta));

  _y = new RooFormulaVar(TString(GetName()) + ".y",
                         TString(GetTitle()) + " y",
                         "@0+@1*sin(@2*pi/180.0)", 
                         RooArgList(*_y0,*_rho,*_theta));

 
  setupBlinding();
  setIsValid(kFALSE);
}

void BdkPdfDKDalitz::recalcX0() {
  // Recalculate X0 after reading in new Dalitz shape parameters:
  dalitzAmp()->calDDbarNorm();
  
  _x0->setVal(pdfType()->dalitzAmp()->x0());
  
  cout << GetName() << ": Recalculated x0 = " << _x0->getVal() << endl;
}

void BdkPdfDKDalitz::setPolarCoords(double x, double y) {
  // Set the values of rho and theta from x and y:
  //  recalcX0();

  _rho->setVal(sqrt(sqr(x - _x0->getVal()) + sqr(y - _y0->getVal())));

  double thetaVal = atan2(y - _y0->getVal(), x - _x0->getVal()) / BdkAbsDDalitzAmp::DEGTORAD;
  if (thetaVal < 0) thetaVal += 360;
  if (thetaVal > 360) thetaVal -= 360;
  _theta->setVal(thetaVal);

  cout << GetName() << ": Setting polar coordinates rho="
       << _rho->getVal() << ", theta=" << _theta->getVal() 
       << " from x=" << x << ", y=" << y << endl;
}

// Set CP parameters (either cartesian or polar) from rB, gamma and delta
// A phase flip gamma -> -gamma is done for D0
void BdkPdfDKDalitz::setCPparams(Double_t rB, Double_t gamma, Double_t delta)
{
  if (flavor()==BdkDalitz::D0) gamma *= (-1);
  
  cout << GetName() << ": Setting CP parameters from rB=" << rB
       << ", gamma=" << gamma << ", delta=" << delta << endl;
  
  Double_t x = rB*cos((delta+gamma)*BdkAbsDDalitzAmp::DEGTORAD);
  Double_t y = rB*sin((delta+gamma)*BdkAbsDDalitzAmp::DEGTORAD);
  
  if (_coord==POLAR) {
    setPolarCoords(x,y);
  }
  else {
    cout << GetName()
         << ": Setting cartesian coordinates x=" << x << ", y=" << y
         << endl;
    ((RooRealVar*)_x)->setVal(x);
    ((RooRealVar*)_y)->setVal(y);    
  }
}


void BdkPdfDKDalitz::setupBlinding()
{
  if (_blindMode) {
    // Unblinded variables
    delete _xUnblind;
    delete _yUnblind;

    _xUnblind = new RooUnblindPrecision(TString(GetName())+".xUnblind",
                                        TString(GetTitle())+" x unblinded",
                                        "TheBlindingString",
                                        0, 0.1, *_x, *_blindMode); 
    
    _yUnblind = new RooUnblindPrecision(TString(GetName())+".yUnblind",
                                        TString(GetTitle())+" y unblinded",
                                        "TheBlindingString",
                                        0, 0.1, *_y, *_blindMode);
  }
}


// Build the PDF:
void BdkPdfDKDalitz::createPdf()
{
  RooAbsReal* x;
  RooAbsReal* y;

  if (_blindMode) {
    x = _xUnblind;
    y = _yUnblind;
  }
  else {
    x = _x;
    y = _y;
  }

  if (verbose().Contains("c")) {
    cout << GetName() << ".createPdf(): flavor = " << flavor()
         << ", DdecMode = " << getDdecMode()
         <<      endl;
  }

  RooAbsPdf * thePdf = 0;
  
  if (_ampBar==0) thePdf = new BdkDKDalitz(TString(GetName())+".pdf",
                                           TString(GetTitle())+" Pdf", 
                                           *_m12, *_m13, *x, *y,
                                           flavor(), getDdecMode(),
                                           _amp);
  
  else thePdf = new BdkDKNonCDalitz(TString(GetName())+".pdf",
                                    TString(GetTitle())+" Pdf", 
                                    *_m12, *_m13, *x, *y, *_deltaD,
                                    flavor(), getDdecMode(),
                                    _amp, _ampBar);

  // Save pointer to dalitzAmp for reuse on subsequent createPdf() calls
  _amp = ((BdkDKDalitz*)thePdf)->dalitzAmp();

  ((BdkDKDalitz*)thePdf)->setEfficiencyFunc(_effFunc);

  setPdf(*thePdf); 
  setIsValid(kTRUE);
}


//----------------------------------
BdkAbsDDalitzAmp * BdkPdfDKDalitz::dalitzAmp() {
  BdkDKDalitz * thePdf = (BdkDKDalitz *)getPdf();
  return (thePdf ? thePdf->dalitzAmp() : 0);
}

//----------------------------------
BdkAbsDDalitzAmp * BdkPdfDKDalitz::dalitzAmpBar() {
  return _ampBar;
}


//-------------------------------------
// Send verbosity flag to dalitzAmp:
void BdkPdfDKDalitz::setVerbose(const char * val) {
  BdkPdfAbsBase::setVerbose(val); // do base class action
  if (verbose().Contains("+")) {  // set to amp
    pdfType()->dalitzAmp()->setVerbose(val);
  }
}
