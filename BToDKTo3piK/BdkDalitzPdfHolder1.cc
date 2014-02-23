/*************************
* class BdkDalitzPdfHolder1
**************************/

#include <iostream>
using std::cout;
using std::endl;


#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooDataHist.hh"
#include "BToDKTo3piK/BdkDalitzPdfHolder1.hh"

ClassImp (BdkDalitzPdfHolder1)

//------------------------------------------
BdkDalitzPdfHolder1::BdkDalitzPdfHolder1() : 
  _ownedDPiXHist(0),
  _ownedSigBadD0Hist(0)
{
}

BdkDalitzPdfHolder1::BdkDalitzPdfHolder1(const char * theName,
					 const char * theDesc,
					 RooAbsReal & m12,
					 RooAbsReal & m13,
					 BdkDalitzBase::Flavor flavor,
                                         BdkPdfDKDalitz::COORD coord,
                                         BdkDalitzEff* effDKGoodD,
                                         BdkDalitzEff* effOther,
                                         const char * DPiXfileName,
                                         const char * sigBadD0fileName,
                                         RooCategory& blindMode,
                                         BdkDalitzPdfHolder1* pdfSource,
                                         Bool_t flipQqGoodDFlavor) :
  _ownedDPiXHist(0),
  _ownedSigBadD0Hist(0)
{

  init(theName, theDesc, m12, m13, flavor, coord,
       effDKGoodD, effOther,
       DPiXfileName, sigBadD0fileName, blindMode,
       pdfSource, flipQqGoodDFlavor);
}  

BdkDalitzPdfHolder1::~BdkDalitzPdfHolder1() {
  delete _ownedDPiXHist;
  delete _ownedSigBadD0Hist;
}  

//--------------------------------------------
void BdkDalitzPdfHolder1::init(const char * theName,
                               const char * theDesc,
			       RooAbsReal & m12,
			       RooAbsReal & m13,
			       BdkDalitzBase::Flavor flavor, 
                               BdkPdfDKDalitz::COORD coord,
                               BdkDalitzEff* effDKGoodD,
                               BdkDalitzEff* effOther,
                               const char * DPiXfileName,
                               const char * sigBadD0fileName,
                               RooCategory& blindMode,
                               BdkDalitzPdfHolder1* pdfSource,
                               Bool_t flipQqGoodDFlavor) {

  SetNameTitle(theName, theDesc);
  _m12 = &m12;
  _m13 = &m13;
  _flavor = flavor;
  _effDKGoodD = effDKGoodD;
  _effOther = effOther;
  _blindMode = &blindMode;

  BdkDDalitzAmp* externalAmp = 0;
  if (pdfSource != 0) externalAmp = (BdkDDalitzAmp*)pdfSource->sigGoodD0Type().dalitzAmp();

  _sigGoodD0.init(TString(GetName()) + ".sigGoodD0", 
		  TString(GetTitle()) + " sigGoodD0",
		  *_m12, *_m13, _flavor, coord, _blindMode, externalAmp);

  _sigGoodD0.setEfficiencyFunc(effDKGoodD);

  _DpiGoodD0.init(TString(GetName()) + ".DpiGoodD0", 
		  TString(GetTitle()) + " DpiGoodD0",
		  *_m12, *_m13, _flavor,(BdkDDalitzAmp*)_sigGoodD0.dalitzAmp());

  _DpiGoodD0.setEfficiencyFunc(effDKGoodD);

  _BBgoodD0.init(TString(GetName()) + ".BBGoodD0", 
		 TString(GetTitle()) + " BBGoodDo",
		 *_m12, *_m13, _flavor,(BdkDDalitzAmp*)_sigGoodD0.dalitzAmp());

  _BBgoodD0.setEfficiencyFunc(effDKGoodD);

  BdkDalitzBase::Flavor qqGoodDFlavor = _flavor;

  if (flipQqGoodDFlavor) {
    if (_flavor==BdkDalitzBase::D0) qqGoodDFlavor = BdkDalitzBase::D0BAR;
    else qqGoodDFlavor = BdkDalitzBase::D0;
  }
  _qqGoodD0.init(TString(GetName()) + ".qqGoodD0", 
		 TString(GetTitle()) + " qqGoodD0",
		 *_m12, *_m13, qqGoodDFlavor,(BdkDDalitzAmp*)_sigGoodD0.dalitzAmp());

  _qqGoodD0.setEfficiencyFunc(effDKGoodD);

  // The incoherent PDFs all have "spin-0" rhos. The resonant
  // components are instantiated once for _DpiBadD0, and then reused
  // in the other bad D PDFs:
  const int incRhoSpin = 0;

  // If we have a PDF source use the same DpiBad PDF
  _DpiBadD0.init(TString(GetName()) + ".DpiBadD0", 
                 TString(GetTitle()) + " DpiBadD0",
                 *_m12, *_m13, _flavor, incRhoSpin,
                 pdfSource ? &pdfSource->DpiBadD0Type() : 0);

  _DpiBadD0.setEfficiencyFunc(effOther);

  _DKX.init(TString(GetName()) + ".DKX", 
            TString(GetTitle()) + " DKX",
            _DpiBadD0);

  _DKX.setEfficiencyFunc(effOther);

  _BBbadD0.init(TString(GetName()) + ".BBBadD0", 
		TString(GetTitle()) + " BBBadD0",
                _DpiBadD0);
  
  _BBbadD0.setEfficiencyFunc(effOther);

  _qqBadD0.init(TString(GetName()) + ".qqBadD0", 
		TString(GetTitle()) + " qqBadD0",
                _DpiBadD0);

  _qqBadD0.setEfficiencyFunc(effOther);


  initSigBadD0Hist(sigBadD0fileName);
  initDPiXHist(DPiXfileName);
}

//--------------------------------------------
void BdkDalitzPdfHolder1::initDPiXHist(const char * filename) {
  // Build DPiX HistPdf:  
  if (0 == filename) {
    return;
  }

  TFile f(filename, "READ");
  if(f.IsZombie()) {
    cout << " unable to open " << filename << endl;
    return;
  }

  _ownedDPiXHist = (TH2*)f.Get("hist_dpix");
  _ownedDPiXHist->SetDirectory(gROOT);
  f.Close();

  setDPiXHist(*_ownedDPiXHist);
}


//--------------------------------------------
void BdkDalitzPdfHolder1::setDPiXHist(const TH2& h2) 
{
  if (BdkDalitzBase::D0 == _flavor)
    _DPiX.init(TString(GetName())+".DPiX", 
	       TString(GetTitle())+" DPiX",
	       *_m12, *_m13, BdkDalitzBase::D0, h2);
  else
    _DPiX.init(TString(GetName())+".DPiX", 
	       TString(GetTitle())+" DPiX",
	       *_m13, *_m12, BdkDalitzBase::D0BAR, h2);
  
  // Restore uniform binning:
  if (dynamic_cast<RooRealVar*>(_m12)) ((RooRealVar*)_m12)->setBins(100);
  if (dynamic_cast<RooRealVar*>(_m13)) ((RooRealVar*)_m13)->setBins(100);
}



//--------------------------------------------
void BdkDalitzPdfHolder1::initSigBadD0Hist(const char * filename)
{
   // Build sigBadD0 HistPdf:  
  if (0 == filename) {
    return;
  }

  TFile f(filename, "READ");
  if(f.IsZombie()) {
    cout << " unable to open " << filename << endl;
    return;
  }

  _ownedSigBadD0Hist = (TH2*)f.Get("hist_dkbadd");
  _ownedSigBadD0Hist->SetDirectory(gROOT);
  f.Close();

  setSigBadD0Hist(*_ownedSigBadD0Hist);
}


//--------------------------------------------
void BdkDalitzPdfHolder1::setSigBadD0Hist(const TH2& h2) 
{
  if (BdkDalitzBase::D0 == _flavor)
    _sigBadD0.init(TString(GetName())+".sigBadD0", 
		   TString(GetTitle())+" sigBadD0",
		   *_m12, *_m13, BdkDalitzBase::D0, h2);
  else
    _sigBadD0.init(TString(GetName())+".sigBadD0", 
		   TString(GetTitle())+" sigBadD0",
		   *_m13, *_m12, BdkDalitzBase::D0, h2);
  
  // Restore uniform binning:
  if (dynamic_cast<RooRealVar*>(_m12)) ((RooRealVar*)_m12)->setBins(100);
  if (dynamic_cast<RooRealVar*>(_m13)) ((RooRealVar*)_m13)->setBins(100);
}
