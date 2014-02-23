/*************************
* class BdkPdfHistPdfHolder
**************************/

#include <iostream>
using std::cout;
using std::endl;


#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include "TString.h"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooDataHist.hh"
#include "BToDKTo3piK/BdkPdfHistPdfHolder.hh"

ClassImp (BdkPdfHistPdfHolder)

//------------------------------------------
BdkPdfHistPdfHolder::BdkPdfHistPdfHolder() 
{
  initHistNames();
}

BdkPdfHistPdfHolder::BdkPdfHistPdfHolder(const char * theName,
					 const char * theDesc,
					 RooRealVar * var,
					 const char * filename) 
{
  initHistNames();
  init(theName, theDesc, var, filename );
}  

BdkPdfHistPdfHolder::~BdkPdfHistPdfHolder() {
  delete _sigBadD0Hist;
  delete _sigGoodD0Hist;
  delete _DpiBadD0Hist;
  delete _DpiGoodD0Hist;
  delete _DPiXHist;
  delete _DKXHist;
  delete _BBbadD0Hist;
  delete _BBgoodD0Hist;
  delete _qqBadD0Hist;
  delete _qqGoodD0Hist;
}  

void BdkPdfHistPdfHolder::init(const char * theName,
                               const char * theDesc,
			       RooRealVar * var,
                               const char * filename) {  
  SetNameTitle(theName, theDesc);
//  nameInit(theName, theDesc);  
  setPdfHist(filename, var);
}

void BdkPdfHistPdfHolder::initHistNames() {
  setHistName(BdkEvtTypes::SIG_BAD_D, "DKBadD");
  setHistName(BdkEvtTypes::SIG_GOOD_D, "DKGoodD");
  setHistName(BdkEvtTypes::DPi_BAD_D, "DPiBadD");
  setHistName(BdkEvtTypes::DPi_GOOD_D, "DPiGoodD");
  setHistName(BdkEvtTypes::DPiX, "DPiX");
  setHistName(BdkEvtTypes::DKX, "DKX");
  setHistName(BdkEvtTypes::BB_BAD_D, "BBBadD");
  setHistName(BdkEvtTypes::BB_GOOD_D, "BBGoodD");
  setHistName(BdkEvtTypes::QQ_BAD_D, "qqBadD");
  setHistName(BdkEvtTypes::QQ_GOOD_D, "qqGoodD");	
}    

void BdkPdfHistPdfHolder::setPdfHist(const char * filename, RooRealVar * var) {
  if (0 == filename || 0 == var) {
    return;
  }
   
  // Build PDFs:
  
  TFile* file = new TFile(filename, "READ");
  if(0 == file ) {
    cout << " unable to open " << filename << endl;
    return;
  }
  
  TH1D * sigBd   = (TH1D*) file->Get(_histNames[BdkEvtTypes::SIG_BAD_D]);
  TH1D * sigGd   = (TH1D*) file->Get(_histNames[BdkEvtTypes::SIG_GOOD_D]);
  TH1D * DpiBd   = (TH1D*) file->Get(_histNames[BdkEvtTypes::DPi_BAD_D]);
  TH1D * DpiGd   = (TH1D*) file->Get(_histNames[BdkEvtTypes::DPi_GOOD_D]);
  TH1D * DpiX    = (TH1D*) file->Get(_histNames[BdkEvtTypes::DPiX]);
  TH1D * DKX     = (TH1D*) file->Get(_histNames[BdkEvtTypes::DKX]);
  TH1D * GBbd    = (TH1D*) file->Get(_histNames[BdkEvtTypes::BB_BAD_D]);
  TH1D * GBgd    = (TH1D*) file->Get(_histNames[BdkEvtTypes::BB_GOOD_D]);
  TH1D * qqBd    = (TH1D*) file->Get(_histNames[BdkEvtTypes::QQ_BAD_D]);
  TH1D * qqGd    = (TH1D*) file->Get(_histNames[BdkEvtTypes::QQ_GOOD_D]);
  
  RooDataHist * _sigBadD0Hist = 
    new RooDataHist("sigBd", "bad sig. ", *var, sigBd);
  
  _sigBadD0.init(TString(GetName())+".sigBadD0PdfHist", 
		 TString(GetTitle())+" sigBadD0PdfHist", *var, *_sigBadD0Hist);
  
  //
  RooDataHist * _sigGoodD0Hist = 
    new RooDataHist("sigGd", "good sig. ", *var, sigGd);
  
  _sigGoodD0.init(TString(GetName())+".sigGoodD0PdfHist", 
		  TString(GetTitle())+" siggGoodD0PdfHist", 
		  *var, *_sigGoodD0Hist);
  //
  RooDataHist * _DpiBadD0Hist = 
    new RooDataHist("DpiBd", "Bad Dpi ", *var, DpiBd);
  
  _DpiBadD0.init(TString(GetName())+".DpiBadD0PdfHist", 
		 TString(GetTitle())+" DpiBadD0PdfHis", *var, *_DpiBadD0Hist);
  //
  RooDataHist * _DpiGoodD0Hist = 
    new RooDataHist("DpiGd", "Good Dpi ", *var, DpiGd);
  
  _DpiGoodD0.init(TString(GetName())+".DpiGoodD0PdfHist", 
		  TString(GetTitle())+" DpiGoodD0PdfHis", *var, *_DpiGoodD0Hist);
  //
  RooDataHist * _DPiXHist = 
    new RooDataHist("DPiXPdf", "DPiX ", *var, DpiX);
  
  _DPiX.init(TString(GetName())+".DPiXPdfHist", 
	      TString(GetTitle())+" DPiXPdfHist", *var, *_DPiXHist);
  //
  RooDataHist * _DKXHist = 
    new RooDataHist("DKXPdf", "DKX ", *var, DKX);
  
  _DKX.init(TString(GetName())+".DKXPdfHist", 
		  TString(GetTitle())+" DKXPdfHist", 
		  *var, *_DKXHist);
  //
  RooDataHist * _BBbadD0Hist = 
    new RooDataHist("", "bad Gen.BB ", *var, GBbd);
  
  _BBbadD0.init(TString(GetName())+".BBbadD0PdfHist", 
		TString(GetName())+" BBbadD0PdfHist", *var, *_BBbadD0Hist);
  //
  RooDataHist * _BBgoodD0Hist = 
    new RooDataHist("GBgd", "good Gen.BB ", *var, GBgd);
  
  _BBgoodD0.init(TString(GetName())+".BBgoodD0PdfHist", 
		 TString(GetName())+" BBgoodD0PdfHist", *var, *_BBgoodD0Hist);
  //
  RooDataHist * _qqBadD0Hist = 
    new RooDataHist("qqBd", "Bad cont. ", *var, qqBd);
  
  _qqBadD0.init(TString(GetName())+".qqBadD0PdfHist", 
		TString(GetName())+" qqBadD0PdfHist", *var, *_qqBadD0Hist);
  //
  RooDataHist * _qqGoodD0Hist = 
    new RooDataHist("qqGd", "good cont. ", *var, qqGd);
  
  _qqGoodD0.init(TString(GetName())+".qqGoddD0PdfHist",
		 TString(GetTitle())+" qqGoodD0PdfHist", *var, *_qqGoodD0Hist);


  file->Close();
  delete file;
}



