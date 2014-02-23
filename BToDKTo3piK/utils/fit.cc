// Does a fit, implementing global flags:

#include "TString.h"

#include "RooFitCore/RooFitResult.hh"
#include "RooFitCore/RooAbsData.hh"
#include "RooFitCore/RooNLLVar.hh"
#include "RooFitCore/RooMinuit.hh"

#include "BToDKTo3piK/BdkPdfAbsBase.hh"
#include "BToDKTo3piK/BdkPdfOnRes.hh"

#include "../BToDKTo3piK/globals/globals.hh"
#include "../BToDKTo3piK/utils/printFitResult.cc"



// Check if an unblinded fit is allowed for the dataset
void setupBlinding(RooAbsData& dat)
{
  TString s1(dat.GetName());
  blindMode->setLabel("unblinded");
  /*  We are unblinded now!
  for (int i=0; i<dataTreeList->GetSize(); i++) {
    if (s1==TString(dataTreeList->At(i)->GetName()))
      blindMode->setLabel("blinded");
  }
  */

  cout << "==============================================================="<<endl;
  cout << "                Doing a "<<blindMode->getLabel()<<" fit!"<<endl;
  if (blindMode->getIndex()==0)
    cout<<"   This better be Monte Carlo or you have your eyes closed !  "<<endl;
  cout << "==============================================================="<<endl;
}


// Fit using a RooAbsPdf
RooFitResult * fit(RooAbsPdf & pdf, 
		   RooAbsData & dat,
		   RooAbsReal * penalty = 0,
		   Bool_t usepenaltyOnly = kFALSE) {

  if (doFit) {
   
    // Only allow unblinded fits on MC for now

    setupBlinding(dat);

    //    fitResult = pdf.getPdf()->fitTo(dat, (const char *)fitOption,
    //                                    (const char *)optOption);

    /////////// START: Copied from RooAbsPdf::fitTo() ///////////////////
    //
    // This is necessary because we want to be able to plot contours
    // using gMinuit->Contour() after the fit. For this the RooNLLVar must
    // not be deleted. Unfortunatelly this is the case with RooAbsPdf::fitTo().
    //
    TString fopt(fitOption);
    TString oopt(optOption);
    fopt.ToLower() ;
    oopt.ToLower() ;

    Bool_t extended = fopt.Contains("e") ;  
    Bool_t saveRes  = fopt.Contains("r") ;
    Bool_t cOpt     = oopt.Contains("p") || // for backward compatibility
                      oopt.Contains("c") ;
    Bool_t blindfit   = fopt.Contains("b") ;  


    Int_t  ncpu = 1 ;
    if (oopt.Contains("2")) ncpu=2 ;
    if (oopt.Contains("3")) ncpu=3 ;
    if (oopt.Contains("4")) ncpu=4 ;
    if (oopt.Contains("5")) ncpu=5 ;
    if (oopt.Contains("6")) ncpu=6 ;
    if (oopt.Contains("7")) ncpu=7 ;
    if (oopt.Contains("8")) ncpu=8 ;
    if (oopt.Contains("9")) ncpu=9 ;

    // Construct NLL
    if (nllVar) delete nllVar;
    nllVar = new RooNLLVar("nll","-log(likelihood)",pdf,dat,
                           RooArgSet(),extended,0,ncpu) ;

    RooAbsReal * finalNLL = nllVar;
    if (0 != penalty) {
      if (kFALSE == usepenaltyOnly) {
	finalNLL = new RooFormulaVar("nll+penalty", "nll+penalty",
				     TString(nllVar->GetName()) + "+" +
				     TString(penalty->GetName()), 
				     RooArgSet(*nllVar, *penalty));
	cout << "fit(): Fitting with penalty+nll" << endl;
      }
      else {
	finalNLL = penalty;
	cout << "fit(): Fitting with penalty only" << endl;
      }
    }

    // Minimize NLL
    if (minuit) delete minuit;
    minuit = new RooMinuit(*finalNLL) ;
    if(blindfit) minuit->setPrintLevel(-1);
    if (cOpt) minuit->optimizeConst(1) ;
    minuit->fit(fopt) ;
  
    // Optionally return fit result
    if (saveRes) fitResult = minuit->save() ;
    else fitResult = 0;

    //////////////////////// END //////////////////////////

    printFitResult(fitResult);
    return fitResult;
  }

  return 0;
}


// Fit using a BdkPdfAbsBase
RooFitResult * fit(BdkPdfAbsBase & pdf, 
		   RooAbsData & dat,
		   RooAbsReal * penalty = 0,
		   Bool_t usepenaltyOnly = kFALSE) {

  cout << "fit(): Fitting parameters of PDF \"" << pdf.GetName() 
       << "\":" << endl;
  pdf.parametersFree().Print("V");
  cout << "fit(): fitOption = \"" << fitOption << "\"" << endl
       << "       optOption = \"" << optOption << "\"" << endl;

  return fit(*pdf.getPdf(), dat, penalty, usepenaltyOnly);
}


// Do both fits (script interface to BdkPdfOnRes::fit())
// Only returns x/y fit result
// For yield fit result use pdfOnResDK.yieldFitResult() after fit
RooFitResult * fit(BdkPdfOnRes& pdf, 
                   RooAbsData& dat, 
                   Bool_t usePenalty, Bool_t usePenaltyOnly = kFALSE)
{
  if (doFit) {
    setupBlinding(dat);
    pdf.fit(dat, usePenalty, usePenaltyOnly);
    fitResult = pdf.xyFitResult();

    if (fitResult) printFitResult(fitResult);
    return fitResult;
  }
  return 0;
}


// Only do the x/y fit (script interface to BdkPdfOnRes::fit())
RooFitResult * fit(BdkPdfOnRes& pdf, 
                   RooAbsData& dat, 
                   RooFitResult& yieldFitResult,
                   Bool_t usePenalty = kTRUE, Bool_t usePenaltyOnly = kFALSE)
{
  if (doFit) {
    setupBlinding(dat);
    pdf.fit(dat, yieldFitResult, usePenalty, usePenaltyOnly);
    fitResult = pdf.xyFitResult();

    if (fitResult) printFitResult(fitResult);
    return fitResult;
  }
  return 0;
}



void useYieldFitVars()
{
  pdfOnResDK.useDE();
  pdfOnResDK.useNnCont();
  pdfOnResDK.useNnComb();
  pdfOnResDK.useDalitz(false); 

  pdfOnResDK.useMd(false);
  pdfOnResDK.useMes(false); 

  pdfOnResDK.getPdf();
}


void useXyFitVars()
{
  pdfOnResDK.useDE();
  pdfOnResDK.useNnCont();
  pdfOnResDK.useNnComb(false);
  pdfOnResDK.useDalitz(); 

  pdfOnResDK.useMd(false);
  pdfOnResDK.useMes(false); 

  pdfOnResDK.getPdf();
}

void useBothFitVars()
{
  pdfOnResDK.useDE();
  pdfOnResDK.useNnCont();
  pdfOnResDK.useNnComb();
  pdfOnResDK.useDalitz();

  pdfOnResDK.useMd(false);
  pdfOnResDK.useMes(false); 

  pdfOnResDK.getPdf();
}
