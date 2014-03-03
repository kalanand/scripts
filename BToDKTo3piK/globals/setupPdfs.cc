// defines the global PDFs and their parameters

#ifndef SETUPPDFS_CC
#define SETUPPDFS_CC

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooRandom.hh"

#include "BToDKTo3piK/BdkDalitzBase.hh"
#include "BToDKTo3piK/BdkPdf2Gauss.hh"
#include "BToDKTo3piK/BdkPdfSum.hh"

#include "../BToDKTo3piK/globals/globals.hh"
#include "../BToDKTo3piK/utils/printEff.cc"

TString paramsPath = "../BToDKTo3piK/params/";

// forward
BdkPdfAbsBase* BdkPdfDalitzRes(const char* name, const char* title, RooRealVar& x);


// Read the parameters of pdfOnResDK
// Set fixAllExceptNumEvts to false if you want the original fit configuration
// that was used to obtain the individual shapes.
// Otherwise only parameters in numEvts.par are allowed to float.
void readOnResDKPar(Bool_t fixAllExceptNumEvts = kTRUE) {
 
  RooArgSet* pars = 0;
  TIterator* iter = 0;
  RooAbsArg *arg = 0;  

  pdfOnResDK.fixAll(false);    // unfix all parameters

  // Set all parameters to a random value.
  // We compare this value to the one after reading all.par to make sure
  // we initialized all parameters.
  pars = pdfOnResDK.getPdf()->getDependents(pdfOnResDK.parameters());
  iter = pars->createIterator();
  Double_t rnd = RooRandom::uniform();
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())=="RooRealVar") {
      RooRealVar *r = (RooRealVar*)arg;
      if (r == mtotal) continue;
      r->setVal(r->getMin()+(fabs(r->getMax())+fabs(r->getMin()))*rnd);
    }
  }
  // Make a backup of these random values
  RooArgSet* parsBeforeInit = (RooArgSet*)pdfOnResDK.parameters().snapshot(kFALSE);
  delete pars;
  delete iter;

  // Read the parameters      
  pdfOnResDK.parameters().readFromFile(paramsPath+"all.par");

  // This is a workaround for a RooFit bug in reading/writing non-existing fit 
  // http://babar-hn.slac.stanford.edu:5090/HyperNews/get/roothelp/1123.html
  // Remove all fit ranges of floating real parameters
  pars = (RooArgSet*)pdfOnResDK.parametersFree().snapshot(kFALSE);
  iter = pars->createIterator();  
  while ((arg = (RooAbsArg*)iter->Next())) {
    TString className = TString(arg->ClassName());
    if (className=="RooRealVar")
      ((RooRealVar*)arg)->removeRange();
    else if (className=="RooCategory")
      ((RooCategory*)arg)->setConstant();
  }
  // Now read the parameters again
  pdfOnResDK.parameters().readFromFile(paramsPath+"all.par");

  if (fixAllExceptNumEvts) {
    // Fix all parameters and only allow the ones in numEvts.par to float
    pdfOnResDK.fixAll();
    pdfOnResDK.parameters().readFromFile(paramsPath+"numEvts.par");
  }

  // Set nBB manually from the value in chains.hh
  pdfOnResDK.nBB()->setVal(N_BB);
  pdfOnResDK.nBB()->setError(N_BB*N_BB_ERR);
  pdfOnResDK.nBB()->setConstant();

  
  // Check if we initialized all parmeters
  pars = pdfOnResDK.getPdf()->getDependents(pdfOnResDK.parameters());
  iter = pars->createIterator();
  cout << "readOnResDKPar(): The following parameters might not have been initialized correctly:"<<endl;
  while ((arg = (RooAbsArg*)iter->Next())) {
    if (TString(arg->ClassName())=="RooRealVar") {      
      if (((RooRealVar*)arg)->getVal()==
          ((RooRealVar*)parsBeforeInit->find(arg->GetName()))->getVal()) {
        ((RooRealVar*)arg)->Print();
      }
    }
  }
  
  delete parsBeforeInit;
  delete pars;
  delete iter;
}


void setupPdfs() {
  cout << "--- START setupPdfs() ---" << endl;

  // Set verbosity first:
  pdfOnResDK.setVerbose(verbosePdf);

  // We don't use blinding for pdfOnResDK.
  // But the dalitzHolder have their own blinding for x and y.

  BdkPdfOnRes::Options opt = BdkPdfOnRes::none;
  
  // link good and bad signal asymmetries (for systematics only!)
  // opt = BdkPdfOnRes::linkGoodBadSigAsym;

  pdfOnResDK.init("pdfOnResDK", "PDF for On Res B->DK",
		  250000, *Hdtrkchge, 
		  *Deltae, *d0mass, *mes, *qprime, *dprime,
                  0, 0,                        // uses the owned analytical PDFs
                  //&nnHolder, &bknnHolder,    // uses the histPDFs
		  &dalitzHolderN, &dalitzHolderP,
                  opt);
		  
  // We never use these components
  pdfOnResDK.useMd(false);
  pdfOnResDK.useMes(false);
  
  // Make sure all other components are in
  pdfOnResDK.useDalitz(true);
  pdfOnResDK.useDE(true);
  pdfOnResDK.useNnComb(true);
  pdfOnResDK.useNnCont(true);

  // Initialize parameters
  readOnResDKPar();

  // Disable N2 by default
  pdfOnResDK.useNnComb(false);

  // Print out the list of efficiency functions used
  cout << "Using the following efficiency functions:"<<endl;
  printEff();
  
  dDalitz.init("dDalitz", "dDalitz", *m12, *m13);
  dkDalitz.init("dkDalitz", "dkDalitz", *m12, *m13);
  bgdDalitz.init("bgdDalitz", "bgdDalitz", *m12, *m13);;


  // setup Dalitz resoluton PDFs
  dalitzRes12 = BdkPdfDalitzRes("dalitzRes12","m12 Dalitz resolution",*m12);
  dalitzRes13 = BdkPdfDalitzRes("dalitzRes13","m13 Dalitz resolution",*m13);
  dalitzRes12->parameters().readFromFile("../BToDKTo3piK/params/dalitzRes.par");
  dalitzRes13->parameters().readFromFile("../BToDKTo3piK/params/dalitzRes.par");
  dalitzRes12->fixAll();
  dalitzRes13->fixAll();
  
  cout << "--- END setupPdfs() ---" << endl;
}


// 4 Gaussians with same mean value
BdkPdfAbsBase* BdkPdfDalitzRes(const char* name, const char* title, RooRealVar& x)
{
  BdkPdf2Gauss* pdfa = new BdkPdf2Gauss(TString(name)+"a",TString(title)+" A",x);
  BdkPdf2Gauss* pdfb = new BdkPdf2Gauss(TString(name)+"b",TString(title)+" B",x);
  //BdkPdfPolyn* lin = new BdkPdfPolyn(TString(name)+"_lin",TString(title)+" linear",x,1,1);

  TList list;
  list.Add(pdfa);
  list.Add(pdfb);
  //  list.Add(lin);

  //  BdkPdfSum* pdf = new BdkPdfSum(name,title,*pdfa,*pdfb);
  BdkPdfSum* pdf = new BdkPdfSum(name,title,list,BdkPdfSum::FULL);
  
  pdfa->linkResolution(0,0,0,pdfa->b1());
  pdfb->linkResolution(0,pdfa->b1(),0,pdfa->b1());

  return pdf;
}

#endif
