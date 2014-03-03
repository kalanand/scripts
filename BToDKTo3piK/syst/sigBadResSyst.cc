// $Id: sigBadResSyst.cc,v 1.3 2006/09/25 23:31:13 fwinkl Exp $
// Script to evaluate the bad signale shape systematic error
// 1) Generate good D signal
// 2) Smear it by the reconstruction resolution obtained from flat signal MC
// 3) Obtain a new hist PDF for the data fit
//
// createDKBadDHist() creates the histogram
// sigBadResSyst() evaluates the systematics assigned with that

// Can also be used to evaluate the systematics due to a wrong kaon in
// the bad signal shape. See where the TFile is openend in the function below.

void sigBadResSyst()
{
  readCut = cutSigReg;
  data = read(dataTree);

  // reference fit
  cout << "-----------------------------------------------------"<<endl;
  cout << "   Fitting with original DKBadD shape"<<endl;
  cout << "-----------------------------------------------------"<<endl;
  useBothFitVars();
  readOnResDKPar();

  fit(pdfOnResDK, *data, true);  
  RooArgSet* ref = pdfOnResDK.fitResult();
  ref->setName("ref");


  // alternative fit
  cout << "-----------------------------------------------------"<<endl;
  cout << "   Fitting with MC based DKBadD shape"<<endl;
  cout << "-----------------------------------------------------"<<endl;
  useBothFitVars();
  readOnResDKPar();

  // Switch to MC based DKBadD histogram
  //  TFile f("../BToDKTo3piK/params/syst/hist_dkbadd_toy.root");

  // Switch to DKBadD histogram with good K's / bad D's only
  //  TFile f("../BToDKTo3piK/params/syst/hist_dkbadd_goodk.root");

  // Switch to DKBadD histogram with bad K events replaced by signal
  TFile f("../BToDKTo3piK/params/syst/hist_dkbadd_replacek.root");

  TH2* h = (TH2*)f.Get("hist_dkbadd");
  h->SetDirectory(gROOT);
  f.Close();
  dalitzHolderN.setSigBadD0Hist(*h);
  dalitzHolderP.setSigBadD0Hist(*h);

  fit(pdfOnResDK, *data, true);  
  RooArgSet* alt = pdfOnResDK.fitResult();
  alt->setName("alt");

  ref->Print("v");
  alt->Print("v");

  removeAllRanges(*ref);
  sub(*ref, *alt);
  ref->Print("v");
}


// Requires getDalitzParams.cc
// Need to change eff -> effDKBadD in setupPdfHolder.cc
void createDKBadDHist() 
{
  const int nEvents = 15000;

  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/params/syst/effDKBadD-norm.par");

  readCut = cutSigReg+cutDKBadD;
  RooDataSet* resData = read(sigFlatTree);
  
  cout << "Generating "<<nEvents/2<<" DKGoodD B- events..."<<endl;
  Hdtrkchge->setLabel("-");
  data = dalitzHolderN.sigGoodD0Type().generate(nEvents/2);
  data->addColumn(*Hdtrkchge);

  cout << "Generating "<<nEvents/2<<" DKGoodD B+ events..."<<endl;
  Hdtrkchge->setLabel("+");
  RooDataSet* dataP = dalitzHolderP.sigGoodD0Type().generate(nEvents/2);
  dataP->addColumn(*Hdtrkchge);

  data->append(*dataP);

  cout << "Smearing events..."<<endl;
  RooDataSet* smear = smearData(data, resData);
  cout << smear->numEntries() << " events in Dalitz plot after smearing."<<endl;

  // The unsquared variables are needed by getDalitzParams
  RooFormulaVar sqrt12("d0pppupmass","d0pppupmass","sqrt(@0)",RooArgList(*m12));
  RooFormulaVar sqrt13("d0ppmupmass","d0ppmupmass","sqrt(@0)",RooArgList(*m13));

  smear->addColumn(sqrt12);
  smear->addColumn(sqrt13);

  TCanvas* can = new TCanvas("smear","smear",800,400);
  can->Divide(2,1);
  can->cd(1);
  data->tree().Draw("m13:m12");
  can->cd(2);
  smear->tree().Draw("m13:m12");  


  TCut cut = "";
  cut.SetName("smear");  // This name is used as a filename in the following function

  getDalitzParamsDKBadD(&smear->tree(),cut);
} 



RooDataSet* smearData(RooDataSet* data, RooDataSet* resData)
{
  RooDataSet* d = (RooDataSet*)data->emptyClone();
  for (int i=0; i<data->numEntries(); i++) {
    RooArgSet* set = (RooArgSet*)data->get(i);
    Double_t r12 = set->getRealValue("m12",-99);
    Double_t r13 = set->getRealValue("m13",-99);
    smearValue(resData, r12, r13);
    // Only take smeared value if it is inside the Dalitz
    if (dalitzCfg->inDalitz(r12,r13)) {
      set->setRealValue("m12",r12);
      set->setRealValue("m13",r13);
      d->add(*set);
    }
  }
  return d;
}



// Smear m12 and m13 by the resolution found in data
void smearValue(RooDataSet* resData, Double_t& m12, Double_t& m13)
{
  // size of bins to randomly pick resolution from
  const Double_t binSize = 0.3;
  
  Int_t N = resData->numEntries();

  while (1) {
    // Get random event
    const RooArgSet* set = resData->get(gRandom->Integer(N));    
    Double_t r12 = set->getRealValue("m12",-99);
    Double_t r13 = set->getRealValue("m13",-99);

    // Check if this event is within (binSize x binSize) of (m12,m13)
    if ((fabs(r12-m12)<binSize) && (fabs(r13-m13)<binSize)) {
      // Get generated Dalitz values
      Double_t true12 = pow(set->getRealValue("d0pppmcmass",-99),2);
      Double_t true13 = pow(set->getRealValue("d0ppmmcmass",-99),2);
      // Smear by resolution
      m12 += (r12-true12);
      m13 += (r13-true13);
      break;
    }
  }
}
