// $Id: fixedFracSyst.cc,v 1.6 2006/07/03 17:36:06 fwinkl Exp $
// Obtain the systematic error due to the fixed fit fractions and other constants.
// Run this in batch with:  bbrroot -q -b setup.cc ".x fixedFracSyst.cc"

void fixedFracSyst(TTree* tree = dataTree,
                   TCut addCut = "")
{
  TList fitResults;

  useBothFitVars();
  readOnResDKPar();

  // list of variables to vary in fit by their error
  RooArgSet vars("vars");
  vars.add(*pdfOnResDK.qqGoodD0Frac());
  vars.add(*pdfOnResDK.BBGoodD0Frac());
  vars.add(*pdfOnResDK.DpiBadD0Frac());
  vars.add(*pdfOnResDK.sigBadD0Frac());

  // set 50% errors on all of the above
  setErrors(vars,0.5);

  // now add some more variables
  pdfOnResDK.DKXFrac()->setError(pdfOnResDK.DKXFrac()->getVal()*0.25);  // 25% on DKX
  vars.add(*pdfOnResDK.DKXFrac());

  vars.add(*pdfOnResDK.absoluteEff());
  vars.add(*pdfOnResDK.nBB());
  vars.add(*pdfOnResDK.brBtoDK());
  vars.add(*pdfOnResDK.brBtoDK());
  vars.add(*pdfOnResDK.brDto3piOverK2pi());


  // save true values
  RooArgSet* trueValues = (RooArgSet*)vars.snapshot(false);
  
  // Read data
  readCut = cutSigReg+addCut;
  data = read(tree);

  // reference fit
  cout << "-----------------------------------------------------"<<endl;
  cout << "   Fitting with original parameters"<<endl;
  cout << "-----------------------------------------------------"<<endl;
  vars.Print("v");

  fit(pdfOnResDK, *data, true);
  RooArgSet* result = pdfOnResDK.fitResult();
  result->setName("ref");
  fitResults.Add(result);

  // 0 = (value+=error), 1 = (value-=error)
  for (int i=0; i<2; i++) {
    TIterator* iter = vars.createIterator();
    RooRealVar* r;
    while (r = (RooRealVar*)iter->Next()) {

      useBothFitVars();
      readOnResDKPar();
      // reset all variables and errors
      vars = *trueValues;  

      // change value by +- error
      if (i==0) r->setVal(r->getVal()+r->getError());
      else r->setVal(r->getVal()-r->getError());

      // fit and save fit result
      cout << "-----------------------------------------------------"<<endl;
      cout << "   Fitting with varied "<<r->GetName() << endl;
      cout << "-----------------------------------------------------"<<endl;
      vars.Print("v");

      fit(pdfOnResDK, *data, true);
      RooArgSet* result = pdfOnResDK.fitResult();
      TString name = r->GetName();
      name += i;
      result->setName(name);
      fitResults.Add(result);
    }
  }

  // save fit results
  TFile f("fixedFracSyst.root","recreate");
  fitResults.Write();
  vars.Write();
  f.Close();

  printFixedFracSyst("fixedFracSyst.root");
}


// set the error on vars to fraction of the value
void setErrors(const RooArgSet& vars, Double_t fraction)
{
  TIterator* iter = vars.createIterator();

  RooRealVar* r;
  while (r = (RooRealVar*)iter->Next()) {
    r->setError(r->getVal()*fraction);
  }
  delete iter;
}


void printFixedFracSyst(const char* file, 
                        Bool_t printLatex = kFALSE,
                        const char* fmt = "%.4f",
                        Int_t nSyst = 999)
{
  // Variables to evaluate systematics for
  RooArgSet syst;
  syst.add(*pdfOnResDK.sigGoodD0NumEvts());
  syst.add(*pdfOnResDK.typeAsym(BdkEvtTypes::SIG_GOOD_D));
  syst.add(pdfOnResDK.cpParams());

  TFile f(file);
  RooArgSet* vars = (RooArgSet*)f.Get("vars");

  // Get reference fit
  RooArgSet* ref = (RooArgSet*)f.Get("ref");
  removeAllRanges(*ref);

  RooAbsCollection* total = ref->snapshot(false);
  total->setName("total");
  zero(*total);

  TIterator* iter = vars->createIterator();
  RooRealVar* r;
  int n = 0;
  // Loop over all variables in set
  while (r = (RooRealVar*)iter->Next()) {
    if (++n > nSyst) continue;

    RooArgSet* result0 = (RooArgSet*)f.Get(r->GetName()+TString("0"));
    RooArgSet* result1 = (RooArgSet*)f.Get(r->GetName()+TString("1"));

    RooAbsCollection* diff0 = ref->snapshot(false);
    RooAbsCollection* diff1 = ref->snapshot(false);

    sub(*diff0, *result0); sqr(*diff0);
    sub(*diff1, *result1); sqr(*diff1);

    RooRealVar half("half","",0.5);
    add(*diff0, *diff1); multiply(*diff0, half); sqrt(*diff0);
    
    /*
    // calculate relative error in %
    divide(*diff0, *ref);
    RooRealVar hundred("100","",100);
    multiply(*diff0, hundred); sqr(*diff0); sqrt(*diff0);
    */

    RooAbsCollection* printSet = diff0->selectCommon(syst);
    if (printLatex) {
      cout << setw(30) << r->GetName() << " & " << printOneLine(*printSet,fmt)
           << " \\\\"<<endl;
      //      printSet->printLatex(Columns(syst.getSize()),Format("e"));
    }
    else {
      cout << "Systematic error due to "<<r->GetName()<<" variation:"<<endl;
      printSet->Print("v");
    }

    // Add for total systematic error
    sqr(*diff0);
    add(*total, *diff0);

    delete printSet;
    delete diff0;
    delete diff1;
  }
  sqrt(*total);
  if (printLatex) cout << setw(30) << "total" << " & "
                       << printOneLine(*total->selectCommon(syst),fmt)
                       << endl;
  else total->selectCommon(syst)->Print("v");

  delete iter;
  delete total;
  f.Close();
}


