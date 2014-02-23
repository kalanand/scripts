// $Id: printMcStatSyst.cc,v 1.2 2006/08/04 22:27:54 fwinkl Exp $
// Script to print the results from mcStatSyst.cc

void printMcStatSyst(TString dir = "dataFits/mcStatFits/",
		     Bool_t verbose = kFALSE)
{  
  // read all fit results
  TList* list = readAllFiles(dir);

  // Variables to evaluate systematics for
  RooArgList systList;
  systList.add(*pdfOnResDK.typeAsym(BdkEvtTypes::SIG_GOOD_D));
  systList.add(*pdfOnResDK.sigGoodD0NumEvts());
  systList.add(pdfOnResDK.cpParams());

  // Get reference fit
  RooAbsCollection* ref = ((RooArgSet*)list->FindObject("ref"))->selectCommon(systList);
  cout << "-----------------------------------------------------"<<endl;
  cout << "Reference fit"<<endl;
  if (ref) ref->Print("v");
  else return;
  cout << "-----------------------------------------------------"<<endl;

  removeAllRanges(*ref);

  RooAbsCollection* sum = ref->snapshot(false);
  RooAbsCollection* mean = ref->snapshot(false);
  zero(*sum);
  zero(*mean);

  TIterator* iter = list->MakeIterator();
  RooArgSet* set;
  int N = 0;
  TMatrixDSym V(systList.getSize());
  while (set = (RooArgSet*)iter->Next()) {

    if (TString(set->GetName())=="ref") continue;
    RooAbsCollection* diff = ref->snapshot(false);
    diff->setName(set->GetName());

    add(*mean, *set->selectCommon(systList));

    // Add differences to reference squared
    N++;
    sub(*diff, *set->selectCommon(systList));   // diff is the reference
    removeAllErrors(*diff);
    if (verbose) diff->Print("v");
    
    for (int n=0; n<V.GetNrows(); n++) {    
      for (int m=0; m<=n; m++) {
        V(n,m) += ((RooArgSet*)diff)->getRealValue(systList[n].GetName()) *
                  ((RooArgSet*)diff)->getRealValue(systList[m].GetName());
      }
    }
    
    sqr(*diff);
    add(*sum, *diff);
  
    delete diff;
  }
  
  for (int n=0; n<V.GetNrows(); n++) {    
    for (int m=0; m<=n; m++) {
      V(n,m) /= 2.0;
      V(m,n) = V(n,m);
    }
  }

  RooRealVar invN("invN","1/N",1.0/N);
  multiply(*mean,invN);
  cout << "-----------------------------------------------------"<<endl;
  cout << "Average over "<< N << " fit results"<<endl;
  mean->Print("v");
  cout << "-----------------------------------------------------"<<endl;

  RooRealVar half("half","0.5",0.5);
  multiply(*sum,half);
  sqrt(*sum);
  cout << "-----------------------------------------------------"<<endl;
  cout << "Systematics from "<< N << " fit results"<<endl;
  removeAllErrors(*sum);
  sum->Print("v");
  cout << "-----------------------------------------------------"<<endl;

  cout << "Covariance matrix:"<<endl;
  V.Print();

  delete iter;
}


void printFitDiffs(TString dir = "dataFits/mcStatFits/",
                   const RooRealVar* var)
{
  if (var==0) return;
  
  TList* list = readAllFiles(dir);
  RooArgSet* refSet = (RooArgSet*)list->FindObject("ref");
  RooRealVar* ref = (RooRealVar*)refSet->find(var->GetName());
  
  RooArgList* sortedDiff = sortFits(list,*var,*ref);
  sortedDiff->Print("v");
  
  delete list;
}

void printAllFitDiffs(TString dir = "dataFits/mcStatFits/")
{
  TList* list = readAllFiles(dir);
  RooArgSet* refSet = (RooArgSet*)list->FindObject("ref");

  TIterator* iter = refSet->createIterator();
  RooRealVar* var;
  while (var = (RooRealVar*)iter->Next()) {
    cout << "------------------------------------------------------"<<endl;
    cout << "Differences for "<<var->GetName()<<endl;
    cout << "------------------------------------------------------"<<endl;
    RooRealVar* ref = (RooRealVar*)refSet->find(var->GetName());
    sortFits(list,*var,*ref)->Print("v");
  }
}

// Return a list with all RooArgSets from fit files in dir
TList* readAllFiles(TString dir)
{
  TString fileBase = dir+"/mcStatSyst-";

  TList* setList = new TList();
  
  for (int t=0; t<BdkEvtTypes::NTYPES; t++) {
    TString filename = fileBase;
    filename += t;
    filename += ".root";

    TFile f(filename);
    if (f.IsZombie()) {
      cout << "Cannot open "<<filename<<endl;
      continue;
    }

    TList* keys = f.GetListOfKeys();
    TIterator* iter = keys->MakeIterator();
    TKey* k;
    int n = 0;
    while (k = (TKey*)iter->Next()) {
      RooArgSet* set = (RooArgSet*)k->ReadObj();
      //      if (TString(set.GetName())=="ref") set->Print("v");
      setList->Add(set);
      n++;
    }
    delete iter;
    f.Close();
    cout << n << " fit results read from "<<filename<<endl;
  }
  cout << setList->GetEntries() << " fit results in total."<<endl;
  return setList;
}


// Return sorted list of fits with descending absolute difference to ref
RooArgList* sortFits(TList* fitList, 
                     const RooRealVar& var, const RooRealVar& ref,
                     Bool_t absDiff = kTRUE)
{
  RooArgList* list = new RooArgList();

  // Make a list of all fit results for var
  TIterator* iter = fitList->MakeIterator();
  RooArgSet* set;
  while (set = (RooArgSet*)iter->Next()) {
    //if (TString(set->GetName())=="ref") continue;
    RooRealVar* rFit = (RooRealVar*)set->find(var.GetName());
    if (rFit==0) {
      cout << "Could not find "<<var.GetName()<<" in "<<set->GetName() << endl;
      continue;
    }

    double diff = ref.getVal() - rFit->getVal();
    if (absDiff) diff = fabs(diff);

    BdkRooRealVar* r = new BdkRooRealVar(set->GetName(),set->GetTitle(),diff);
    list->addOwned(*r);    
  }

  // BdkRooRealVar knows how to compare itself numerically to other RooRealVars.
  // That's why the following sort() sorts numerically instead of alpha-numerically.
  list->sort(true);   // true = largest value first

  return list;
}
