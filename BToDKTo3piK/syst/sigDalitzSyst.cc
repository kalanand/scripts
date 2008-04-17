// $Id: sigDalitzSyst.cc,v 1.3 2006/06/19 17:58:29 fwinkl Exp $
// Obtain the systematic error due to uncertainties in the signal Dalitz model
// Run this in batch with:  bbrroot -q -b setup.cc ".x sigDalitzSyst.cc"

void sigDalitzSyst(TTree* tree = dataTree,
                   TCut addCut = "")
{
  TString path = "../BToDKTo3piK/params/syst/";

  // List of signal Dalitz models
  const int nParFiles = 9;
  TString parFile[nParFiles] = {"dstar_PPP_Group0.Bdk.par",
				"dstar_PPP_Group1.Bdk.par",
				"dstar_PPP_Group2.Bdk.par",
				"dstar_PPP_Group3.Bdk.par",
				"dstar_PPP_Group4.Bdk.par",
				"dstar_PPP_Group5.Bdk.par",
				"dstar_PPP_Group6.Bdk.par",
				"dstar_PPP_Group7.Bdk.par",
				"dstar_PPP_removeOmegaF2NRPW.Bdk.par"};

  RooFitResult* fitResults[nParFiles];

  // Read data
  readCut = cutSigReg+addCut;
  data = read(tree);

  // Do yield fit
  fitOption = "mer";
  useYieldFitVars();

  !!!! need to think about this one again !!!!
  RooFitResult* yieldFit = fit(pdfOnResDK,*data);

  // Do all x/y fits
  for (int i=0; i<nParFiles; i++) {
    cout << "==============================================================="<<endl;
    cout << "        x/y fit number "<<i<<endl;
    cout << "==============================================================="<<endl;

    useXyFitVars();

    TString filename = path+parFile[i];

    const BdkDDalitzAmp& amp = *dalitzHolderN.sigGoodD0Type().dalitzAmp();
    cout << "Reading resonance parameters from "<<filename<<":"<<endl;
    pdfOnResDK.parameters().readFromFile(filename);

    // Print amplitudes and phases for all resonances
    for (int r=0; r<amp.nComps(); r++) {
      cout << amp.nameRes(r) <<": "
           << "amp = " << amp.ampRes(r)->getVal()
           << ", phase = " << amp.phaseRes(r)->getVal()
           << ", mass = " << amp.massRes(r)->getVal()
           << ", width = " << amp.gammaRes(r)->getVal()
           << endl;
    }

    fit(pdfOnResDK, *data, *yieldFit);
    fitResults[i] = pdfOnResDK.xyFitResult();    
  }

  // Print and save all fit results
  TFile f("sigDalitzSyst.root","recreate");
  for (int i=0; i<nParFiles; i++) {    
    cout << "Fit result for fit number "<<i<<endl;
    if (fitResults[i]) {
      fitResults[i]->Print();
      
      TString name = "sigDalitzSyst_";
      name += i;
      RooArgSet set(fitResults[i]->floatParsFinal());
      set.writeToFile(name+".par"); 

      fitResults[i]->Write(name);
    }
  }    
  f.Close();

  printSigDalitzSyst("sigDalitzSyst.root");
}

// Print and calculate systematics compared to fit with number "ref"
void printSigDalitzSyst(const char* file, 
                        Bool_t printLatex = false,
                        Int_t ref = 0)
{
  TFile f(file);
  TList fitResults;
  int i = 0;
  while (1) {
    TString name = "sigDalitzSyst_";
    name += i;
    RooFitResult* result = (RooFitResult*)f.Get(name);
    if (result==0) break;
    fitResults.Add(result);
    i++;
  }
  f.Close();

  if (ref >= fitResults.GetSize()) {
    cout << "Invalid reference number."<<endl;
    return;
  }

  // Print fit results sorted by fit parameter
  TIterator* iter = ((RooFitResult*)fitResults.At(0))->floatParsFinal().createIterator();
  RooAbsArg* arg;
  while (arg = (RooAbsArg*)iter->Next()) {
    cout << endl;
    for (int i=0; i<fitResults.GetSize(); i++) {
      RooRealVar* r = (RooRealVar*)((RooFitResult*)fitResults.At(i))->floatParsFinal().find(arg->GetName());
      if (printLatex) {
        TString s;
        s.Form("%s $%.3f\\pm%.3f$",r->GetName(),r->getVal(),r->getError());
        cout << s << endl;
      }
      else r->Print();
    }
  }
  delete iter;

  // Calculate systematic errors
  RooArgList* sum = (RooArgList*)((RooFitResult*)fitResults.At(ref))->floatParsFinal().snapshot(false);
  zero(*sum);
  for (int i=0; i<fitResults.GetSize(); i++) {
    if (i==ref) continue;   // skip reference
    RooArgSet* diff = ((RooFitResult*)fitResults.At(ref))->floatParsFinal().snapshot(false);

    // Add squared differences to reference
    sub(*diff,((RooFitResult*)fitResults.At(i))->floatParsFinal());
    sqr(*diff);
    add(*sum, *diff);

    delete diff;
  }

  sqrt(*sum);
  cout <<endl<< "Differences to fit "<<ref<<" added in quadrature:"<<endl;
  sum->Print("v");
  delete sum;
}
