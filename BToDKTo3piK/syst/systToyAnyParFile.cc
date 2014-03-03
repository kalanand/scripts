// $Id: systToyAnyParFile.cc,v 1.1 2007/04/16 09:27:56 fwinkl Exp $
// Toy MC fits for fitting same dataset with different parameters

void systToyAnyParFile(const char* parFile, int nExp, int job,
                       const char* title="systToyAnyParFile", const char* resultDir = "./")
{  
  //  incVar(pdfOnResDK.parameters(), "dalitzHolderN.sigGoodD0.pdf.dalitzAmp.F0_1370_amp", 0.03);
  useBothFitVars();
  cout << "Reading file "<<parFile<<endl;
  pdfOnResDK.parameters().readFromFile(parFile);
  pdfOnResDK.parameters().Print("v");

  // Renormalize all amplitudes
  //  BdkDDalitzAmp::normalizeAll();
  
  // Renormalize B amplitudes
  dalitzHolderN.sigGoodD0Type().dalitzAmp()->calDDbarNorm(1e7);

  // Save new parameters (including normalization constants)
  RooArgSet* newParams = pdfOnResDK.parameters().snapshot(false);

  RooDataSet* fitResults = 0;

  TString name(title);
  name += "_";
  name += job;

  // Initialize seed with job number
  setRandomGenSeed(job);
  
  for (int i=0; i<nExp; i++) {
    // generate with original paramters
    // Make sure we generate data for all variables
    useBothFitVars();
    readOnResDKPar();

    gROOT->cd();
    data = pdfOnResDK.generate(pdfOnResDK.totalNumEvts(), true);

    // fit with original parameters
    fit(pdfOnResDK, *data, true);
    pdfOnResDK.yieldFitResult()->Print();
    pdfOnResDK.xyFitResult()->Print();

    RooArgSet* result = pdfOnResDK.fitResult();
    
    // fit with new parameters
    pdfOnResDK.parameters() = *newParams;
    fit(pdfOnResDK, *data, true);
    pdfOnResDK.yieldFitResult()->Print();
    pdfOnResDK.xyFitResult()->Print();
    
    RooArgSet* result2 = pdfOnResDK.fitResult();
    addSuffix(*result2,"_2");
    
    // merge the two results
    result->addOwned(*result2);
    // make sure errors are saved
    result->setAttribAll("StoreError",kTRUE);
    result->setAttribAll("StoreAsymError",kTRUE);
 
    // add results to dataset
    if (fitResults==0) fitResults = new RooDataSet(name,name,*result);
    fitResults->add(*result);
  }

  // Save all fit results
  TString filename(resultDir);
  filename += "/";
  filename += name;
  filename += ".root";

  TFile f(filename,"recreate");
  fitResults->Write();
  f.Close();
}


RooDataSet* getSystToyAnyParFileData(int nJobs,
                                     const char* title="systToyAnyParFile",
                                     const char* resultDir = "./")
{
  RooDataSet* data = 0;
  for (int i=1; i<=nJobs; i++) {
    TString name(title);
    name += "_";
    name += i;
    TString filename = TString(resultDir)+name+".root";

    cout << "Reading "<<filename<<endl;
    TFile f(filename);
    if (f.IsZombie()) {
      cout << "Cannot open "<<filename<<endl;
      continue;
    }
    RooDataSet* d = (RooDataSet*)f.Get(name);
    if (data==0) data = new RooDataSet(*d);
    else data->append(*d);
    f.Close();
  }
  cout << "Read "<<data->numEntries()<<" entries."<<endl;
  return data;
}
