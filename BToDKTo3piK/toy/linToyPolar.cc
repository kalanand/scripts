// $Id: linToyPolar.cc,v 1.2 2006/07/11 21:43:30 fwinkl Exp $
// Script for toy MC studies of fit linearity in rho and theta

RooDataSet* proto;
RooAbsPdf* pdf;

// BdkBatchMCStudy doesn't know how to deal with pointers
RooCategory kcharge(*Hdtrkchge,"*Hdtrkchge");


// setup the PDF (parameters are read from parFile if supplied)
void setupPdf(const char* parFile = 0) 
{ 
  proto = 0;   // don't need it for pdfOnResDK toys

  pdf = pdfOnResDK.getPdf();
  readOnResDKPar();

  if (parFile) pdfOnResDK.parameters().readFromFile(parFile);

  //  pdfOnResDK.parameters().
  //  readFromFile("../BToDKTo3piK/params/cross/plotNLL_RhoPhases2.par");

  pdfOnResDK.parameters().Print("v");

  pdfOnResDK.setNLLYieldsSystBit(0);
  pdfOnResDK.setNsigAsymFromXY();

  // Set correct names for BdkBatchMCStudy
  pdf->SetName("pdfOnResDK");
  m12->SetName("*m12");
  m13->SetName("*m13");
  Deltae->SetName("*Deltae");
  qprime->SetName("*qprime");
  dprime->SetName("*dprime");
}

// Create prototype dataset
RooDataSet* createProto(Int_t nEvents = 0)
{ 
  if (nEvents<=0) nEvents = pdfOnResDK.totalNumEvts();

  Double_t SumPos = 0;
  Double_t SumNeg = 0;
  
  for(Int_t i=0; i<BdkEvtTypes::NTYPES; i++) {
    SumPos += pdfOnResDK.numEvtPos(i)->getVal();
    SumNeg += pdfOnResDK.numEvtNeg(i)->getVal();
  }
  
  Int_t nEvts1 = (Int_t) (nEvents*SumPos/(SumPos+SumNeg));
  Int_t nEvts2 = (Int_t) (nEvents*SumNeg/(SumPos+SumNeg));

  // Create the prototype dataset with the RooCategory
  RooDataSet* proto = new RooDataSet("proto","",RooArgSet(*Hdtrkchge));

  cout << "Creating prototype dataset for "<<nEvts1<< " B+ events."<<endl;
  for (int i=0;i<nEvts1;i++) {    
    Hdtrkchge->setLabel("+");
    proto->add(*Hdtrkchge);
  }

  cout << "Creating prototype dataset for "<<nEvts2<< " B- events."<<endl;
  for (int i=0;i<nEvts2;i++) {    
    Hdtrkchge->setLabel("-");
    proto->add(*Hdtrkchge);
  }

  return proto;
}


// Create the setup script for the different toy MCs
// Example:
// batchScripts(-0.2,0.2,0.05,-0.2,0.2,0.05,1000,5000,100,"./awg/toyMC/linTest")
// in case you are using pdfOnResDK set nEvtPerSample = pdfOnResDK.totalNumEvts()
void batchScripts(Double_t rhoMin, Double_t rhoMax, Double_t rhoStep,
                  Double_t thetaMin, Double_t thetaMax, Double_t thetaStep,
                  Int_t fitsPerJob, Int_t nSamples, Int_t nEvtPerSample,
                  const char* dir)
{
  setupPdf();   // dummy
  BdkBatchMCStudy *toyMC = new BdkBatchMCStudy(pdfOnResDK, pdfOnResDK,
                                               RooArgSet(*m12,*m13,*Deltae,*qprime,*dprime),
                                               "e","me",proto,RooArgSet(kcharge));
  toyMC->setSubmitOption("-q xlong -C 0");

  TString dirFile = TString(dir)+"/dirs.txt";
  ofstream ofDirs;
  ofDirs.open(dirFile);
  
  for (Double_t rho=rhoMin; rho<=rhoMax; rho+=rhoStep) {
    for (Double_t theta=thetaMin; theta<=thetaMax; theta+=thetaStep) {
      
      // Create the directory
      TString srho, stheta;
      srho.Form("%.2f",rho);
      stheta.Form("%.2f",theta);
      TString path = TString(dir)+TString("/toy_")+srho+"_"+stheta;
      if (gSystem->mkdir(path)<0) {
        cout << "Cannot create "<<path<<endl;
        return;
      }
      // Write directory name to text file
      ofDirs << path <<endl;

      // Initialize PDF and create par file     
      dalitzHolderP.sigGoodD0Type().rho()->setVal(rho);
      dalitzHolderP.sigGoodD0Type().theta()->setVal(theta);
      dalitzHolderN.sigGoodD0Type().rho()->setVal(rho);
      dalitzHolderN.sigGoodD0Type().theta()->setVal(theta);

      dalitzHolderP.sigGoodD0Type().rho()->setError(0.1);
      dalitzHolderP.sigGoodD0Type().theta()->setError(1);
      dalitzHolderN.sigGoodD0Type().rho()->setError(0.1);
      dalitzHolderN.sigGoodD0Type().theta()->setError(1);

      // The rho- phase: 
      /*
      RooRealVar* phaseN = (RooRealVar*)pdfOnResDK.parameters().
        find("dalitzHolderN.sigGoodD0.pdf.dalitzAmp.Rho-_phase");
      RooRealVar* phaseP = (RooRealVar*)pdfOnResDK.parameters().
        find("dalitzHolderP.sigGoodD0.pdf.dalitzAmp.Rho-_phase");

      phaseN->setVal(theta);
      phaseP->setVal(theta);
      */

      TString parFile = path+"/pdf.par";
      pdfOnResDK.parameters().writeToFile(parFile);

      // Create setup script
      TString toySetup = path+"/toy-setup.cc";
      cout <<"Creating setup script "<<toySetup<<endl;
      ofstream of;
      of.open(toySetup,ios_base::out);
      of << "void toy_setup() {"<<endl
         << "gROOT->ProcessLine(\".x setup.cc\");"<<endl
         << "gROOT->ProcessLine(\".L linToyPolar.cc\");"<<endl
         << "setupPdf(\""<<parFile<<"\");"<<endl
         << "}" <<endl;
      of.close();

      toyMC->createScripts(0,path+"/fitMC",toySetup,fitsPerJob,nSamples,nEvtPerSample);
    }
  }
  ofDirs.close();
}
  

// Merge all toy MCs listed in dirFile
void mergeLinToy(const char* dirFile)
{
  ifstream ifdirs;
  ifdirs.open(dirFile);
  while (!ifdirs.eof()) {
    TString dir;
    ifdirs >> dir;   
    if (dir=="") continue;
    // Merge toy files
    cout << "Merging "<<dir<<endl;
    TString files(dir+"/fitMC.files");
    TString out(dir+TString("/fitMC.root"));
    mergeBatchMCStudy(files, out, true);
  }
  ifdirs.close();
}


// Merge all toy MCs listed in dirFile in batch queue
void mergeLinToyBatch(const char* dirFile, Int_t mergesPerJob = 10)
{
  ifstream ifdirs;
  ifdirs.open(dirFile);

  // Submit script
  ofstream sub;
  TString mergeSub(TString(dirFile)+"-merge-sub");
  sub.open(mergeSub);

  int n = 0;
  TString merge;
  ofstream mergeScript;    
  while (!ifdirs.eof()) {
    TString dir;
    ifdirs >> dir;   
    if (dir=="") continue;
    // Create merge script

    if ((n % mergesPerJob) == 0) {
      if (n>0) {
        mergeScript << "}" <<endl;
        mergeScript.close();
        sub << "bsub -q long -o "<<merge
            <<".log bbrroot -b -l -q setup.cc " << merge <<endl;
      }

      merge = dir+"/../merge-";
      merge += n/mergesPerJob;
      merge += ".cc";
      cout << "Creating merge script "<<merge<<endl;
      mergeScript.open(merge);
      mergeScript << "{" <<endl;
    }

    mergeScript << "mergeBatchMCStudy(\"" << dir <<"/fitMC.files\",\""
                << dir << "/fitMC.root\");" <<endl;
    n++;

  }  
  mergeScript.close();
  sub << "bsub -q long -o "<<merge
      <<".log bbrroot -b -l -q setup.cc " << merge <<endl;

  ifdirs.close();
  sub.close();
  cout << "Submit script written to "<<mergeSub<<endl;
}

// Plot and fit all toy MCs listed in dirFile
// Save results of 
void fitLinToy(const char* dirFile)              
{

  const Bool_t fitPull = kTRUE;
  const Bool_t fitErr = kFALSE;

  // Loop over list of directories
  ifstream ifdirs;
  ifdirs.open(dirFile);
  while (!ifdirs.eof()) {
    TString dir;
    ifdirs >> dir;   
    if (dir=="") continue;

    // Open file with toy MC results
    TFile* f = new TFile(dir+"/fitMC.root");
    RooDataSet* data = (RooDataSet*)f->Get("fitParData");
    
    TCanvas* can = new TCanvas("can","Signal toy MC",1200,600);
    can->Divide(4,2);

    RooRealVar rhoP(*pdfOnResDK.rhoPlus());
    RooRealVar rhoN(*pdfOnResDK.rhoMinus());
    RooRealVar thetaP(*pdfOnResDK.thetaPlus());
    RooRealVar thetaN(*pdfOnResDK.thetaMinus());
    rhoP.SetTitle("rho+");
    rhoN.SetTitle("rho-");
    thetaP.SetTitle("theta+");
    thetaN.SetTitle("theta-");

    RooArgSet set;
    // Plot rho+/-
    plotToy(data, rhoP, fitPull, fitErr, 0, 0.2, 40, -4, 4, 40, can, 1);    
    if (toyGaussPull) {
      set.add(*(RooAbsArg*)toyGaussPull->b()->Clone("rhoPpull.b"));
      set.add(*(RooAbsArg*)toyGaussPull->s()->Clone("rhoPpull.s"));
    }    
    plotToy(data, rhoN, fitPull, fitErr, 0, 0.2, 40, -4, 4, 40, can, 5);    
    if (toyGaussPull) {
      set.add(*(RooAbsArg*)toyGaussPull->b()->Clone("rhoNpull.b"));
      set.add(*(RooAbsArg*)toyGaussPull->s()->Clone("rhoNpull.s"));
    }    

    // Plot theta +/-
    plotToy(data, thetaP, fitPull, fitErr, 0, 60, 40, -4, 4, 40, can, 3);
    if (toyGaussPull) {
      set.add(*(RooAbsArg*)toyGaussPull->b()->Clone("thetaPpull.b"));
      set.add(*(RooAbsArg*)toyGaussPull->s()->Clone("thetaPpull.s"));      
    }
    plotToy(data, thetaN, fitPull, fitErr, 0, 60, 40, -4, 4, 40, can, 7);
    if (toyGaussPull) {
      set.add(*(RooAbsArg*)toyGaussPull->b()->Clone("thetaNpull.b"));
      set.add(*(RooAbsArg*)toyGaussPull->s()->Clone("thetaNpull.s"));      
    }


    // Get the peak of the error distributions
    RooRealVar rhoPerrMax("rhoPerrMax","",0); set.add(rhoPerrMax);
    RooRealVar rhoNerrMax("rhoNerrMax","",0); set.add(rhoNerrMax);
    RooRealVar thetaPerrMax("thetaPerrMax","",0); set.add(thetaPerrMax);
    RooRealVar thetaNerrMax("thetaNerrMax","",0); set.add(thetaNerrMax);

    TH1D herr("herr","",100,0,1);    
    data->tree().Project("herr",rhoP.GetName()+TString("err"));
    rhoPerrMax.setVal(herr.GetXaxis()->GetBinCenter(herr.GetMaximumBin()));
    data->tree().Project("herr",rhoN.GetName()+TString("err"));
    rhoNerrMax.setVal(herr.GetXaxis()->GetBinCenter(herr.GetMaximumBin()));
    data->tree().Project("herr",thetaP.GetName()+TString("err"));
    thetaPerrMax.setVal(herr.GetXaxis()->GetBinCenter(herr.GetMaximumBin()));
    data->tree().Project("herr",thetaN.GetName()+TString("err"));
    thetaNerrMax.setVal(herr.GetXaxis()->GetBinCenter(herr.GetMaximumBin()));
    

    // Get the median
    RooRealVar rhoPerrMedian("rhoPerrMedian","",0); set.add(rhoPerrMedian);
    RooRealVar rhoNerrMedian("rhoNerrMedian","",0); set.add(rhoNerrMedian);
    RooRealVar thetaPerrMedian("thetaPerrMedian","",0); set.add(thetaPerrMedian);
    RooRealVar thetaNerrMedian("thetaNerrMedian","",0); set.add(thetaNerrMedian);

    rhoPerrMedian.setVal(getMedian(data->tree(),rhoP.GetName()+TString("err")));
    rhoNerrMedian.setVal(getMedian(data->tree(),rhoN.GetName()+TString("err")));
    thetaPerrMedian.setVal(getMedian(data->tree(),thetaP.GetName()+TString("err")));
    thetaNerrMedian.setVal(getMedian(data->tree(),thetaN.GetName()+TString("err")));

    TString parFile = dir+"/fit.par";
    set.writeToFile(parFile);

    can->SaveAs(dir+"/fitMC.eps");
    delete can;
  }
  ifdirs.close();
}

void printLinToy(const char* dirFile, const char* var)
{
  const Int_t digits = 1;
  
  // Loop over list of directories
  ifstream ifdirs;
  ifdirs.open(dirFile);

  Double_t x;
  ostringstream ss;
  ostringstream error;
  ostringstream errorMax;
  ostringstream header;
  bool headerDone = false;
  bool first = true;
  int cols = 1;
  header << "x/y &";
  while (!ifdirs.eof()) {
    TString dir;
    ifdirs >> dir;   
    if (dir=="") continue;

    RooRealVar b(TString(var)+"pull.b","",1);
    RooRealVar s(TString(var)+"pull.s","",1);
    RooRealVar errMax(TString(var)+"errMax","",1);
    RooRealVar errMedian(TString(var)+"errMedian","",1);

    RooArgSet set(b,s,errMax,errMedian);
    set.readFromFile(dir+"/fit.par");
    dalitzHolderP.sigGoodD0Type().parameters().readFromFile(dir+"/pdf.par");
    if (first) {
       x = dalitzHolderP.sigGoodD0Type().x()->getVal();
       ss << "\\multirow{3}{10mm}{"<< x <<"}";
       ss << "& $\\mu$";
       first = false;
    }   
    if (dalitzHolderP.sigGoodD0Type().x()->getVal()!=x) {
      x = dalitzHolderP.sigGoodD0Type().x()->getVal();
      headerDone = true;
      ss << "\\"<<"\\" << endl
         << "& $\\sigma$"
         << error.str() << "\\"<<"\\" << endl
         << "& $e$"
         << errorMax.str() << "\\"<<"\\" <<"\\hline"<< endl
         << "\\multirow{3}{10mm}{"<< x <<"}"
         << "& $\\mu$";
      error.str("");
      errorMax.str("");
    }
    if (!headerDone) {
      header << " & "<< dalitzHolderP.sigGoodD0Type().y()->getVal();
      cols++;
    }
    TString str;
    str.Form("$%.2f\\pm %.2f$",b.getVal(),b.getError());
    ss <<" & "<< str;
    str.Form("$%.2f\\pm %.2f$",s.getVal(),s.getError());
    error << " & "<< str;
    str.Form("$%.2f$",errMax.getVal());
    errorMax << " & " << str;
    str.Form("$%.2f$",errMedian.getVal());
    errorMax << " / " << str;
  }
  ifdirs.close();
  ss << "\\"<<"\\" << endl << "& $\\sigma$" << error.str()<<"\\"<<"\\" << endl
     << "& $e$" << errorMax.str();
  cout << "\\begin{tabular}{";
  for (int i=0;i<=cols;i++) {
    cout <<"c";
    if (i<2) cout <<"|";
  }
  cout <<"}"<<endl;
  cout <<"\\hline\\hline"<<endl;
  cout << header.str()<<"\\"<<"\\"<<endl;
  cout <<"\\hline"<<endl;
  cout << ss.str()<<"\\"<<"\\"<<endl;
  cout <<"\\hline\\hline"<<endl;
  cout <<"\\end{tabular}"<<endl;
}
