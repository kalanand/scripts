// $Id: linToyTest.cc,v 1.3 2006/05/19 17:34:42 fwinkl Exp $
// Script for toy MC studies of fit linearity in x and y

RooDataSet* proto;
RooAbsPdf* pdf;
BdkPdfAbsBase& bdkPdf = pdfOnResDK;   // dalitzHolderP.sigGoodD0Type();

// BdkBatchMCStudy doesn't know how to deal with pointers
RooCategory kcharge(*Hdtrkchge,"*Hdtrkchge");


// setup the PDF (parameters are read from parFile if supplied)
void setupPdf(const char* parFile = 0) 
{
  if (parFile) bdkPdf.parameters().readFromFile(parFile);
  pdf = bdkPdf.getPdf();
  bdkPdf.parameters().Print("v");

  proto = createProto();

  // Set correct names for BdkBatchMCStudy
  pdf->SetName("*pdf");
  m12->SetName("*m12");
  m13->SetName("*m13");
  Deltae->SetName("*Deltae");
  qprime->SetName("*qprime");
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
void batchScripts(Double_t xmin, Double_t xmax, Double_t xstep,
                  Double_t ymin, Double_t ymax, Double_t ystep,
                  Int_t fitsPerJob, Int_t nSamples, Int_t nEvtPerSample,
                  const char* dir)
{
  setupPdf();   // dummy
  BdkBatchMCStudy *toyMC = new BdkBatchMCStudy(*pdf,*pdf,
                                               //RooArgSet(*m12,*m13),"","m");
                                               RooArgSet(*m12,*m13,*Deltae,*qprime),"","me",proto,RooArgSet(kcharge));

  TString dirFile = TString(dir)+"/dirs.txt";
  ofstream ofDirs;
  ofDirs.open(dirFile);
  
  for (Double_t x=xmin; x<=xmax; x+=xstep) {
    for (Double_t y=ymin; y<=ymax; y+=ystep) {
      
      // Create the directory
      TString sx, sy;
      sx.Form("%.2f",x);
      sy.Form("%.2f",y);
      TString path = TString(dir)+TString("/toy_")+sx+"_"+sy;
      if (gSystem->mkdir(path)<0) {
        cout << "Cannot create "<<path<<endl;
        return;
      }
      // Write directory name to text file
      ofDirs << path <<endl;

      // Initialize PDF and create par file      
      dalitzHolderP.sigGoodD0Type().x()->setVal(x);
      dalitzHolderP.sigGoodD0Type().y()->setVal(y);
      dalitzHolderN.sigGoodD0Type().x()->setVal(x);
      dalitzHolderN.sigGoodD0Type().y()->setVal(y);

      dalitzHolderP.sigGoodD0Type().x()->setError(0.01);
      dalitzHolderP.sigGoodD0Type().y()->setError(0.01);
      dalitzHolderN.sigGoodD0Type().x()->setError(0.01);
      dalitzHolderN.sigGoodD0Type().y()->setError(0.01);

      TString parFile = path+"/pdf.par";
      bdkPdf.parameters().writeToFile(parFile);

      // Create setup script
      TString toySetup = path+"/toy-setup.cc";
      cout <<"Creating setup script "<<toySetup<<endl;
      ofstream of;
      of.open(toySetup,ios_base::out);
      of << "void toy_setup() {"<<endl
         << "gROOT->ProcessLine(\".x setup.cc\");"<<endl
         << "gROOT->ProcessLine(\".L linToyTest.cc\");"<<endl
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
    mergeBatchMCStudy(files, out);
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
    
    TCanvas* can = new TCanvas("can","Signal toy MC",500,500);
    can->Divide(2,2);

    RooRealVar Bpx(*dalitzHolderP.sigGoodD0Type().x());
    RooRealVar Bpy(*dalitzHolderP.sigGoodD0Type().y());
    Bpx.SetTitle("B+ x");
    Bpy.SetTitle("B+ y");

    RooArgSet set;
    // Plot x
    plotToy(data, Bpx, fitPull, fitErr, 0, 1, 40, -4, 4, 40, can, 1);    
    if (toyGaussPull) {
      set.add(*(RooAbsArg*)toyGaussPull->b()->Clone("xpull.b"));
      set.add(*(RooAbsArg*)toyGaussPull->s()->Clone("xpull.s"));
    }    
    // Plot y
    plotToy(data, Bpy, fitPull, fitErr, 0, 1, 40, -4, 4, 40, can, 3);
    if (toyGaussPull) {
      set.add(*(RooAbsArg*)toyGaussPull->b()->Clone("ypull.b"));
      set.add(*(RooAbsArg*)toyGaussPull->s()->Clone("ypull.s"));      
    }

    // Get the peak of the error distributions
    RooRealVar xerrMax("xerrMax","",0);
    RooRealVar yerrMax("yerrMax","",0);

    TH1D herr("herr","",100,0,1);    
    data->tree().Project("herr",Bpx.GetName()+TString("err"));
    xerrMax.setVal(herr.GetXaxis()->GetBinCenter(herr.GetMaximumBin()));
    data->tree().Project("herr",Bpy.GetName()+TString("err"));
    yerrMax.setVal(herr.GetXaxis()->GetBinCenter(herr.GetMaximumBin()));
    set.add(xerrMax);
    set.add(yerrMax);

    // Get the median
    RooRealVar xerrMedian("xerrMedian","",0);
    RooRealVar yerrMedian("yerrMedian","",0);
    xerrMedian.setVal(getMedian(data->tree(),Bpx.GetName()+TString("err")));
    yerrMedian.setVal(getMedian(data->tree(),Bpy.GetName()+TString("err")));
    set.add(xerrMedian);
    set.add(yerrMedian);

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
