// $Id: fullToy.cc,v 1.25 2007/03/02 09:51:10 fwinkl Exp $
// Full toy MC study

RooDataSet* proto;
RooAbsPdf* pdf;

// BdkBatchMCStudy doesn't know how to deal with pointers
RooCategory kcharge(*Hdtrkchge,"*Hdtrkchge");

// setup the PDF
void setupPdf() {

  // proto = createProto();
  proto = 0;   // don't need it for pdfOnResDK toys

  pdf = pdfOnResDK.getPdf();
  readOnResDKPar();

  //  pdfOnResDK.parameters().readFromFile("numEvts.par");
  //  pdfOnResDK.parameters().readFromFile("numEvts10.par");
  //pdfOnResDK.parameters().readFromFile("numEvts-noBkg.par");
  //pdfOnResDK.parameters().readFromFile("numEvts-noBkg10.par");
  //  pdfOnResDK.parameters().readFromFile("numEvts-noBkg100.par");
  //pdfOnResDK.parameters().readFromFile("numEvts-noBkg-XYp.par");

  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/results/fitData.par");

  //  pdfOnResDK.setNLLYieldsSystBit(0);
  pdfOnResDK.setNsigAsymFromXY();

  pdfOnResDK.parameters().Print("v");

  // Set correct names for BdkBatchMCStudy
  //  pdf->SetName("*pdf");
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

// Create batch scripts
// Example: batchScripts(50,0,5000,"awg/toyMC/fullToy-x")
void batchScripts(int fitsPerJob, const char* toyFile = "fullToy.root",
                  int nSamples = 0,
                  const char* toyDir = "awg/toy/fullToy-x")
{ 
  setupPdf();

  // Copy this script to toy directory
  TString src = "../BToDKTo3piK/toy/fullToy.cc";
  TString dest = toyDir+TString("/fullToy.cc");
  cout << "Copying " << src << " --> " << dest << endl;
  gSystem->CopyFile(src,dest,true);

  // Create batch setup script
  TString batchSetup = toyDir+TString("/fullToy-batch.cc");
  ofstream of;
  of.open(batchSetup);
  of << "void fullToy_batch() {" << endl
     << "gROOT->ProcessLine(\".x setup.cc\");" << endl
     << "gROOT->ProcessLine(\".L "<< toyDir <<"/fullToy.cc" << "\");" << endl
     << "setupPdf();}" << endl;
  of.close();

  // one-stage fit
  //  BdkBatchMCStudy *toyMC = new BdkBatchMCStudy(*pdf,*pdf,RooArgSet(*m12,*m13,*Deltae,*qprime),
  //                                               "","me",proto,RooArgSet(kcharge));

  // two-stage fit
  BdkBatchMCStudy *toyMC = new BdkBatchMCStudy(pdfOnResDK, pdfOnResDK,
                                               RooArgSet(*m12,*m13,*Deltae,*qprime,*dprime),
                                               "e","me",proto,RooArgSet(kcharge));
  //					       "e","e",proto,RooArgSet(kcharge));
  toyMC->setSubmitOption("-q kanga -C 0");

  Int_t nEvtPerSample = pdfOnResDK.totalNumEvts();
  TString pathPrefix = toyDir+TString("/fitMC");
  toyMC->createScripts(toyFile,pathPrefix,batchSetup,fitsPerJob,nSamples,nEvtPerSample);
}
  


void plotFullToy(const char* toyFile, Bool_t plotCPonly = kFALSE,
                 TString filePrefix = "")
{

  const Double_t pullMin = -4;
  const Double_t pullMax = 4;
  const Int_t pullBins = 40;

  assignTitles();

  RooArgList vars;
  if (!plotCPonly) {
    vars.add(*pdfOnResDK.totBBNumEvts());
    vars.add(*pdfOnResDK.DpiGoodD0NumEvts());
    vars.add(*pdfOnResDK.typeAsym(BdkEvtTypes::SIG_GOOD_D));
    vars.add(*pdfOnResDK.qqBadD0NumEvts());
    vars.add(*pdfOnResDK.sigGoodD0NumEvts());
    vars.add(*pdfOnResDK.DPiXFrac());
  }

  RooArgList xyset;
  xyset.add(*dalitzHolderP.sigGoodD0Type().x());
  xyset.add(*dalitzHolderP.sigGoodD0Type().y());
  xyset.add(*dalitzHolderN.sigGoodD0Type().x());
  xyset.add(*dalitzHolderN.sigGoodD0Type().y());

  RooArgList polarSet;
  polarSet.add(*dalitzHolderP.sigGoodD0Type().rho());
  polarSet.add(*dalitzHolderP.sigGoodD0Type().theta());
  polarSet.add(*dalitzHolderN.sigGoodD0Type().rho());
  polarSet.add(*dalitzHolderN.sigGoodD0Type().theta());

  vars.add(pdfOnResDK.cpParams());

  //  if (plotxy) vars.add(xyset);
  //  else vars.add(polarSet);
                

  ostringstream os;

  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.06,"x");
  gStyle->SetLabelSize(0.06,"y");
  gStyle->SetNdivisions(508,"x");
  gStyle->SetNdivisions(506,"y");

  for (int i=0; i<vars.getSize(); i++) {
    TFile* f = new TFile(toyFile);
    RooDataSet* data = (RooDataSet*)f->Get("fitParData");
    gROOT->cd();

    TCanvas* can = 0;

    Bool_t fitPull = kTRUE;
    Bool_t fitErr = kTRUE;
    Bool_t fitVar = kTRUE;
    Double_t varMin = -1;
    Double_t varMax = -2;
    Int_t varBins = 40;
    Double_t errMin = -1;
    Double_t errMax = -2;
    Int_t errBins = 40;
    // Separate settings for theses
    if ((&vars[i]==dalitzHolderN.sigGoodD0Type().theta()) ||
        (&vars[i]==dalitzHolderP.sigGoodD0Type().theta())) {
      fitErr = kFALSE;
      /*      varMin = -1.5;
      varMax = 1.5;
      varBins = 30;*/
      errMin = 0;
      errMax = 60;
      errBins = 40;
    }
    /*
    else if (&vars[i]==pdfOnResDK.totBBNumEvts()) {
      errMin = 60;
      errMax = 72;
    }
    else if (&vars[i]==pdfOnResDK.qqBadD0NumEvts()) {
      errMin = 65;
      errMax = 72;
    }
    */
    if (data->get()->contains(vars[i])) {
      can = plotToy2(data, (RooRealVar&)vars[i], fitPull, fitErr, fitVar,
                     errMin, errMax, errBins, pullMin, pullMax, pullBins,
                     varMin, varMax, varBins);
      
      os << vars[i].GetTitle();
      TString line1, line2;
      
      BdkPdfGauss* pdfs[] = {toyGaussPull,toyGaussErr,toyGaussVar};
      for (int j=0; j<3; j++) {
        if (pdfs[j]) {
          TString str1, str2;
          str1.Form("$%.3f\\pm %.3f$",pdfs[j]->b()->getVal(),pdfs[j]->b()->getError());
          str2.Form("$%.3f\\pm %.3f$",pdfs[j]->s()->getVal(),pdfs[j]->s()->getError());
          line1 += TString(" & ")+str1;
          line2 += TString(" & ")+str2;
        }
        else {
          line1 += TString(" & -");
          line2 += TString(" & -");
        }
      }
      os << line1 << endl<< line2<<endl;
      
      f->Close();
      delete f;
      TString file = filePrefix+"fullToy-";
      file += vars[i].GetName();
      can->SaveAs(file+".root");
      can->SaveAs(file+".eps");  
    }
  }
  cout << endl << "Pull fits:"<<endl;
  cout << os.str();
  TString file = filePrefix+"fullToy-fits.txt";
  ofstream of(file);
  of << os.str();
  of.close();
  cout << "... saved to "<<file<<endl;
}

void plotFullToyNLL(const char* toyFile, Bool_t showDataNLL = kFALSE)
{
  const Double_t BRdataNLL = -20505.9;
  const Double_t CPdataNLL = -19068.5;

  RooRealVar nll("NLL","BR-fit NLL",0);
  RooRealVar xynll("xyNLL","CP-fit NLL",0);

  TFile f(toyFile);  
  RooDataSet* data = (RooDataSet*)f.Get("fitParData");
  RooDataSet* red = data->reduce(RooArgSet(nll,xynll),"xyNLL<0 && NLL<0");
  gROOT->cd();

  RooPlot* p1 = nll.frame(AutoRange(*red),Bins(20));
  RooPlot* p2 = xynll.frame(AutoRange(*red),Bins(20));
  data->plotOn(p1,MarkerSize(0.6),XErrorSize(0));
  data->plotOn(p2,MarkerSize(0.6),XErrorSize(0));

  p1->SetTitle("");
  p2->SetTitle("");
  p1->GetYaxis()->SetTitle("");
  p2->GetYaxis()->SetTitle("");
  p1->GetXaxis()->SetNdivisions(504);
  p2->GetXaxis()->SetNdivisions(504);

  TCanvas* can = new TCanvas("can","NLL",800,400);
  can->Divide(2,1);
  can->cd(1);
  p1->Draw();
  TArrow* a1 = new TArrow(BRdataNLL,p1->GetMaximum(),
                          BRdataNLL,p1->GetMinimum(),0.03,"-|>");
  a1->SetLineStyle(kDashed);
  a1->SetLineColor(kBlue);
  a1->SetFillColor(kBlue);
  a1->Draw();

  can->cd(2);
  p2->Draw();
  TArrow* a2 = new TArrow(CPdataNLL,p2->GetMaximum(),
                          CPdataNLL,p2->GetMinimum(),0.03,"-|>");
  a2->SetLineStyle(kDashed);
  a2->SetLineColor(kBlue);
  a2->SetFillColor(kBlue);
  a2->Draw();

  can->SaveAs("fullToy-NLL.root");
  can->SaveAs("fullToy-NLL.eps");
  delete red;
}


void assignTitles() {
  pdfOnResDK.totBBNumEvts()->SetTitle("N(totBB)");
  pdfOnResDK.DpiGoodD0NumEvts()->SetTitle("N(D#pi_{D})");
  pdfOnResDK.typeAsym(BdkEvtTypes::SIG_GOOD_D)->SetTitle("A(DK_{D})");
  pdfOnResDK.qqBadD0NumEvts()->SetTitle("N(qqBadD)");
  pdfOnResDK.sigGoodD0NumEvts()->SetTitle("N(DK_{D})");
  pdfOnResDK.DPiXFrac()->SetTitle("R(D#pi X)");
  if (dalitzHolderN.sigGoodD0Type().coord()==BdkPdfDKDalitz::CART) {
    dalitzHolderP.sigGoodD0Type().x()->SetTitle("x+");
    dalitzHolderP.sigGoodD0Type().y()->SetTitle("y+");
    dalitzHolderN.sigGoodD0Type().x()->SetTitle("x-");
    dalitzHolderN.sigGoodD0Type().y()->SetTitle("y-");
  }
  else {
    dalitzHolderP.sigGoodD0Type().rho()->SetTitle("#rho+");
    dalitzHolderP.sigGoodD0Type().theta()->SetTitle("#theta+");
    dalitzHolderN.sigGoodD0Type().rho()->SetTitle("#rho-");
    dalitzHolderN.sigGoodD0Type().theta()->SetTitle("#theta-");
  }
}

